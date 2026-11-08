// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2026 Zidan Zhang <zhangzidan@gmail.com>
//
// fix poly/react (see fix_poly_react.h). Every check runs in four steps:
//
//   1. each processor tries its own beads: unimolecular reactions, and pairs with neighbors
//      within the cutoff, decided by the rule engine; the reactions that pass are proposals;
//   2. the proposals of all processors are gathered on every processor;
//   3. every processor goes through them in the same order (random, or by distance with
//      "nearest") and accepts one if neither bead has reacted in this check yet, so each bead
//      reacts at most once per check; a bond within a molecule is refused if the reaction is not
//      intramolecular, also when earlier bonds of the same check joined the two molecules;
//   4. each processor carries out the accepted reactions on the beads it owns: types,
//      functionalities, the bond, the special lists and the angles that the bond closes.
//
// Molecules are followed with a union-find over atom IDs that every processor keeps in full.
// rebuild_special_one(), create_angles() and dedup() are adapted from fix bond/create (LAMMPS,
// MC package, GPL-2.0, Sandia Corporation).

#include "fix_poly_react.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "label_map.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "random_mars.h"
#include "update.h"

#include "polyreact/engine.hpp"
#include "polyreact/reaction_file.hpp"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <tuple>
#include <unordered_set>

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {

// remove duplicates of copy[nstart, nstop) among copy[0, nstop); from fix bond/create
int dedup(int nstart, int nstop, tagint *copy)
{
  int m = nstart;
  while (m < nstop) {
    int i;
    for (i = 0; i < m; i++)
      if (copy[i] == copy[m]) {
        copy[m] = copy[nstop - 1];
        nstop--;
        break;
      }
    if (i == m) m++;
  }
  return nstop;
}

}    // namespace

/* ---------------------------------------------------------------------- */

FixPolyReact::FixPolyReact(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), random_(nullptr), list_(nullptr), random_partner_(false), nearest_(false),
    functionality_index_(-1), max_per_check_(-1), cutoff_(0.0), cutsq_(0.0), commflag_(0),
    special_depth_(3), angles_(false), log_(nullptr), total_(0), last_(0), tally_{},
    stats_step_(-1), stats_bonds_(-1), stats_{}, nangles_(0), overflow_(0)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix poly/react", error);
  if (atom->molecular != Atom::MOLECULAR || !atom->avec->bonds_allow)
    error->all(FLERR, "Fix poly/react requires a molecular atom style with bonds");

  // the reaction file, read by processor 0 and parsed by all

  std::string text;
  int size = 0;
  if (comm->me == 0) {
    std::ifstream in(arg[3], std::ios::binary);
    if (in.is_open()) {
      text.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
      size = static_cast<int>(text.size());
    } else
      size = -1;
  }
  MPI_Bcast(&size, 1, MPI_INT, 0, world);
  if (size < 0) error->all(FLERR, 3, "Fix poly/react cannot open the reaction file {}", arg[3]);
  text.resize(size);
  MPI_Bcast(text.data(), size, MPI_CHAR, 0, world);
  try {
    std::istringstream in(text);
    file_ = std::make_unique<polyreact::ReactionFile>(polyreact::parse_reaction_file(in, arg[3]));
  } catch (std::exception &e) {
    error->all(FLERR, 3, "Fix poly/react: {}", e.what());
  }

  const int seed = utils::inumeric(FLERR, arg[4], false, lmp);
  if (seed <= 0) error->all(FLERR, 4, "Fix poly/react seed must be positive");

  // keywords; the defaults come from the [general] section

  nevery = static_cast<int>(file_->general.interval);
  max_per_check_ = file_->general.max_per_interval;
  nearest_ = file_->general.nearest;
  std::string logname;
  for (int iarg = 5; iarg < narg; iarg += 2) {
    const std::string key = arg[iarg];
    if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix poly/react " + key, error);
    const std::string value = arg[iarg + 1];
    if (key == "every") {
      nevery = utils::inumeric(FLERR, value, false, lmp);
    } else if (key == "partner") {
      if (value == "all")
        random_partner_ = false;
      else if (value == "random")
        random_partner_ = true;
      else
        error->all(FLERR, iarg + 1, "Fix poly/react partner must be all or random");
    } else if (key == "functionality") {
      functionality_name_ = utils::strmatch(value, "^i_") ? value.substr(2) : value;
    } else if (key == "max") {
      max_per_check_ = utils::bnumeric(FLERR, value, false, lmp);
      if (max_per_check_ < 0) error->all(FLERR, iarg + 1, "Fix poly/react max must be >= 0");
    } else if (key == "log") {
      logname = value;
    } else
      error->all(FLERR, iarg, "Unknown fix poly/react keyword {}", key);
  }
  if (nevery <= 0) error->all(FLERR, "Fix poly/react needs a positive interval between checks");

  // bead types: atom type labels, or atom types given as numbers

  std::map<std::string, int, std::less<>> types;
  for (const std::string &name : file_->type_names()) {
    int t = -1;
    if (atom->labelmapflag) t = atom->lmap->find_type(name, Atom::ATOM);
    if (t < 0 && utils::is_integer(name)) t = utils::inumeric(FLERR, name, false, lmp);
    if (t < 1 || t > atom->ntypes)
      error->all(FLERR, 3, "Fix poly/react: bead type {} is neither an atom type label nor an atom type",
                 name);
    types.emplace(name, t - 1);
  }
  try {
    engine_ = std::make_unique<polyreact::Engine>(*file_, types, atom->ntypes);
  } catch (std::exception &e) {
    error->all(FLERR, 3, "Fix poly/react: {}", e.what());
  }

  const auto nreactions = file_->reactions.size();
  bond_type_.assign(nreactions, 0);
  angle_type_.assign(nreactions, 0);
  for (std::size_t r = 0; r < nreactions; ++r) {
    const polyreact::Reaction &reaction = file_->reactions[r];
    if (!reaction.active) continue;
    const polyreact::Kind kind = reaction.equation.kind;
    if (kind == polyreact::Kind::Dissociation)
      error->all(FLERR, 3, "Fix poly/react: [reaction_{}]: dissociation is not supported yet",
                 reaction.name);
    if (kind != polyreact::Kind::Association || reaction.virtual_bond) continue;
    const polyreact::Group *group = file_->find_group(reaction.group);
    if (!group || !group->bond_type)
      error->all(FLERR, 3, "Fix poly/react: [reaction_{}] creates bonds, but its group has no bond_type",
                 reaction.name);
    bond_type_[r] = *group->bond_type;
    if (bond_type_[r] < 1 || bond_type_[r] > atom->nbondtypes)
      error->all(FLERR, 3, "Fix poly/react: [reaction_{}]: bond type {} does not exist", reaction.name,
                 bond_type_[r]);
    angle_type_[r] = group->angle_type.value_or(0);
    if (angle_type_[r] != 0 && !atom->avec->angles_allow)
      error->all(FLERR, 3, "Fix poly/react: [reaction_{}] creates angles, but the atom style has none",
                 reaction.name);
    if (angle_type_[r] < 0 || angle_type_[r] > atom->nangletypes)
      error->all(FLERR, 3, "Fix poly/react: [reaction_{}]: angle type {} does not exist",
                 reaction.name, angle_type_[r]);
  }

  // functionalities are needed if any reaction looks at them or changes them

  bool windows = false;
  for (const polyreact::Rule &rule : engine_->rules()) {
    const int slots = rule.kind == polyreact::Kind::Unimolecular ? 1 : 2;
    for (int s = 0; s < slots; ++s)
      if (rule.min[s] != INT_MIN || rule.max[s] != INT_MAX || rule.delta[s] != 0) windows = true;
  }
  if (windows && functionality_name_.empty())
    error->all(FLERR, "Fix poly/react: the reactions use functionalities; give the integer per-atom "
               "vector that holds them with the keyword functionality");

  cutoff_ = engine_->max_cutoff();
  cutsq_ = cutoff_ * cutoff_;

  force_reneighbor = 1;
  next_reneighbor = -1;
  vector_flag = 1;
  size_vector = NHEAD + static_cast<int>(nreactions);
  global_freq = 1;
  extvector = 0;
  restart_global = 1;
  comm_forward = std::max(2, 1 + atom->maxspecial);
  counts_.assign(nreactions, 0);

  random_ = new RanMars(lmp, seed + comm->me);

  if (comm->me == 0) {
    for (const std::string &note : file_->notes) utils::logmesg(lmp, "fix poly/react: {}\n", note);
    if (!logname.empty()) {
      log_ = fopen(logname.c_str(), "w");
      if (!log_) error->one(FLERR, "Fix poly/react cannot open the log file {}", logname);
      utils::print(log_, "# step reaction center partner distance intramolecular\n");
    }
  }
}

/* ---------------------------------------------------------------------- */

FixPolyReact::~FixPolyReact()
{
  delete random_;
  if (log_) fclose(log_);
}

/* ---------------------------------------------------------------------- */

int FixPolyReact::setmask()
{
  return POST_INTEGRATE;
}

/* ---------------------------------------------------------------------- */

void FixPolyReact::init()
{
  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, Error::NOLASTLINE, "Fix poly/react does not support run style respa");
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, Error::NOLASTLINE, "Fix poly/react requires an atom map");

  if (!functionality_name_.empty()) {
    int flag = -1, cols = -1;
    functionality_index_ = atom->find_custom(functionality_name_.c_str(), flag, cols);
    if (functionality_index_ < 0 || flag != 0 || cols != 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix poly/react needs the integer per-atom vector i_{} for the functionalities",
                 functionality_name_);
  }

  // the special lists hold what Special::build() would hold for these special_bonds weights

  if (force->special_onefive)
    error->all(FLERR, Error::NOLASTLINE, "Fix poly/react does not support special_bonds one/five");
  const double *lj = force->special_lj;
  const double *coul = force->special_coul;
  if (lj[2] == 1.0 && coul[2] == 1.0 && lj[3] == 1.0 && coul[3] == 1.0)
    special_depth_ = 1;
  else if (lj[3] == 1.0 && coul[3] == 1.0)
    special_depth_ = 2;
  else
    special_depth_ = 3;
  if ((force->special_angle || force->special_dihedral) && comm->me == 0)
    error->warning(FLERR, "Fix poly/react keeps all 1-3 and 1-4 neighbors of new bonds, as if "
                   "special_bonds angle no dihedral no");

  angles_ = std::any_of(angle_type_.begin(), angle_type_.end(), [](int t) { return t > 0; });
  if (angles_ && !force->angle) {
    if (comm->me == 0)
      error->warning(FLERR, "Fix poly/react creates no angles, since no angle style is defined");
    angles_ = false;
  }

  for (const auto &fix : modify->get_fix_list()) {
    if (fix == this) continue;
    if (utils::strmatch(fix->style, "^bond/") || utils::strmatch(fix->style, "^poly/react")) {
      if (comm->me == 0)
        error->warning(FLERR, "Fix poly/react does not see bonds that fix {} changes; molecules "
                       "and intramolecular reactions may be wrong", fix->style);
      break;
    }
  }

  const int maxspecial = atom->maxspecial;
  comm_forward = std::max(2, 1 + maxspecial);
  copy_.resize(static_cast<std::size_t>(maxspecial) * maxspecial + maxspecial);

  // an occasional neighbor list with the largest reaction cutoff

  if (cutoff_ > 0.0) {
    const double skin = neighbor->skin;
    const double cutghost =
        force->pair ? std::max(force->pair->cutforce + skin, comm->cutghostuser) : comm->cutghostuser;
    if (cutoff_ + skin > cutghost)
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix poly/react cutoff plus skin {} exceeds the ghost atom range {}; increase it "
                 "with comm_modify cutoff", cutoff_ + skin, cutghost);
    if (neighbor->style == Neighbor::MULTI)
      error->all(FLERR, Error::NOLASTLINE, "Fix poly/react requires neighbor style bin or nsq");
    // random partners need all neighbors of an atom; otherwise a half list with newton on holds
    // every pair on one processor only
    const int flags = random_partner_ ? NeighConst::REQ_FULL : NeighConst::REQ_NEWTON_ON;
    auto *req = neighbor->add_request(this, flags | NeighConst::REQ_OCCASIONAL);
    req->set_cutoff_fixed(cutoff_ + skin);
  }
}

/* ---------------------------------------------------------------------- */

void FixPolyReact::init_list(int /*id*/, NeighList *ptr)
{
  list_ = ptr;
}

/* ---------------------------------------------------------------------- */

void FixPolyReact::setup(int /*vflag*/)
{
  build_molecules();
}

/* ----------------------------------------------------------------------
   union-find of all atoms over the current bonds, the same on every processor
------------------------------------------------------------------------- */

void FixPolyReact::build_molecules()
{
  const int nlocal = atom->nlocal;
  const tagint *tag = atom->tag;
  tagint maxtag = 0;
  for (int i = 0; i < nlocal; ++i) maxtag = std::max(maxtag, tag[i]);
  tagint maxall = 0;
  MPI_Allreduce(&maxtag, &maxall, 1, MPI_LMP_TAGINT, MPI_MAX, world);
  parent_.resize(static_cast<std::size_t>(maxall) + 1);
  std::iota(parent_.begin(), parent_.end(), tagint(0));
  degree_.assign(parent_.size(), 0);

  // every bond once: with newton_bond off, both atoms store it
  std::vector<tagint> mine;
  for (int i = 0; i < nlocal; ++i)
    for (int k = 0; k < atom->num_bond[i]; ++k) {
      if (atom->bond_type[i][k] <= 0) continue;    // turned off
      if (!force->newton_bond && tag[i] > atom->bond_atom[i][k]) continue;
      mine.push_back(tag[i]);
      mine.push_back(atom->bond_atom[i][k]);
    }

  const int nprocs = comm->nprocs;
  const int n = static_cast<int>(mine.size());
  std::vector<int> counts(nprocs), displs(nprocs);
  MPI_Allgather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, world);
  bigint total = 0;
  for (int p = 0; p < nprocs; ++p) {
    displs[p] = static_cast<int>(total);
    total += counts[p];
  }
  if (total > INT_MAX) error->all(FLERR, Error::NOLASTLINE, "Fix poly/react: too many bonds");
  std::vector<tagint> bonds(static_cast<std::size_t>(total));
  MPI_Allgatherv(mine.data(), n, MPI_LMP_TAGINT, bonds.data(), counts.data(), displs.data(),
                 MPI_LMP_TAGINT, world);
  for (std::size_t k = 0; k + 1 < bonds.size(); k += 2) {
    unite(bonds[k], bonds[k + 1]);
    ++degree_[bonds[k]];
    ++degree_[bonds[k + 1]];
  }
  stats_step_ = -1;
}

/* ---------------------------------------------------------------------- */

tagint FixPolyReact::root(tagint t)
{
  if (t >= static_cast<tagint>(parent_.size())) {    // atoms added since the setup
    const auto old = parent_.size();
    parent_.resize(static_cast<std::size_t>(t) + 1);
    std::iota(parent_.begin() + static_cast<std::ptrdiff_t>(old), parent_.end(),
              static_cast<tagint>(old));
    degree_.resize(parent_.size(), 0);
  }
  while (parent_[t] != t) {
    parent_[t] = parent_[parent_[t]];
    t = parent_[t];
  }
  return t;
}

/* ---------------------------------------------------------------------- */

void FixPolyReact::unite(tagint a, tagint b)
{
  a = root(a);
  b = root(b);
  if (a != b) parent_[std::max(a, b)] = std::min(a, b);
}

/* ---------------------------------------------------------------------- */

void FixPolyReact::post_integrate()
{
  if (update->ntimestep % nevery) return;

  // current positions, types and functionalities of the ghost atoms

  comm->forward_comm();
  commflag_ = COMM_STATE;
  comm->forward_comm(this, 2);

  const int nlocal = atom->nlocal;
  const int *mask = atom->mask;
  const int *type = atom->type;
  const tagint *tag = atom->tag;
  double **x = atom->x;

  // 1. proposals of my atoms

  std::vector<Proposal> mine;
  if (engine_->has_unimolecular()) {
    for (int i = 0; i < nlocal; ++i) {
      if (!(mask[i] & groupbit)) continue;
      const auto out = engine_->single({type[i] - 1, functionality(i)}, *random_);
      if (!out) continue;
      Proposal p{};
      p.key = nearest_ ? 0.0 : random_->uniform();
      p.center = tag[i];
      p.reaction = static_cast<int>(out->reaction);
      p.type[0] = out->center.type + 1;
      p.functionality[0] = out->center.functionality;
      mine.push_back(p);
    }
  }

  if (cutoff_ > 0.0) {
    neighbor->build_one(list_);
    const int inum = list_->inum;
    const int *ilist = list_->ilist;
    const int *numneigh = list_->numneigh;
    int **firstneigh = list_->firstneigh;

    for (int ii = 0; ii < inum; ++ii) {
      const int i = ilist[ii];
      if (!(mask[i] & groupbit)) continue;
      const int ti = type[i] - 1;
      if (!engine_->reactive(ti)) continue;
      const int *jlist = firstneigh[i];
      const int jnum = numneigh[i];

      if (random_partner_) {
        // one neighbor in the group within the cutoff, chosen at random (reservoir sampling)
        int chosen = -1, seen = 0;
        double chosen_rsq = 0.0;
        for (int jj = 0; jj < jnum; ++jj) {
          const int j = jlist[jj] & NEIGHMASK;
          if (!(mask[j] & groupbit) || tag[j] == tag[i]) continue;
          const double dx = x[i][0] - x[j][0], dy = x[i][1] - x[j][1], dz = x[i][2] - x[j][2];
          const double rsq = dx * dx + dy * dy + dz * dz;
          if (rsq >= cutsq_) continue;
          ++seen;
          if (random_->uniform() * seen < 1.0) {
            chosen = j;
            chosen_rsq = rsq;
          }
        }
        if (chosen >= 0) propose(i, chosen, chosen_rsq, mine);
        continue;
      }

      // every pair once (half list, newton on), with the closest image of the other atom
      for (int jj = 0; jj < jnum; ++jj) {
        const int j = jlist[jj] & NEIGHMASK;
        if (tag[i] == tag[j] || !(mask[j] & groupbit) || !engine_->may_react(ti, type[j] - 1)) continue;
        const double dx = x[i][0] - x[j][0], dy = x[i][1] - x[j][1], dz = x[i][2] - x[j][2];
        const double rsq = dx * dx + dy * dy + dz * dz;
        if (rsq >= cutsq_) continue;
        if (domain->closest_image(i, atom->map(tag[j])) != j) continue;
        propose(i, j, rsq, mine);
      }
    }
  }

  // 2. all proposals on every processor

  const int nprocs = comm->nprocs;
  const int n = static_cast<int>(mine.size() * sizeof(Proposal));
  std::vector<int> counts(nprocs), displs(nprocs);
  MPI_Allgather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, world);
  bigint bytes = 0;
  for (int p = 0; p < nprocs; ++p) {
    displs[p] = static_cast<int>(bytes);
    bytes += counts[p];
  }
  if (bytes > INT_MAX) error->all(FLERR, Error::NOLASTLINE, "Fix poly/react: too many reactions in one check");
  std::vector<Proposal> all(static_cast<std::size_t>(bytes) / sizeof(Proposal));
  MPI_Allgatherv(mine.data(), n, MPI_BYTE, all.data(), counts.data(), displs.data(), MPI_BYTE, world);

  // 3. the same order and the same choices on every processor

  std::sort(all.begin(), all.end(), [](const Proposal &a, const Proposal &b) {
    return std::tie(a.key, a.center, a.partner, a.reaction) <
        std::tie(b.key, b.center, b.partner, b.reaction);
  });

  std::unordered_set<tagint> reacted;
  created_.clear();
  last_ = 0;
  bigint nbonds = 0;
  for (const Proposal &p : all) {
    if (max_per_check_ >= 0 && last_ >= max_per_check_) break;
    if (reacted.count(p.center) || (p.partner && reacted.count(p.partner))) continue;
    bool intra = false;
    if (p.bond) {
      intra = root(p.center) == root(p.partner);
      if (intra && !file_->reactions[p.reaction].intramolecular) continue;
    }
    reacted.insert(p.center);
    if (p.partner) reacted.insert(p.partner);
    ++last_;
    ++counts_[p.reaction];
    if (p.bond) {
      // as in the lattice modules: a bead without bonds is a monomer
      const bool free_center = degree_[p.center] == 0, free_partner = degree_[p.partner] == 0;
      if (intra)
        ++tally_[INTRAMOLECULAR];
      else if (free_center && free_partner)
        ++tally_[MONOMER_MONOMER];
      else if (!free_center && !free_partner)
        ++tally_[POLYMER_POLYMER];
      else
        ++tally_[POLYMER_MONOMER];
      unite(p.center, p.partner);
      ++degree_[p.center];
      ++degree_[p.partner];
      ++nbonds;
      created_.emplace(std::minmax(p.center, p.partner), angle_type_[p.reaction]);
    }

    // 4. my atoms

    apply(p);
    if (log_)
      utils::print(log_, "{} {} {} {} {:.6g} {}\n", update->ntimestep,
                   file_->reactions[p.reaction].name, p.center, p.partner, p.distance, intra ? 1 : 0);
  }
  if (log_) fflush(log_);
  total_ += last_;
  atom->nbonds += nbonds;

  // ghost atoms get the new types at the reneighboring that follows right away

  if (last_ == 0) return;
  next_reneighbor = update->ntimestep;
  if (nbonds == 0) return;

  commflag_ = COMM_SPECIAL;
  comm->forward_comm(this);
  update_topology();
}

/* ----------------------------------------------------------------------
   try the pair of my atom I and atom J
------------------------------------------------------------------------- */

void FixPolyReact::propose(int i, int j, double rsq, std::vector<Proposal> &out)
{
  const tagint *tag = atom->tag;
  const int *type = atom->type;
  polyreact::PairContext context;
  context.distance = std::sqrt(rsq);
  context.bonded = bonded(i, tag[j]);
  context.same_molecule = root(tag[i]) == root(tag[j]);
  context.same_residue = same_residue(i, j);
  const auto o = engine_->pair({type[i] - 1, functionality(i)}, {type[j] - 1, functionality(j)},
                               context, *random_);
  if (!o) return;
  Proposal p{};
  p.key = nearest_ ? context.distance : random_->uniform();
  p.distance = context.distance;
  p.center = tag[i];
  p.partner = tag[j];
  p.reaction = static_cast<int>(o->reaction);
  p.bond = o->creates_bond ? 1 : 0;
  p.type[0] = o->center.type + 1;
  p.type[1] = o->partner.type + 1;
  p.functionality[0] = o->center.functionality;
  p.functionality[1] = o->partner.functionality;
  out.push_back(p);
}

/* ----------------------------------------------------------------------
   carry out an accepted reaction on the atoms I own
------------------------------------------------------------------------- */

void FixPolyReact::apply(const Proposal &p)
{
  const tagint ids[2] = {p.center, p.partner};
  const int sides = p.partner ? 2 : 1;
  for (int s = 0; s < sides; ++s) {
    const int i = atom->map(ids[s]);
    if (i < 0 || i >= atom->nlocal) continue;
    atom->type[i] = p.type[s];
    if (functionality_index_ >= 0) atom->ivector[functionality_index_][i] = p.functionality[s];
    if (p.bond) add_bond(i, ids[1 - s], bond_type_[p.reaction]);
  }
}

/* ----------------------------------------------------------------------
   new bond of my atom I with atom OTHER: bond storage and 1-2 special neighbor,
   as in fix bond/create; the other atom does the same, whatever processor owns it
------------------------------------------------------------------------- */

void FixPolyReact::add_bond(int i, tagint other, int btype)
{
  if (!force->newton_bond || atom->tag[i] < other) {
    int &num = atom->num_bond[i];
    if (num == atom->bond_per_atom)
      error->one(FLERR, Error::NOLASTLINE,
                 "New bond from fix poly/react exceeds the bonds per atom; use extra/bond/per/atom");
    atom->bond_type[i][num] = btype;
    atom->bond_atom[i][num] = other;
    num++;
  }

  // the bond partner becomes the first 1-2 neighbor; drop it from the 1-3 and 1-4 neighbors,
  // which rebuild_special_one() rebuilds

  tagint *slist = atom->special[i];
  int *ns = atom->nspecial[i];
  int n1 = ns[0], n2 = ns[1], n3 = ns[2];
  int m;
  for (m = n1; m < n3; m++)
    if (slist[m] == other) break;
  if (m < n3) {
    for (int k = m; k < n3 - 1; k++) slist[k] = slist[k + 1];
    n3--;
    if (m < n2) n2--;
  }
  if (n3 == atom->maxspecial)
    error->one(FLERR, Error::NOLASTLINE,
               "New bond from fix poly/react exceeds the special list size; use extra/special/per/atom");
  for (m = n3; m > n1; m--) slist[m] = slist[m - 1];
  slist[n1] = other;
  ns[0] = n1 + 1;
  ns[1] = n2 + 1;
  ns[2] = n3 + 1;
}

/* ----------------------------------------------------------------------
   rebuild the special lists of my atoms near a new bond and create the angles it closes
------------------------------------------------------------------------- */

void FixPolyReact::update_topology()
{
  std::unordered_set<tagint> ends;
  for (const auto &bond : created_) {
    ends.insert(bond.first.first);
    ends.insert(bond.first.second);
  }

  const tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  nangles_ = 0;
  overflow_ = 0;
  for (int i = 0; i < atom->nlocal; ++i) {
    bool influenced = ends.count(tag[i]) > 0;
    const int n = special_depth_ == 3 ? nspecial[i][1] : nspecial[i][0];
    for (int k = 0; k < n && !influenced; ++k) influenced = ends.count(special[i][k]) > 0;
    if (!influenced) continue;
    rebuild_special_one(i);
    if (angles_) create_angles(i);
  }

  int overflow = 0;
  MPI_Allreduce(&overflow_, &overflow, 1, MPI_INT, MPI_MAX, world);
  if (overflow)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix poly/react created too many angles per atom; use extra/angle/per/atom");
  if (angles_) {
    int all = 0;
    MPI_Allreduce(&nangles_, &all, 1, MPI_INT, MPI_SUM, world);
    if (!force->newton_bond) all /= 3;
    atom->nangles += all;
  }
}

/* ----------------------------------------------------------------------
   rebuild the 1-3 and 1-4 neighbors of atom M from the 1-2 neighbors of its 1-2 and 1-3
   neighbors, as far as the special_bonds weights need them; from fix bond/create
------------------------------------------------------------------------- */

void FixPolyReact::rebuild_special_one(int m)
{
  const tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  tagint *copy = copy_.data();

  const auto neighbors_of = [&](int from, int to, int count) {
    for (int i = from; i < to; i++) {
      const int n = atom->map(copy[i]);
      if (n < 0)
        error->one(FLERR, Error::NOLASTLINE, "Fix poly/react needs ghost atoms from further away");
      for (int j = 0; j < nspecial[n][0]; j++)
        if (special[n][j] != tag[m]) copy[count++] = special[n][j];
    }
    return count;
  };

  int cn1 = 0;
  for (int i = 0; i < nspecial[m][0]; i++) copy[cn1++] = special[m][i];
  int cn2 = cn1;
  if (special_depth_ >= 2) {
    cn2 = dedup(cn1, neighbors_of(0, cn1, cn1), copy);
    if (cn2 > atom->maxspecial)
      error->one(FLERR, Error::NOLASTLINE, "Special list size exceeded in fix poly/react");
  }
  int cn3 = cn2;
  if (special_depth_ >= 3) {
    cn3 = dedup(cn2, neighbors_of(cn1, cn2, cn2), copy);
    if (cn3 > atom->maxspecial)
      error->one(FLERR, Error::NOLASTLINE, "Special list size exceeded in fix poly/react");
  }

  nspecial[m][0] = cn1;
  nspecial[m][1] = cn2;
  nspecial[m][2] = cn3;
  memcpy(special[m], copy, cn3 * sizeof(tagint));
}

/* ----------------------------------------------------------------------
   angle type for a new angle with the bond A-B: that of the reaction that created the bond,
   0 if the bond is new but its group gives no angle type, -1 if the bond is not new
------------------------------------------------------------------------- */

int FixPolyReact::created_angle_type(tagint a, tagint b) const
{
  const auto it = created_.find(std::minmax(a, b));
  return it == created_.end() ? -1 : it->second;
}

/* ----------------------------------------------------------------------
   create the angles of atom M that contain a new bond; from fix bond/create
   newton_bond on: atom M is the central atom; off: atom M is any of the 3 atoms
------------------------------------------------------------------------- */

void FixPolyReact::create_angles(int m)
{
  const tagint *tag = atom->tag;
  int **nspecial = atom->nspecial;
  tagint **special = atom->special;
  int &num_angle = atom->num_angle[m];

  const auto add = [&](tagint i1, tagint i2, tagint i3) {
    const int t1 = created_angle_type(i1, i2), t2 = created_angle_type(i2, i3);
    const int t = t1 > 0 ? t1 : t2;
    if (t <= 0) return;
    if (num_angle == atom->angle_per_atom) {
      overflow_ = 1;
      return;
    }
    atom->angle_type[m][num_angle] = t;
    atom->angle_atom1[m][num_angle] = i1;
    atom->angle_atom2[m][num_angle] = i2;
    atom->angle_atom3[m][num_angle] = i3;
    num_angle++;
    nangles_++;
  };

  // atom M as the central atom: pairs of its 1-2 neighbors
  const tagint *s2list = special[m];
  const int n2 = nspecial[m][0];
  for (int i = 0; i < n2; i++)
    for (int j = i + 1; j < n2; j++) add(s2list[i], tag[m], s2list[j]);
  if (force->newton_bond) return;

  // atom M as an end atom: 1-2 neighbors of its 1-2 neighbors
  for (int i = 0; i < n2; i++) {
    const tagint i2 = s2list[i];
    const int i2local = atom->map(i2);
    if (i2local < 0)
      error->one(FLERR, Error::NOLASTLINE, "Fix poly/react needs ghost atoms from further away");
    for (int j = 0; j < nspecial[i2local][0]; j++) {
      const tagint i3 = special[i2local][j];
      if (i3 != tag[m]) add(tag[m], i2, i3);
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixPolyReact::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;
  if (commflag_ == COMM_STATE) {
    for (int i = 0; i < n; i++) {
      const int j = list[i];
      buf[m++] = ubuf(atom->type[j]).d;
      buf[m++] = ubuf(functionality(j)).d;
    }
    return m;
  }
  for (int i = 0; i < n; i++) {
    const int j = list[i];
    const int ns = atom->nspecial[j][0];
    buf[m++] = ubuf(ns).d;
    for (int k = 0; k < ns; k++) buf[m++] = ubuf(atom->special[j][k]).d;
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixPolyReact::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  const int last = first + n;
  if (commflag_ == COMM_STATE) {
    for (int i = first; i < last; i++) {
      atom->type[i] = static_cast<int>(ubuf(buf[m++]).i);
      const int f = static_cast<int>(ubuf(buf[m++]).i);
      if (functionality_index_ >= 0) atom->ivector[functionality_index_][i] = f;
    }
    return;
  }
  for (int i = first; i < last; i++) {
    const int ns = static_cast<int>(ubuf(buf[m++]).i);
    atom->nspecial[i][0] = ns;
    for (int k = 0; k < ns; k++) atom->special[i][k] = static_cast<tagint>(ubuf(buf[m++]).i);
  }
}

/* ---------------------------------------------------------------------- */

int FixPolyReact::functionality(int i) const
{
  return functionality_index_ < 0 ? 0 : atom->ivector[functionality_index_][i];
}

/* ---------------------------------------------------------------------- */

bool FixPolyReact::bonded(int i, tagint other) const
{
  const tagint *slist = atom->special[i];
  for (int k = 0; k < atom->nspecial[i][0]; k++)
    if (slist[k] == other) return true;
  return false;
}

/* ----------------------------------------------------------------------
   same residue: the same nonzero molecule ID, i.e. the same initial molecule
------------------------------------------------------------------------- */

bool FixPolyReact::same_residue(int i, int j) const
{
  return atom->molecule_flag && atom->molecule[i] != 0 && atom->molecule[i] == atom->molecule[j];
}

/* ---------------------------------------------------------------------- */

bigint FixPolyReact::bonds_created() const
{
  return tally_[0] + tally_[1] + tally_[2] + tally_[3];
}

/* ----------------------------------------------------------------------
   sizes of the molecules of 2 or more atoms, from the union-find; computed when asked for
   and kept until the step or the bonds change
------------------------------------------------------------------------- */

void FixPolyReact::molecule_statistics()
{
  if (stats_step_ == update->ntimestep && stats_bonds_ == bonds_created()) return;
  stats_step_ = update->ntimestep;
  stats_bonds_ = bonds_created();

  std::vector<tagint> size(parent_.size(), 0);
  for (tagint t = 1; t < static_cast<tagint>(parent_.size()); ++t) ++size[root(t)];
  double n = 0.0, s1 = 0.0, s2 = 0.0, largest = 0.0;
  for (const tagint s : size) {
    if (s < 2) continue;
    const auto l = static_cast<double>(s);
    n += 1.0;
    s1 += l;
    s2 += l * l;
    largest = std::max(largest, l);
  }
  stats_[0] = n;
  stats_[1] = n > 0.0 ? s1 / n : 0.0;
  stats_[2] = s1 > 0.0 ? s2 / s1 : 0.0;
  stats_[3] = largest;
  stats_[4] = s1 > largest ? (s2 - largest * largest) / (s1 - largest) : 0.0;
}

/* ----------------------------------------------------------------------
   1 reactions in all, 2 in the last check, 3 bonds created, of which 4 monomer-monomer,
   5 polymer-monomer, 6 polymer-polymer, 7 within a molecule (rings); molecules of 2 or more
   atoms: 8 number, 9 Mn, 10 Mw, 11 largest, 12 Mw without the largest (reduced degree of
   polymerization); 13... reactions per reaction section
------------------------------------------------------------------------- */

double FixPolyReact::compute_vector(int n)
{
  if (n == 0) return static_cast<double>(total_);
  if (n == 1) return static_cast<double>(last_);
  if (n == 2) return static_cast<double>(bonds_created());
  if (n < 7) return static_cast<double>(tally_[n - 3]);
  if (n < NHEAD) {
    molecule_statistics();
    return stats_[n - 7];
  }
  return static_cast<double>(counts_[n - NHEAD]);
}

/* ---------------------------------------------------------------------- */

void FixPolyReact::write_restart(FILE *fp)
{
  std::vector<double> list{static_cast<double>(total_)};
  for (const bigint t : tally_) list.push_back(static_cast<double>(t));
  list.push_back(static_cast<double>(counts_.size()));
  for (const bigint c : counts_) list.push_back(static_cast<double>(c));
  if (comm->me == 0) {
    const int size = static_cast<int>(list.size() * sizeof(double));
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list.data(), sizeof(double), list.size(), fp);
  }
}

/* ---------------------------------------------------------------------- */

void FixPolyReact::restart(char *buf)
{
  const auto *list = reinterpret_cast<const double *>(buf);
  total_ = static_cast<bigint>(list[0]);
  for (std::size_t k = 0; k < tally_.size(); ++k) tally_[k] = static_cast<bigint>(list[1 + k]);
  const auto n = static_cast<std::size_t>(list[5]);
  if (n == counts_.size()) {
    for (std::size_t r = 0; r < n; ++r) counts_[r] = static_cast<bigint>(list[6 + r]);
  } else if (comm->me == 0) {
    error->warning(FLERR, "Fix poly/react: the restart has {} reactions, the reaction file {}; "
                   "the counts per reaction start from zero", n, counts_.size());
  }
}

/* ---------------------------------------------------------------------- */

double FixPolyReact::memory_usage()
{
  return static_cast<double>((parent_.capacity() + copy_.capacity()) * sizeof(tagint) +
                             degree_.capacity() * sizeof(int));
}
