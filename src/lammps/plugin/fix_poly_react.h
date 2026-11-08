// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2026 Zidan Zhang <zhangzidan@gmail.com>
//
// fix poly/react: reactions between coarse-grained beads according to a reaction file (see
// ../README.md). The topology updates (bond storage, special lists, angles) follow fix bond/create
// of LAMMPS (MC package, GPL-2.0).
//
//   fix ID group-ID poly/react file seed keyword value ...
//
//   every N               steps between checks (default: interval of the file)
//   partner all|random    all: every pair within the cutoff is tried once (default)
//                         random: every bead that can take part in a pair reaction tries one
//                         random neighbor within the cutoff, so a pair of such beads has two
//                         chances per check, as in the lattice models
//   functionality name    integer per-atom vector (i_name) with the functionalities (default: 0)
//   max N                 most reactions per check (default: max_per_interval of the file)
//   log file              write every reaction to a file

#ifdef FIX_CLASS
// clang-format off
FixStyle(poly/react,FixPolyReact);
// clang-format on
#else

#ifndef LMP_FIX_POLY_REACT_H
#define LMP_FIX_POLY_REACT_H

#include "fix.h"

#include <array>
#include <cstdio>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace polyreact {
class Engine;
struct ReactionFile;
}  // namespace polyreact

namespace LAMMPS_NS {

class FixPolyReact : public Fix {
 public:
  FixPolyReact(class LAMMPS *, int, char **);
  ~FixPolyReact() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void post_integrate() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double compute_vector(int) override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  double memory_usage() override;

  /// A reaction found by one processor; all processors see all of them before any is carried out.
  struct Proposal {
    double key;  ///< priority in the conflict resolution: random, or the distance with "nearest"
    double distance;
    tagint center;
    tagint partner;  ///< 0 for a unimolecular reaction
    int reaction;    ///< index into the reaction file
    int bond;        ///< 1 if a bond forms
    int type[2];     ///< new LAMMPS types of center and partner
    int functionality[2];
  };

 private:
  enum { COMM_STATE = 1, COMM_SPECIAL = 2 };

  std::unique_ptr<polyreact::ReactionFile> file_;
  std::unique_ptr<polyreact::Engine> engine_;
  class RanMars *random_;
  class NeighList *list_;
  bool random_partner_;
  bool nearest_;
  std::string functionality_name_;
  int functionality_index_;
  bigint max_per_check_;
  double cutoff_, cutsq_;
  int commflag_;
  int special_depth_;  ///< special lists hold 1-2 (1), 1-2 and 1-3 (2) or up to 1-4 neighbors (3)
  bool angles_;
  FILE *log_;

  std::vector<int> bond_type_;   ///< per reaction: bond type of the new bond (0: none)
  std::vector<int> angle_type_;  ///< per reaction: angle type of the angles it closes (0: none)

  /// Molecules as a union-find over atom IDs and the bonds of every atom, the same on every
  /// processor.
  std::vector<tagint> parent_;
  std::vector<int> degree_;

  // statistics
  enum { MONOMER_MONOMER, POLYMER_MONOMER, POLYMER_POLYMER, INTRAMOLECULAR };
  static constexpr int NHEAD = 12;  ///< elements of the global vector before the counts per reaction
  bigint total_, last_;
  std::array<bigint, 4> tally_;  ///< new bonds by what they connect; a bead without bonds is a monomer
  std::vector<bigint> counts_;   ///< per reaction
  bigint stats_step_, stats_bonds_;
  std::array<double, 5> stats_;  ///< molecules of 2 or more atoms: number, Mn, Mw, largest, Mw without the largest

  // bonds created in the current check: (smaller ID, larger ID) -> angle type
  std::map<std::pair<tagint, tagint>, int> created_;
  std::vector<tagint> copy_;
  int nangles_, overflow_;

  [[nodiscard]] int functionality(int i) const;
  [[nodiscard]] bool bonded(int i, tagint other) const;
  [[nodiscard]] bool same_residue(int i, int j) const;
  tagint root(tagint t);
  void unite(tagint a, tagint b);
  void build_molecules();
  [[nodiscard]] bigint bonds_created() const;
  void molecule_statistics();
  void propose(int i, int j, double rsq, std::vector<Proposal> &out);
  void apply(const Proposal &p);
  void add_bond(int i, tagint other, int type);
  void update_topology();
  void rebuild_special_one(int m);
  void create_angles(int m);
  [[nodiscard]] int created_angle_type(tagint a, tagint b) const;
};

}  // namespace LAMMPS_NS

#endif
#endif
