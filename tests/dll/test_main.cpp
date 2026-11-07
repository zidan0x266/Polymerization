// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Self-contained unit tests (no framework dependency). The most important ones compare every
// time step with a brute-force implementation of the dynamic lattice liquid rules, and the
// reaction stage with the exact per-step reaction probabilities of its two partner rules.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "dll/io/config.hpp"
#include "dll/io/writers.hpp"
#include "dll/model/initial.hpp"
#include "dll/model/system.hpp"
#include "dll/sim/crossing.hpp"
#include "dll/sim/motion.hpp"
#include "dll/sim/reaction.hpp"
#include "dll/sim/run.hpp"
#include "poly/chem/progress.hpp"
#include "poly/chem/scheme.hpp"

namespace {

int g_failures = 0;
int g_checks = 0;

#define CHECK(cond)                                                                                  \
  do {                                                                                               \
    ++g_checks;                                                                                      \
    if (!(cond)) {                                                                                   \
      ++g_failures;                                                                                  \
      std::fprintf(stderr, "CHECK failed: %s  (%s:%d)\n", #cond, __FILE__, __LINE__);                \
    }                                                                                                \
  } while (0)

using namespace dll;

DllSystem make_system(InitialConfiguration kind, int cells, std::size_t chain_length, std::size_t n_chains,
                      std::uint64_t seed) {
  poly::Rng rng(seed);
  return DllSystem(cells, chain_length, initial_chains(kind, cells, chain_length, n_chains, rng));
}

std::vector<std::uint8_t> random_field(std::size_t n, poly::Rng& rng) {
  std::vector<std::uint8_t> field(n);
  for (auto& d : field) d = static_cast<std::uint8_t>(rng.below(12));
  return field;
}

// ------------------------------------------------------------------------------------------------
// Crossing geometry: fixed contacts, invariants on random lattice configurations, and an
// independent check by sampling time and both segments on a grid.

Path reversed(Path p) { return {p.start + p.step, -p.step}; }

double sampled_min_distance2(Path a, Path b, Path c, Path d) {
  const auto at = [](Path p, double t, int k) { return (k == 0 ? p.start.x : k == 1 ? p.start.y : p.start.z) +
                                                       t * (k == 0 ? p.step.x : k == 1 ? p.step.y : p.step.z); };
  double best = 1e30;
  for (int i = 0; i <= 200; ++i) {
    const double t = i / 200.0;
    for (int j = 0; j <= 20; ++j) {
      for (int k = 0; k <= 20; ++k) {
        const double s = j / 20.0, u = k / 20.0;
        double d2 = 0;
        for (int m = 0; m < 3; ++m) {
          const double x = at(a, t, m) + s * (at(b, t, m) - at(a, t, m)) - at(c, t, m) - u * (at(d, t, m) - at(c, t, m));
          d2 += x * x;
        }
        best = std::min(best, d2);
      }
    }
  }
  return best;
}

void test_crossing_geometry() {
  // Two parallel bonds that pass through each other, and a bead passing through a bond.
  const Path a{{0, 0, 0}, {0, 1, 1}}, b{{-1, 0, 1}, {0, 1, 1}}, c{{-1, 1, 2}, {-1, -1, 0}}, d{{0, 1, 1}, {-1, -1, 0}};
  CHECK(segments_touch(a, b, c, d) && segments_touch(d, c, b, a));
  CHECK(segments_touch(reversed(a), reversed(b), reversed(c), reversed(d)));
  CHECK(sampled_min_distance2(a, b, c, d) < 1e-3);
  const Path pa{{0, 0, 0}, {-1, 0, -1}}, pb{{0, 1, -1}, {-1, 0, -1}}, bead{{-1, 0, -1}, {1, 1, 0}};
  CHECK(segments_touch(pa, pb, bead, bead) && segments_touch(pb, pa, bead, bead));
  CHECK(sampled_min_distance2(pa, pb, bead, bead) < 1e-3);
  CHECK(!segments_touch(pa, pb, {{2, 2, 0}, {1, 1, 0}}, {{2, 2, 0}, {1, 1, 0}}));

  // Random valid configurations: four distinct beads, no two on one site before or after the
  // step, no swaps, no bead entering the site of one that stays, bonds kept.
  poly::Rng rng(41);
  std::vector<Vec3i> steps{{0, 0, 0}};
  for (const Vec3i v : FccGeometry::kNeighbors) steps.push_back(v);
  std::vector<Vec3i> near;
  for (int x = -4; x <= 4; ++x) {
    for (int y = -4; y <= 4; ++y) {
      for (int z = -4; z <= 4; ++z) {
        if (((x + y + z) & 1) == 0 && x * x + y * y + z * z <= 16) near.push_back({x, y, z});
      }
    }
  }
  const auto pick = [&](const std::vector<Vec3i>& v) { return v[rng.below(static_cast<std::uint32_t>(v.size()))]; };
  int contacts = 0, invariant = 0, both_ends = 0, sampled = 0, agree = 0;
  for (int it = 0; it < 400000; ++it) {
    const Path p{{}, pick(steps)}, q{pick(steps), pick(steps)};
    const bool point = rng.below(2) == 0;
    const Path r{pick(near), pick(steps)};
    const Path s = point ? r : Path{r.start + pick(steps), pick(steps)};
    if (s.start == r.start && !point) continue;
    if ((q.start + q.step - p.step).norm2() != 2 || q.start.norm2() != 2) continue;
    if (!point && ((s.start - r.start).norm2() != 2 || (s.start + s.step - r.start - r.step).norm2() != 2)) continue;
    const std::vector<Path> all = point ? std::vector<Path>{p, q, r} : std::vector<Path>{p, q, r, s};
    bool valid = true;
    for (std::size_t i = 0; i < all.size(); ++i) {
      for (std::size_t j = i + 1; j < all.size(); ++j) {
        const Vec3i i0 = all[i].start, i1 = i0 + all[i].step, j0 = all[j].start, j1 = j0 + all[j].step;
        if (i0 == j0 || i1 == j1 || (i1 == j0 && j1 == i0) || (all[i].step == Vec3i{} && i0 == j1) ||
            (all[j].step == Vec3i{} && j0 == i1)) {
          valid = false;
        }
      }
    }
    if (!valid) continue;
    const bool touch = segments_touch(p, q, r, s);
    contacts += touch ? 1 : 0;
    const bool same = touch == segments_touch(q, p, r, s) && touch == segments_touch(reversed(p), reversed(q), reversed(r), reversed(s)) &&
                      (point || (touch == segments_touch(r, s, p, q) && touch == segments_touch(p, q, s, r)));
    invariant += same ? 0 : 1;
    const bool moving_bond = (p.step != Vec3i{} && q.step != Vec3i{}) || (!point && r.step != Vec3i{} && s.step != Vec3i{});
    both_ends += (touch && !moving_bond) ? 1 : 0;
    if ((touch || sampled < 300) && sampled < 400) {
      ++sampled;
      const double m = sampled_min_distance2(p, q, r, s);
      agree += (touch ? m < 1e-3 : m > 1e-3) ? 1 : 0;
    }
  }
  CHECK(contacts > 0 && invariant == 0 && both_ends == 0);
  CHECK(sampled > 0 && agree == sampled);
}

// ------------------------------------------------------------------------------------------------
// Brute force: a site is on a loop if following the attempted moves leads back to it after at
// least three steps. In rounds, every loop that would break a bond given the loops accepted at the
// start of the round is removed until none does; then every loop that would make two bonds touch
// (or a bead without bonds touch a bond); and so on until no loop fails. Every pair is tested, in
// every periodic image.
struct Outcome {
  std::vector<Bead> occupant;
  std::vector<Vec3i> unwrapped;
};

Outcome reference_step(const DllSystem& sys, const std::vector<std::uint8_t>& field, bool allow_crossing) {
  const FccGeometry& lat = sys.lattice();
  const auto n = static_cast<Site>(sys.size());
  std::vector<Site> target(n);
  for (Site s = 0; s < n; ++s) target[s] = lat.neighbor(lat.position(s), field[s]);

  std::vector<Site> loop_id(n, n);  // smallest site of the loop, or n
  for (Site s = 0; s < n; ++s) {
    Site u = target[s];
    Site smallest = s;
    std::size_t length = 1;
    while (u != s && length <= n) {
      smallest = std::min(smallest, u);
      u = target[u];
      ++length;
    }
    if (u == s && length >= 3) loop_id[s] = smallest;
  }
  std::set<Site> accepted;
  for (Site s = 0; s < n; ++s) {
    if (loop_id[s] != n) accepted.insert(loop_id[s]);
  }
  const auto moving = [&](Bead b) { return loop_id[sys.site(b)] != n && accepted.contains(loop_id[sys.site(b)]); };
  const auto step_of = [&](Bead b) { return moving(b) ? FccGeometry::kNeighbors[field[sys.site(b)]] : Vec3i{}; };
  const auto new_pos = [&](Bead b) { return sys.unwrapped(b) + step_of(b); };
  std::vector<std::pair<Bead, Bead>> bonds;
  std::vector<Bead> loners;
  for (Bead b = 0; b < n; ++b) {
    for (const Bead p : sys.bonds(b)) {
      if (b < p) bonds.emplace_back(b, p);
    }
    if (sys.bonds(b).empty()) loners.push_back(b);
  }
  // Every image of the segment c-d within reach of the segment a-b (relative to a).
  const auto touches = [&](Bead a, Bead b, Bead c, Bead d) {
    const Vec3i rb = lat.minimum_image(sys.unwrapped(b) - sys.unwrapped(a));
    const Vec3i rc = lat.minimum_image(sys.unwrapped(c) - sys.unwrapped(a));
    const Vec3i rd = rc + lat.minimum_image(sys.unwrapped(d) - sys.unwrapped(c));
    const int L = lat.extent();
    for (int i = -2; i <= 2; ++i) {
      for (int j = -2; j <= 2; ++j) {
        for (int k = -2; k <= 2; ++k) {
          const Vec3i shift{i * L, j * L, k * L};
          const Vec3i c0 = rc + shift;
          if (std::abs(c0.x) > 6 || std::abs(c0.y) > 6 || std::abs(c0.z) > 6) continue;
          if (segments_touch({{}, step_of(a)}, {rb, step_of(b)}, {c0, step_of(c)}, {rd + shift, step_of(d)})) return true;
        }
      }
    }
    return false;
  };
  const auto bond_breakers = [&] {
    std::set<Site> violators;
    for (Site s = 0; s < n; ++s) {
      if (loop_id[s] == n || !accepted.contains(loop_id[s])) continue;
      const Bead b = sys.occupant(s);
      for (const Bead p : sys.bonds(b)) {
        if (lat.minimum_image(new_pos(p) - new_pos(b)).norm2() != 2) violators.insert(loop_id[s]);
      }
    }
    return violators;
  };
  const auto crossers = [&] {
    std::set<Site> violators;
    const auto fail = [&](std::initializer_list<Bead> beads) {
      for (const Bead q : beads) {
        if (moving(q)) violators.insert(loop_id[sys.site(q)]);
      }
    };
    for (const auto [a, b] : bonds) {
      if (!moving(a) && !moving(b)) continue;
      for (const auto [c, d] : bonds) {
        if (c == a || c == b || d == a || d == b || (!moving(c) && !moving(d))) continue;
        if (touches(a, b, c, d)) fail({a, b, c, d});
      }
      for (const Bead w : loners) {
        if (touches(a, b, w, w)) fail({a, b, w});
      }
    }
    return violators;
  };
  // Loops that break a bond go first, round after round; then those that let chains cross.
  for (;;) {
    for (auto v = bond_breakers(); !v.empty(); v = bond_breakers()) {
      for (const Site loop : v) accepted.erase(loop);
    }
    if (allow_crossing) break;
    const auto v = crossers();
    if (v.empty()) break;
    for (const Site loop : v) accepted.erase(loop);
  }
  Outcome out{std::vector<Bead>(n), std::vector<Vec3i>(n)};
  for (Site s = 0; s < n; ++s) {
    const Bead b = sys.occupant(s);
    const bool moves = loop_id[s] != n && accepted.contains(loop_id[s]);
    out.occupant[moves ? target[s] : s] = b;
    out.unwrapped[b] = new_pos(b);
  }
  return out;
}

void test_against_reference() {
  struct Case {
    InitialConfiguration kind;
    int cells;
    std::size_t chain_length, n_chains;
    bool allow_crossing;
  };
  const Case cases[] = {{InitialConfiguration::RandomWalk, 2, 4, 3, false},  {InitialConfiguration::RandomWalk, 3, 6, 9, false},
                        {InitialConfiguration::Aligned, 2, 4, 8, false},     {InitialConfiguration::Aligned, 3, 3, 36, false},
                        {InitialConfiguration::Aligned, 3, 2, 30, false},    {InitialConfiguration::RandomWalk, 4, 8, 0, false},
                        {InitialConfiguration::RandomWalk, 4, 4, 30, false}, {InitialConfiguration::RandomWalk, 3, 6, 9, true},
                        {InitialConfiguration::Aligned, 3, 2, 30, true}};
  poly::Rng rng(31);
  std::size_t moved = 0, rejected = 0, blocked = 0;
  for (const Case& c : cases) {
    DllSystem sys = make_system(c.kind, c.cells, c.chain_length, c.n_chains, 5);
    Motion motion(sys, c.allow_crossing);
    bool same = true;
    for (int step = 0; step < 2000; ++step) {
      const auto field = random_field(sys.size(), rng);
      const Outcome expected = reference_step(sys, field, c.allow_crossing);
      const StepStats stats = motion.step(sys, field);
      moved += stats.moved;
      rejected += stats.rejected;
      blocked += stats.blocked;
      for (Site s = 0; s < sys.size(); ++s) same = same && sys.occupant(s) == expected.occupant[s];
      for (Bead b = 0; b < sys.size(); ++b) same = same && sys.unwrapped(b) == expected.unwrapped[b];
      if (!same) break;
    }
    CHECK(same);
    CHECK(sys.validate().has_value());
  }
  CHECK(moved > 0 && rejected > 0 && blocked > 0);  // every outcome was exercised
}

// ------------------------------------------------------------------------------------------------
void test_uniform_field() {
  // Every bead attempting the same step closes every lattice row into a loop: everything moves.
  DllSystem sys = make_system(InitialConfiguration::Aligned, 3, 3, 30, 1);
  Motion motion(sys);
  for (std::uint8_t k = 0; k < 12; ++k) {
    std::vector<Vec3i> before(sys.size());
    for (Bead b = 0; b < sys.size(); ++b) before[b] = sys.unwrapped(b);
    const StepStats stats = motion.step(sys, std::vector<std::uint8_t>(sys.size(), k));
    CHECK(stats.moved == sys.size() && stats.rejected == 0);
    bool shifted = true;
    for (Bead b = 0; b < sys.size(); ++b) shifted = shifted && sys.unwrapped(b) == before[b] + FccGeometry::kNeighbors[k];
    CHECK(shifted);
  }
  CHECK(sys.validate().has_value());
}

void test_invariants() {
  poly::Rng rng(77);
  for (const auto kind : {InitialConfiguration::RandomWalk, InitialConfiguration::Aligned}) {
    // dense solution of tetramers / solution of octamers (a melt of straight chains hardly moves)
    DllSystem sys = kind == InitialConfiguration::Aligned ? make_system(kind, 4, 4, 48, 9) : make_system(kind, 4, 8, 12, 9);
    Motion motion(sys);
    std::size_t moved = 0;
    bool valid = true;
    for (int step = 1; step <= 4000; ++step) {
      moved += motion.step(sys, rng).moved;
      if (step % 200 == 0) valid = valid && sys.validate().has_value();
    }
    CHECK(valid);
    CHECK(moved > 0);
  }
}

// ------------------------------------------------------------------------------------------------
void test_initial_configurations() {
  CHECK(max_aligned_chains(4, 8) == 32 * 1);
  CHECK(max_aligned_chains(4, 4) == 32 * 2);
  CHECK(max_aligned_chains(5, 3) == 50 * 3);
  const DllSystem melt = make_system(InitialConfiguration::Aligned, 4, 4, 64, 1);
  CHECK(melt.chain_beads() == melt.size());  // a fully filled lattice
  CHECK(melt.validate().has_value());
  const DllSystem solution = make_system(InitialConfiguration::RandomWalk, 6, 12, 40, 2);
  CHECK(solution.validate().has_value());
  CHECK(solution.chain_beads() == 480 && solution.size() == 864);
  CHECK(solution.type(479) == BeadType::Chain && solution.type(480) == BeadType::Solvent);
  CHECK(solution.bonds(0).size() == 1 && solution.bonds(1).size() == 2 && solution.bonds(480).empty());
  bool threw = false;
  try {
    (void)make_system(InitialConfiguration::Aligned, 2, 4, 9, 1);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);
}

void test_config() {
  Config cfg;
  CHECK(cfg.set("initial", " Aligned").has_value() && cfg.initial == InitialConfiguration::Aligned);
  CHECK(!cfg.set("initial", "lattice").has_value());
  CHECK(!cfg.set("chainlength", "4").has_value());
  CHECK(cfg.set("reaction_partner", "All_Neighbors").has_value() && cfg.reaction_partner == ReactionPartner::AllNeighbors);
  CHECK(!cfg.set("reaction_partner", "nearest").has_value());
  cfg.apply_overrides({"cells=6", "n_units=300", "--chain_length=12", "n_chains=40", "run.seed=3", "dump_every=10",
                       "mechanism=FRP", "k_init=0.01", "radical_fraction=0.02", "extra_snapshots=0.255"});
  CHECK(cfg.mechanism == "frp" && cfg.n_units == 300);
  cfg.validate();
  const auto dir = std::filesystem::temp_directory_path() / "dll_test_config";
  std::filesystem::create_directories(dir);
  {
    std::ofstream out(dir / "c.ini");
    out << cfg.to_ini();
  }
  CHECK(Config::from_file(dir / "c.ini").to_ini() == cfg.to_ini());
  std::filesystem::remove_all(dir);
  bool threw = false;
  try {
    Config bad;
    bad.cells = 2;  // 32 sites for the default number of units
    bad.validate();
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);
}


// ---- phase 2: units, reactions, output ----------------------------------------------------------

using poly::chem::CompetitionRule;
using poly::chem::RateConstants;
using poly::chem::ReactionProgress;
using poly::chem::Species;

constexpr UnitState kMonomer = UnitState::monovinyl_monomer();
constexpr UnitState kCrosslinker = UnitState::divinyl_monomer();
constexpr UnitState kRadical = UnitState::primary_radical();
constexpr UnitState kDormant = UnitState::dormant_initiator();

// Units on a ring of six sites along (1, 1, 0) in a box of 3 cells, plus one off the ring.
void test_units_and_bonds() {
  std::vector<Vec3i> points;
  for (int k = 0; k < 6; ++k) points.push_back({k, k, 0});
  points.push_back({1, 0, 1});
  const std::vector<UnitState> states(points.size(), kMonomer);
  const std::vector<Vec3i> chain{{0, 2, 0}, {1, 3, 0}, {2, 4, 0}};
  DllSystem sys(3, 3, chain, points, states);
  CHECK(sys.validate().has_value());
  CHECK(sys.n_chains() == 1 && sys.n_units() == 7 && sys.solute_beads() == 10 && sys.first_unit() == 3);
  CHECK(sys.type(2) == BeadType::Chain && sys.type(3) == BeadType::Unit && sys.type(9) == BeadType::Unit &&
        sys.type(10) == BeadType::Solvent);
  const auto u = [&](int k) { return sys.first_unit() + static_cast<Bead>(k); };

  sys.bond(u(0), u(1));  // two monomers
  sys.bond(u(2), u(3));
  sys.bond(u(1), u(2));  // two molecules
  sys.bond(u(4), u(5));
  sys.bond(u(3), u(4));
  sys.bond(u(5), u(0));  // closes the ring
  sys.bond(u(0), u(6));  // a lone unit joins
  const ConnectivityTally& t = sys.connectivity();
  CHECK(t.monomer_monomer == 3 && t.polymer_polymer == 2 && t.intramolecular == 1 && t.polymer_monomer == 1);
  CHECK(sys.molecule(u(6)) == sys.molecule(u(3)) && sys.molecule(0) == sys.molecule(2) && sys.molecule(0) != sys.molecule(u(0)));
  CHECK(sys.bonded(u(0), u(5)) && !sys.bonded(u(0), u(3)) && sys.bonds(u(0)).size() == 3);
  CHECK(sys.structural_counts()[0] == 7);
  CHECK(sys.validate().has_value());

  const auto throws = [](auto&& make) {
    try {
      make();
    } catch (const std::exception&) {
      return true;
    }
    return false;
  };
  CHECK(throws([&] { sys.bond(u(0), u(1)); }));  // already bonded
  CHECK(throws([&] { sys.bond(0, u(2)); }));     // chain beads do not react
  CHECK(throws([&] { DllSystem(3, 3, chain, std::vector<Vec3i>{{0, 2, 0}}, std::vector<UnitState>{kMonomer}); }));
  const std::vector<BondPair> far{{0, 1}};
  CHECK(throws([&] { DllSystem(3, 0, {}, std::vector<Vec3i>{{0, 0, 0}, {2, 2, 0}}, std::vector<UnitState>{kMonomer, kMonomer}, far); }));
  const std::vector<BondPair> near{{0, 1}};
  const DllSystem pair(3, 0, {}, std::vector<Vec3i>{{0, 0, 0}, {1, 1, 0}}, std::vector<UnitState>{kRadical, kMonomer}, near);
  CHECK(pair.validate().has_value() && pair.bonded(0, 1) && pair.connectivity().monomer_monomer == 0);
}

// Fraction of trials in which one reaction stage makes at least one reaction, for units in the
// given states at the given points of a box of 3 cells.
double stage_probability(std::string_view mechanism, RateConstants k, ReactionPartner partner, const std::vector<Vec3i>& points,
                         const std::vector<UnitState>& states, int trials, std::uint64_t seed) {
  const auto scheme = poly::chem::make_scheme(mechanism, k, CompetitionRule::Partition);
  poly::Rng rng(seed);
  int hits = 0;
  for (int t = 0; t < trials; ++t) {
    DllSystem sys(3, 0, {}, points, states);
    Reactor reactor(sys, *scheme, partner);
    ReactionProgress progress;
    progress.reset(scheme->category_labels().size());
    hits += reactor.step(sys, rng, progress) > 0 ? 1 : 0;
  }
  return static_cast<double>(hits) / trials;
}

void test_reaction_probabilities() {
  // With k the rate constant and 1/12 the chance to sample a given neighbor:
  //   random_neighbor: radical + monomer react if the radical samples the monomer or the monomer
  //   samples the radical, p = k/12 + (1 - k/12) k/12; two radicals likewise, 1 - (1 - k/12)^2.
  //   all_neighbors: every neighboring pair is tried once, p = k.
  const double k = 0.6, q = k / 12;
  const int trials = 40000;
  const std::vector<Vec3i> pair{{0, 0, 0}, {1, 1, 0}};
  const auto near = [&](double p, double expected) {
    return std::fabs(p - expected) < 5 * std::sqrt(expected * (1 - expected) / trials);
  };
  const RateConstants kp{0, k, 0, 0}, kt{0, 0, 0, k}, ki{0.3, 0, 0, 0};
  const std::vector<UnitState> rm{kRadical, kMonomer}, rr{kRadical, kRadical}, lone{kDormant};
  CHECK(near(stage_probability("atrp", kp, ReactionPartner::RandomNeighbor, pair, rm, trials, 1), q + (1 - q) * q));
  CHECK(near(stage_probability("atrp", kp, ReactionPartner::AllNeighbors, pair, rm, trials, 2), k));
  CHECK(near(stage_probability("atrp", kp, ReactionPartner::RandomNeighbor, pair, {kMonomer, kRadical}, trials, 3), q + (1 - q) * q));
  CHECK(near(stage_probability("frp", kt, ReactionPartner::RandomNeighbor, pair, rr, trials, 4), 1 - (1 - q) * (1 - q)));
  CHECK(near(stage_probability("frp", kt, ReactionPartner::AllNeighbors, pair, rr, trials, 5), k));
  for (const auto partner : {ReactionPartner::RandomNeighbor, ReactionPartner::AllNeighbors}) {
    CHECK(near(stage_probability("frp", ki, partner, {{0, 0, 0}}, lone, trials, 6), 0.3));
  }
  CHECK(stage_probability("atrp", kp, ReactionPartner::AllNeighbors, pair, {kMonomer, kCrosslinker}, 100, 7) == 0.0);

  // A radical passed on during a stage acts only in the next stage.
  const auto atrp = poly::chem::make_scheme("atrp", RateConstants{0, 1, 1, 0}, CompetitionRule::Partition);
  DllSystem chain(3, 0, {}, std::vector<Vec3i>{{0, 0, 0}, {1, 1, 0}, {2, 2, 0}}, std::vector<UnitState>{kRadical, kMonomer, kMonomer});
  Reactor reactor(chain, *atrp, ReactionPartner::AllNeighbors);
  ReactionProgress progress;
  progress.reset(atrp->category_labels().size());
  progress.count_unit(kRadical);
  progress.count_unit(kMonomer);
  progress.count_unit(kMonomer);
  poly::Rng rng(8);
  CHECK(reactor.step(chain, rng, progress) == 1 && chain.bonded(0, 1) && !chain.bonded(1, 2));
  CHECK(chain.state(1).has_radical() && !chain.state(0).has_radical());
  CHECK(reactor.step(chain, rng, progress) == 1 && chain.bonded(1, 2) && chain.state(2).has_radical());
  CHECK(reactor.step(chain, rng, progress) == 0 && progress.propagations == 2 && progress.conversion() == 1.0);
  CHECK(progress.category_counts[1] == 1 && progress.category_counts[2] == 1);  // primary + monomer, chain + monomer
}

// Bookkeeping that must hold at every step of a reacting system that also moves.
bool consistent(const DllSystem& sys, const ReactionProgress& p, bool atrp, std::size_t initiators) {
  if (!sys.validate()) return false;
  std::uint64_t radicals = 0, consumed = 0, bonds = 0;
  for (const UnitState s : sys.states()) {
    radicals += s.radicals;
    consumed += static_cast<std::uint64_t>(poly::chem::max_vinyls(s.species) - s.vinyls) + (s.species == Species::Initiator ? 1 : 0);
  }
  for (Bead b = sys.first_unit(); b < sys.solute_beads(); ++b) bonds += sys.bonds(b).size();
  if (radicals != p.live_radicals || consumed != p.consumed_vinyls || bonds != 2 * (p.propagations + p.terminations)) return false;
  if (atrp) return radicals == initiators;
  return radicals + 2 * p.terminations == p.initiations;
}

void test_reacting_systems() {
  for (const std::string mechanism : {"atrp", "frp"}) {
    for (const auto partner : {ReactionPartner::RandomNeighbor, ReactionPartner::AllNeighbors}) {
      const bool atrp = mechanism == "atrp";
      const auto scheme = poly::chem::make_scheme(mechanism, atrp ? RateConstants{0, 0.5, 0.25, 0} : RateConstants{0.02, 0.5, 0.25, 0.9},
                                                  CompetitionRule::Legacy);
      poly::Rng rng(atrp ? 11 : 12);
      const auto states = scheme->initial_composition(600, 0.9, 0.03);
      const auto chains = initial_chains(InitialConfiguration::RandomWalk, 6, 8, 10, rng);
      DllSystem sys(6, 8, chains, random_free_points(6, chains, states.size(), rng), states);
      ReactionProgress progress;
      progress.reset(scheme->category_labels().size());
      for (const UnitState s : sys.states()) progress.count_unit(s);
      Motion motion(sys);
      Reactor reactor(sys, *scheme, partner);
      bool ok = consistent(sys, progress, atrp, 18);
      for (int step = 0; step < 3000 && ok; ++step) {
        motion.step(sys, rng);
        reactor.step(sys, rng, progress);
        if (step % 50 == 0) ok = consistent(sys, progress, atrp, 18);
      }
      CHECK(ok && consistent(sys, progress, atrp, 18));
      CHECK(progress.conversion() > (atrp ? 0.6 : 0.2));
      if (!atrp) CHECK(progress.initiations > 0 && progress.terminations > 0);
      std::uint64_t events = 0;
      for (const auto c : progress.category_counts) events += c;
      CHECK(events == progress.propagations + progress.terminations);
    }
  }
}

// Living polymerization of a monovinyl monomer: every chain grows from one initiator and stays linear.
void test_living_chains() {
  const auto scheme = poly::chem::make_scheme("atrp", RateConstants{0, 0.5, 0, 0}, CompetitionRule::Legacy);
  poly::Rng rng(21);
  const auto states = scheme->initial_composition(864, 1.0, 0.05);  // bulk, 43 initiators
  DllSystem sys(6, 0, {}, random_free_points(6, {}, states.size(), rng), states);
  ReactionProgress progress;
  progress.reset(scheme->category_labels().size());
  for (const UnitState s : sys.states()) progress.count_unit(s);
  Motion motion(sys);
  Reactor reactor(sys, *scheme, ReactionPartner::RandomNeighbor);
  for (int step = 0; step < 20000 && progress.conversion() < 0.95; ++step) {
    motion.step(sys, rng);
    reactor.step(sys, rng, progress);
  }
  CHECK(progress.conversion() >= 0.95);
  std::map<Bead, int> initiators, size;
  bool linear = true;
  for (Bead b = 0; b < sys.solute_beads(); ++b) {
    linear = linear && sys.bonds(b).size() <= 2;
    ++size[sys.molecule(b)];
    if (sys.state(b).species == Species::Initiator) ++initiators[sys.molecule(b)];
  }
  CHECK(linear);
  bool one_each = true;
  for (const auto& [root, n] : size) one_each = one_each && (n == 1 || initiators[root] == 1);
  CHECK(one_each && initiators.size() == 43);
  CHECK(sys.connectivity().polymer_polymer == 0 && sys.connectivity().intramolecular == 0);
}

void test_output_and_restart() {
  const auto scheme = poly::chem::make_scheme("atrp", RateConstants{0, 0.5, 0.5, 0}, CompetitionRule::Legacy);
  poly::Rng rng(31);
  const auto states = scheme->initial_composition(400, 0.8, 0.05);
  const auto chains = initial_chains(InitialConfiguration::RandomWalk, 5, 6, 5, rng);
  DllSystem sys(5, 6, chains, random_free_points(5, chains, states.size(), rng), states);
  ReactionProgress progress;
  progress.reset(scheme->category_labels().size());
  Motion motion(sys);
  Reactor reactor(sys, *scheme, ReactionPartner::RandomNeighbor);
  for (int step = 0; step < 500; ++step) {
    motion.step(sys, rng);
    reactor.step(sys, rng, progress);
  }
  CHECK(progress.propagations > 100);

  const auto dir = std::filesystem::temp_directory_path() / "dll_test_output";
  std::filesystem::create_directories(dir);
  write_configuration(dir / "c.dat", sys);
  const DllSystem back = read_configuration(dir / "c.dat");
  bool same = back.size() == sys.size() && back.solute_beads() == sys.solute_beads() && back.n_chains() == sys.n_chains();
  for (Bead b = 0; same && b < sys.solute_beads(); ++b) {
    auto mine = std::vector<Bead>(sys.bonds(b).begin(), sys.bonds(b).end());
    auto theirs = std::vector<Bead>(back.bonds(b).begin(), back.bonds(b).end());
    std::ranges::sort(mine);
    std::ranges::sort(theirs);
    same = sys.unwrapped(b) == back.unwrapped(b) && mine == theirs &&
           (sys.type(b) != BeadType::Unit || sys.state(b) == back.state(b));
  }
  CHECK(same && back.validate().has_value());

  write_psf(dir / "c.psf", sys);
  write_pdb(dir / "c.pdb", sys);
  std::ifstream psf(dir / "c.psf");
  std::size_t atoms = 0, bonds = 0;
  std::string tag;
  psf >> atoms >> tag;
  for (std::string line; std::getline(psf, line);) {
    if (line.contains("!NBOND")) bonds = std::stoul(line);
  }
  std::size_t expected_bonds = 0;
  for (Bead b = 0; b < sys.solute_beads(); ++b) expected_bonds += sys.bonds(b).size();
  CHECK(atoms == sys.solute_beads() && bonds == expected_bonds);
  std::ifstream pdb(dir / "c.pdb");
  std::size_t atom_lines = 0;
  for (std::string line; std::getline(pdb, line);) atom_lines += line.starts_with("ATOM") ? 1 : 0;
  CHECK(atom_lines == sys.solute_beads());
  std::filesystem::remove_all(dir);
}

void test_driver() {
  const auto dir = std::filesystem::temp_directory_path() / "dll_test_driver";
  std::filesystem::remove_all(dir);
  Config cfg;
  cfg.apply_overrides({"cells=5", "n_units=450", "monovinyl_fraction=0.9", "radical_fraction=0.04", "n_chains=4", "chain_length=6",
                       "k_p1=0.5", "k_p2=0.25", "seed=5", "equilibration_steps=20", "stop_conversion=0.5", "snapshot_start=0.2",
                       "snapshot_step=0.1", "dump_every=100", "samples=2", "directory=" + dir.string()});
  std::ostringstream log;
  const DllSystem last = run(cfg, log);
  CHECK(last.validate().has_value());
  for (const char* name : {"run_config.ini", "equilibrated_2.dat", "final_2.dat", "structural_1.dat", "propagation_1.dat",
                           "reactype_2.dat", "config_1_0.200000.psf", "config_2_0.500000.pdb", "trajectory_1.lammpstrj"}) {
    CHECK(std::filesystem::exists(dir / name));
  }
  CHECK(log.str().contains("stop_conversion reached"));
  // continue from the final configuration of sample 1
  Config again = cfg;
  again.apply_overrides({"initial_configuration=" + (dir / "final_1.dat").string(), "stop_conversion=0.6", "samples=1",
                         "directory=" + (dir / "again").string()});
  std::ostringstream log2;
  CHECK(run(again, log2).validate().has_value() && log2.str().contains("loading configuration"));
  std::filesystem::remove_all(dir);
}

}  // namespace

int main() {
  test_crossing_geometry();
  test_initial_configurations();
  test_uniform_field();
  test_against_reference();
  test_invariants();
  test_config();
  test_units_and_bonds();
  test_reaction_probabilities();
  test_reacting_systems();
  test_living_chains();
  test_output_and_restart();
  test_driver();
  std::cout << g_checks - g_failures << " / " << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
