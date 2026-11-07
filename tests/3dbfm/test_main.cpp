// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Self-contained unit tests (no framework dependency). The most important one checks the
// generic RadicalChainScheme rules against the complete type-pair transition table of the
// original implementation.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "bfm/core/index_set.hpp"
#include "bfm/core/random.hpp"
#include "bfm/io/config.hpp"
#include "bfm/io/writers.hpp"
#include "bfm/lattice/bond_set.hpp"
#include "bfm/lattice/lattice.hpp"
#include "bfm/model/system.hpp"
#include "bfm/reaction/scheme.hpp"
#include "bfm/sim/simulation.hpp"

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

#define CHECK_NEAR(a, b, tol)                                                                        \
  do {                                                                                               \
    ++g_checks;                                                                                      \
    const double _a = (a), _b = (b);                                                                 \
    if (std::fabs(_a - _b) > (tol)) {                                                                \
      ++g_failures;                                                                                  \
      std::fprintf(stderr, "CHECK_NEAR failed: %s = %g vs %s = %g  (%s:%d)\n", #a, _a, #b, _b, __FILE__, __LINE__); \
    }                                                                                                \
  } while (0)

using namespace bfm;

// ------------------------------------------------------------------------------------------------
void test_bond_set() {
  const BondVectorSet set;
  std::set<std::tuple<int, int, int>> distinct;
  std::map<int, int> by_norm2;
  for (int i = 0; i < kNumBondVectors; ++i) {
    const Vec3i v = set.vector(static_cast<BondIndex>(i));
    distinct.insert({v.x, v.y, v.z});
    ++by_norm2[v.norm2()];
    CHECK(set.index_of(v).has_value() && *set.index_of(v) == i);
    const Vec3i r = set.vector(BondVectorSet::reverse(static_cast<BondIndex>(i)));
    CHECK(r == -v);
    for (int d = 0; d < kNumDirections; ++d) {
      const auto s = set.shifted(static_cast<BondIndex>(i), static_cast<Direction>(d));
      const auto expected = set.index_of(v + unit_step(static_cast<Direction>(d)));
      CHECK(s == expected);
    }
  }
  CHECK(distinct.size() == kNumBondVectors);
  CHECK(by_norm2[4] == 6);    // (2,0,0)
  CHECK(by_norm2[5] == 24);   // (2,1,0)
  CHECK(by_norm2[6] == 24);   // (2,1,1)
  CHECK(by_norm2[9] == 30);   // (2,2,1) + (3,0,0)
  CHECK(by_norm2[10] == 24);  // (3,1,0)
  CHECK((!set.index_of(Vec3i{1, 0, 0}).has_value()));
  CHECK((!set.index_of(Vec3i{2, 2, 2}).has_value()));
  CHECK((!set.index_of(Vec3i{3, 3, 3}).has_value()));

  // legacy ordering: the first vector produced by bdibfl() is (2,0,0)
  CHECK((set.vector(0) == Vec3i{2, 0, 0}));

  const auto offsets = set.reaction_offsets();
  CHECK(offsets.size() == kNumReactionOffsets);
  std::set<std::tuple<int, int, int>> distinct_offsets;
  for (const auto& o : offsets) {
    distinct_offsets.insert({o.offset.x, o.offset.y, o.offset.z});
    const int l2 = o.offset.norm2();
    CHECK(l2 >= 4 && l2 <= 6);
    CHECK(set.vector(o.bond_from_partner_to_center) == -o.offset);
  }
  CHECK(distinct_offsets.size() == kNumReactionOffsets);
  // legacy enumeration order: lexicographic in (dx,dy,dz); the first is (-2,-1,-1)
  CHECK((offsets[0].offset == Vec3i{-2, -1, -1}));
}

// ------------------------------------------------------------------------------------------------
void test_random() {
  // Same seed gives the same stream; a different seed gives a different one.
  Xoshiro256pp a(12345), b(12345), c(12346);
  CHECK(a() == b());
  CHECK(a() != c());

  // below2: both halves uniform and independent (chi-square style tolerance on counts)
  Rng rng(42);
  constexpr int n1 = 7, n2 = 54, draws = 378000;
  std::vector<int> h1(n1, 0), h2(n2, 0), joint(n1 * n2, 0);
  for (int i = 0; i < draws; ++i) {
    const auto [x, y] = rng.below2(n1, n2);
    CHECK(x < static_cast<std::uint32_t>(n1));
    CHECK(y < static_cast<std::uint32_t>(n2));
    ++h1[x];
    ++h2[y];
    ++joint[x * n2 + y];
  }
  const auto max_rel_dev = [](const std::vector<int>& h, double expected) {
    double worst = 0;
    for (int v : h) worst = std::max(worst, std::fabs(v - expected) / expected);
    return worst;
  };
  CHECK(max_rel_dev(h1, draws / double(n1)) < 0.02);
  CHECK(max_rel_dev(h2, draws / double(n2)) < 0.04);
  CHECK(max_rel_dev(joint, draws / double(n1 * n2)) < 0.15);  // 1000 expected per cell

  // uniform01 mean and below() bounds
  double sum = 0;
  for (int i = 0; i < 100000; ++i) sum += rng.uniform01();
  CHECK_NEAR(sum / 100000, 0.5, 0.005);
  for (int i = 0; i < 1000; ++i) CHECK(rng.below(1) == 0);
}

// ------------------------------------------------------------------------------------------------
void test_lattice_and_moves() {
  const BondVectorSet set;
  Rng rng(12345);
  BfmSystem sys(PeriodicLattice(16, 16, 16), set);
  std::vector<UnitState> comp(300, UnitState::monovinyl_monomer());
  sys.place_randomly(comp, rng);
  sys.validate();

  // bond a few neighboring pairs via the reaction machinery to exercise bonded moves
  int bonded = 0;
  for (Index id = 0; id < sys.size() && bonded < 60; ++id) {
    for (int slot = 0; slot < kNumReactionOffsets; ++slot) {
      const Index p = sys.partner_at(id, slot);
      if (p == kNone || sys.monomer(id).nbonds >= 2 || sys.monomer(p).nbonds >= 2 || sys.monomer(id).bonded_to(p)) continue;
      sys.bond_through_offset(id, p, slot);
      ++bonded;
      break;
    }
  }
  CHECK(bonded > 10);
  sys.validate();

  std::uint64_t accepted = 0, attempts = 0;
  for (int mcs = 0; mcs < 300; ++mcs) {
    for (Index i = 0; i < sys.size(); ++i) {
      const auto id = static_cast<Index>(rng.below(static_cast<std::uint32_t>(sys.size())));
      accepted += sys.attempt_move(id, rng.random_direction()) ? 1 : 0;
      ++attempts;
    }
  }
  CHECK(accepted > 0 && accepted < attempts);
  sys.validate();  // occupancy, reciprocity and bond geometry all intact after 90k attempts

  // excluded volume: a cube cannot move into an occupied cube
  BfmSystem two(PeriodicLattice(8, 8, 8), set);
  Monomer a, b;
  a.set_anchor({0, 0, 0});
  b.set_anchor({2, 0, 0});
  a.state = b.state = UnitState::monovinyl_monomer();
  two.load({a, b});
  CHECK(two.lattice().occupied_count() == 16);
  CHECK((two.lattice()[1, 1, 1] && !two.lattice()[4, 0, 0]));
  CHECK(!two.attempt_move(0, Direction::PlusX));
  CHECK(!two.attempt_move(1, Direction::MinusX));
  CHECK(two.attempt_move(0, Direction::PlusY));
  CHECK((two.monomer(0).anchor() == Vec3i{0, 1, 0}));
  two.validate();

  // crossing the periodic boundary wraps the anchor and counts the image
  CHECK(two.attempt_move(0, Direction::MinusX));
  CHECK((two.monomer(0).anchor() == Vec3i{7, 1, 0}));
  CHECK((two.unwrapped(0) == Vec3i{-1, 1, 0}));
  CHECK(two.attempt_move(0, Direction::PlusX));
  CHECK((two.unwrapped(0) == Vec3i{0, 1, 0}));
  CHECK(two.attempt_move(1, Direction::PlusZ) && two.attempt_move(1, Direction::PlusZ));
  for (int i = 0; i < 7; ++i) CHECK(two.attempt_move(1, Direction::PlusZ));
  CHECK((two.monomer(1).anchor() == Vec3i{2, 0, 1}));
  CHECK((two.unwrapped(1) == Vec3i{2, 0, 9}));
  two.validate();

  // periodic wrap of unit_anchored_at: unit at x=6 on an 8-lattice is at offset (-2,0,0) from x=0
  BfmSystem wrap(PeriodicLattice(8, 8, 8), set);
  Monomer c, d;
  c.set_anchor({0, 0, 0});
  d.set_anchor({6, 0, 0});
  c.state = d.state = UnitState::monovinyl_monomer();
  wrap.load({c, d});
  bool found = false;
  const auto offsets = set.reaction_offsets();
  const Vec3i minus_two_x{-2, 0, 0};
  for (int slot = 0; slot < kNumReactionOffsets; ++slot) {
    if (offsets[static_cast<std::size_t>(slot)].offset == minus_two_x) {
      found = (wrap.partner_at(0, slot) == 1);
      wrap.bond_through_offset(0, 1, slot);
    }
  }
  CHECK(found);
  wrap.validate();
}

// ------------------------------------------------------------------------------------------------
void test_unit_codes() {
  CHECK(legacy_type_code(UnitState::monovinyl_monomer(), 0) == 0);
  CHECK((legacy_type_code(UnitState{Species::Monovinyl, 1, 0}, 1) == 1));
  CHECK((legacy_type_code(UnitState{Species::Monovinyl, 0, 0}, 2) == 2));
  CHECK(legacy_type_code(UnitState::divinyl_monomer(), 0) == 3);
  CHECK((legacy_type_code(UnitState{Species::Divinyl, 1, 1}, 1) == 4));
  CHECK((legacy_type_code(UnitState{Species::Divinyl, 0, 1}, 1) == 5));
  CHECK((legacy_type_code(UnitState{Species::Divinyl, 1, 0}, 2) == 6));
  CHECK((legacy_type_code(UnitState{Species::Divinyl, 2, 0}, 2) == 7));
  CHECK((legacy_type_code(UnitState{Species::Divinyl, 0, 0}, 4) == 8));
  CHECK(legacy_type_code(UnitState::dormant_initiator(), 0) == 9);
  CHECK(legacy_type_code(UnitState::primary_radical(), 0) == 1);
  CHECK((legacy_type_code(UnitState{Species::Initiator, 0, 0}, 1) == 2));  // consumed initiator
  CHECK(UnitState::divinyl_monomer().is_unreacted_monomer());
  CHECK((!UnitState{Species::Divinyl, 0, 1}.is_unreacted_monomer()));
  CHECK((UnitState{Species::Divinyl, 0, 1}.has_pendant_vinyl()));
  CHECK(legacy_residue_name(9) == "INI");
}

// ------------------------------------------------------------------------------------------------
// Representative unit for each legacy code (nbonds chosen as they occur in a real run).
UnitView view_for(int code, bool fresh_primary = false) {
  switch (code) {
    case 0: return {UnitState::monovinyl_monomer(), 0};
    case 1: return fresh_primary ? UnitView{UnitState::primary_radical(), 0} : UnitView{UnitState{Species::Monovinyl, 1, 0}, 1};
    case 2: return {UnitState{Species::Monovinyl, 0, 0}, 2};
    case 3: return {UnitState::divinyl_monomer(), 0};
    case 4: return {UnitState{Species::Divinyl, 1, 1}, 1};
    case 5: return {UnitState{Species::Divinyl, 0, 1}, 1};
    case 6: return {UnitState{Species::Divinyl, 1, 0}, 2};
    case 7: return {UnitState{Species::Divinyl, 2, 0}, 2};
    case 8: return {UnitState{Species::Divinyl, 0, 0}, 4};
    case 9: return {UnitState::dormant_initiator(), 0};
  }
  return {};
}

int code_after(UnitView v, int radicals_delta, int vinyls_delta, int bonds_delta) {
  UnitState s = v.state;
  s.radicals = static_cast<std::uint8_t>(s.radicals + radicals_delta);
  s.vinyls = static_cast<std::uint8_t>(s.vinyls + vinyls_delta);
  return legacy_type_code(s, v.nbonds + bonds_delta);
}

enum class Cls { None, First, Second, Term, Competition };

struct LegacyRow {
  int a, b;             // legacy montype of mono (center) and nu1 (partner)
  int prod_a, prod_b;   // products of the propagation (or termination when Cls::Term)
  int term_a, term_b;   // products of the termination branch when Cls::Competition
  Cls cls;
};

// Transition table of the original FRP implementation (ATRP: the same without Term/Competition rows).
const std::vector<LegacyRow> kLegacyTable = {
    {0, 1, 1, 2, 0, 0, Cls::First}, {0, 4, 1, 5, 0, 0, Cls::First}, {0, 6, 1, 8, 0, 0, Cls::First}, {0, 7, 1, 6, 0, 0, Cls::First},
    {1, 0, 2, 1, 0, 0, Cls::First}, {1, 3, 2, 4, 0, 0, Cls::First}, {1, 5, 2, 6, 0, 0, Cls::Second},
    {1, 1, 2, 2, 0, 0, Cls::Term},  {1, 6, 2, 8, 0, 0, Cls::Term},  {1, 7, 2, 6, 0, 0, Cls::Term},
    {1, 4, 2, 7, 2, 5, Cls::Competition},
    {3, 1, 4, 2, 0, 0, Cls::First}, {3, 4, 4, 5, 0, 0, Cls::First}, {3, 6, 4, 8, 0, 0, Cls::First}, {3, 7, 4, 6, 0, 0, Cls::First},
    {4, 0, 5, 1, 0, 0, Cls::First}, {4, 3, 5, 4, 0, 0, Cls::First}, {4, 5, 5, 6, 0, 0, Cls::Second},
    {4, 1, 7, 2, 5, 2, Cls::Competition}, {4, 4, 5, 7, 5, 5, Cls::Competition}, {4, 6, 7, 8, 5, 8, Cls::Competition},
    {4, 7, 7, 6, 5, 6, Cls::Competition},
    {5, 1, 6, 2, 0, 0, Cls::Second}, {5, 4, 6, 5, 0, 0, Cls::Second}, {5, 6, 6, 8, 0, 0, Cls::Second}, {5, 7, 6, 6, 0, 0, Cls::Second},
    {6, 0, 8, 1, 0, 0, Cls::First}, {6, 3, 8, 4, 0, 0, Cls::First}, {6, 5, 8, 6, 0, 0, Cls::Second},
    {6, 1, 8, 2, 0, 0, Cls::Term},  {6, 6, 8, 8, 0, 0, Cls::Term},  {6, 7, 8, 6, 0, 0, Cls::Term},
    {6, 4, 8, 7, 8, 5, Cls::Competition},
    {7, 0, 6, 1, 0, 0, Cls::First}, {7, 3, 6, 4, 0, 0, Cls::First}, {7, 5, 6, 6, 0, 0, Cls::Second},
    {7, 1, 6, 2, 0, 0, Cls::Term},  {7, 6, 6, 8, 0, 0, Cls::Term},  {7, 7, 6, 6, 0, 0, Cls::Term},
    {7, 4, 6, 7, 6, 5, Cls::Competition},
};

struct Outcome {
  int propagations{0};
  int terminations{0};
  int prod_a{-1}, prod_b{-1};  // products of the first propagation seen
  int term_a{-1}, term_b{-1};  // products of the first termination seen
};

Outcome sample_pair(const ReactionScheme& scheme, int ca, int cb, int trials, Rng& rng) {
  Outcome out;
  const UnitView va = view_for(ca), vb = view_for(cb);
  for (int t = 0; t < trials; ++t) {
    const auto ev = scheme.propose_pair(0, va, 1, vb, rng);
    if (!ev) continue;
    if (ev->kind == EventKind::Propagation) {
      ++out.propagations;
      const bool a_attacks = (ev->first == 0);
      const int pa = a_attacks ? code_after(va, -1, 0, 1) : code_after(va, +1, -1, 1);
      const int pb = a_attacks ? code_after(vb, +1, -1, 1) : code_after(vb, -1, 0, 1);
      if (out.prod_a < 0) { out.prod_a = pa; out.prod_b = pb; }
    } else if (ev->kind == EventKind::Termination) {
      ++out.terminations;
      if (out.term_a < 0) { out.term_a = code_after(va, -1, 0, 1); out.term_b = code_after(vb, -1, 0, 1); }
    }
  }
  return out;
}

void test_scheme_fidelity() {
  const RateConstants k{0.0, 0.7, 0.3, 0.5};  // k_p2 + k_t = 0.8 < 1 -> legacy p_prop = k_p2
  const auto frp = make_scheme("frp", k, CompetitionRule::Legacy);
  const auto atrp = make_scheme("atrp", k, CompetitionRule::Legacy);
  Rng rng(2024);
  const int trials = 40000;
  const double tol = 0.015;

  std::set<std::pair<int, int>> in_table;
  for (const LegacyRow& row : kLegacyTable) {
    in_table.insert({row.a, row.b});
    const Outcome f = sample_pair(*frp, row.a, row.b, trials, rng);
    const Outcome a = sample_pair(*atrp, row.a, row.b, trials, rng);
    const double pf = static_cast<double>(f.propagations) / trials;
    const double tf = static_cast<double>(f.terminations) / trials;
    const double pa = static_cast<double>(a.propagations) / trials;
    switch (row.cls) {
      case Cls::First:
        CHECK_NEAR(pf, k.k_p1, tol);
        CHECK(f.terminations == 0);
        CHECK_NEAR(pa, k.k_p1, tol);
        CHECK(f.prod_a == row.prod_a && f.prod_b == row.prod_b);
        CHECK(a.prod_a == row.prod_a && a.prod_b == row.prod_b);
        break;
      case Cls::Second:
        CHECK_NEAR(pf, k.k_p2, tol);
        CHECK(f.terminations == 0);
        CHECK_NEAR(pa, k.k_p2, tol);
        CHECK(f.prod_a == row.prod_a && f.prod_b == row.prod_b);
        break;
      case Cls::Term:
        CHECK(f.propagations == 0);
        CHECK_NEAR(tf, k.k_term, tol);
        CHECK(a.propagations == 0 && a.terminations == 0);  // ATRP: no termination
        CHECK(f.term_a == row.prod_a && f.term_b == row.prod_b);
        break;
      case Cls::Competition: {
        const double p = k.k_p2;  // sum < 1
        const double expected_term = (1 - p) * (1 - p) * std::min(1.0, k.k_term + p);
        CHECK_NEAR(pf, p, tol);
        CHECK_NEAR(tf, expected_term, tol);
        CHECK(f.prod_a == row.prod_a && f.prod_b == row.prod_b);
        CHECK(f.term_a == row.term_a && f.term_b == row.term_b);
        // ATRP: the same pair is a plain pendant propagation
        CHECK_NEAR(pa, k.k_p2, tol);
        CHECK(a.terminations == 0);
        CHECK(a.prod_a == row.prod_a && a.prod_b == row.prod_b);
        break;
      }
      case Cls::None: break;
    }
  }
  // every pair not in the table must be inert in both mechanisms
  for (int ca = 0; ca <= 9; ++ca) {
    for (int cb = 0; cb <= 9; ++cb) {
      if (in_table.count({ca, cb})) continue;
      const Outcome f = sample_pair(*frp, ca, cb, 2000, rng);
      const Outcome a = sample_pair(*atrp, ca, cb, 2000, rng);
      CHECK(f.propagations == 0 && f.terminations == 0);
      CHECK(a.propagations == 0 && a.terminations == 0);
    }
  }

  // Paper constants: k_p2 + k_t = 1.2 >= 1 -> p_prop = 0.25/1.2, P(term) = (1-p)^2 * 1
  const RateConstants paper{1e-5, 0.5, 0.25, 0.95};
  const auto frp_paper = make_scheme("frp", paper, CompetitionRule::Legacy);
  const Outcome c = sample_pair(*frp_paper, 1, 4, trials, rng);
  const double p = 0.25 / 1.2;
  CHECK_NEAR(static_cast<double>(c.propagations) / trials, p, tol);
  CHECK_NEAR(static_cast<double>(c.terminations) / trials, (1 - p) * (1 - p), tol);

  // Partition rule: same pair, probabilities renormalized to k_p2/(sum), k_t/(sum)
  const auto frp_part = make_scheme("frp", paper, CompetitionRule::Partition);
  const Outcome d = sample_pair(*frp_part, 1, 4, trials, rng);
  CHECK_NEAR(static_cast<double>(d.propagations) / trials, 0.25 / 1.2, tol);
  CHECK_NEAR(static_cast<double>(d.terminations) / trials, 0.95 / 1.2, tol);

  // Initiation: dormant initiator only, with k_init
  const RateConstants ki{0.3, 0.5, 0.25, 0.95};
  const auto frp_init = make_scheme("frp", ki, CompetitionRule::Legacy);
  int inits = 0;
  for (int t = 0; t < trials; ++t) inits += frp_init->propose_unimolecular(0, view_for(9), rng).has_value() ? 1 : 0;
  CHECK_NEAR(static_cast<double>(inits) / trials, 0.3, tol);
  CHECK(!frp_init->propose_unimolecular(0, view_for(0), rng).has_value());
  CHECK(!atrp->has_unimolecular_step());
  CHECK(frp->has_unimolecular_step());

  // termination adds 15 categories to the 20 propagation categories
  CHECK(frp->category_labels().size() == 35);
  CHECK(atrp->category_labels().size() == 20);
  bool unknown = false;
  try {
    (void)make_scheme("unknown", k, CompetitionRule::Legacy);
  } catch (const std::invalid_argument&) {
    unknown = true;
  }
  CHECK(unknown);
}

// ------------------------------------------------------------------------------------------------
void test_categories() {
  using S = RadicalChainScheme;
  CHECK(S::propagation_category(view_for(1, true), view_for(0)) == 1);
  CHECK(S::propagation_category(view_for(1), view_for(0)) == 2);
  CHECK(S::propagation_category(view_for(4), view_for(0)) == 3);
  CHECK(S::propagation_category(view_for(6), view_for(0)) == 4);
  CHECK(S::propagation_category(view_for(7), view_for(0)) == 5);
  CHECK(S::propagation_category(view_for(1, true), view_for(3)) == 6);
  CHECK(S::propagation_category(view_for(1), view_for(3)) == 7);
  CHECK(S::propagation_category(view_for(1), view_for(4)) == 9);
  CHECK(S::propagation_category(view_for(1, true), view_for(5)) == 10);
  CHECK(S::propagation_category(view_for(4), view_for(3)) == 12);
  CHECK(S::propagation_category(view_for(4), view_for(4)) == 15);
  CHECK(S::propagation_category(view_for(6), view_for(4)) == 17);
  CHECK(S::propagation_category(view_for(7), view_for(5)) == 20);
  CHECK(S::termination_category(view_for(1, true), view_for(1, true)) == 21);
  CHECK(S::termination_category(view_for(1, true), view_for(1)) == 22);
  CHECK(S::termination_category(view_for(1), view_for(1)) == 23);
  CHECK(S::termination_category(view_for(4), view_for(1, true)) == 24);
  CHECK(S::termination_category(view_for(1), view_for(4)) == 25);
  CHECK(S::termination_category(view_for(6), view_for(1)) == 27);
  CHECK(S::termination_category(view_for(7), view_for(1, true)) == 28);
  CHECK(S::termination_category(view_for(4), view_for(4)) == 30);
  CHECK(S::termination_category(view_for(6), view_for(4)) == 31);
  CHECK(S::termination_category(view_for(7), view_for(4)) == 32);
  CHECK(S::termination_category(view_for(6), view_for(6)) == 33);
  CHECK(S::termination_category(view_for(7), view_for(6)) == 34);
  CHECK(S::termination_category(view_for(7), view_for(7)) == 35);
}

// ------------------------------------------------------------------------------------------------
void test_composition() {
  const auto atrp = make_scheme("atrp", RateConstants{0, 0.005, 0.0025, 0}, CompetitionRule::Legacy);
  const auto comp = atrp->initial_composition(106250, 0.96, 0.002);
  std::size_t n_init = 0, n_mono = 0, n_di = 0;
  for (const UnitState& s : comp) {
    if (s.species == Species::Initiator) ++n_init;
    else if (s.species == Species::Monovinyl) ++n_mono;
    else ++n_di;
  }
  CHECK(comp.size() == 106250);
  CHECK(n_init == 212);   // int(106250 * 0.002)
  CHECK(n_mono == 101788);  // 102000 - 212
  CHECK(n_di == 4250);
  CHECK(comp[0] == UnitState::primary_radical());
  const auto frp = make_scheme("frp", RateConstants{1e-5, 0.5, 0.25, 0.95}, CompetitionRule::Legacy);
  CHECK(frp->initial_composition(1000, 0.96, 0.002)[0] == UnitState::dormant_initiator());

  // structural-unit counts of a freshly placed system: legacy codes 1 (live chains), 0, 3
  const BondVectorSet set;
  Rng rng(5);
  BfmSystem sys(PeriodicLattice(32, 32, 32), set);
  sys.place_randomly(atrp->initial_composition(3481, 0.96, 0.002), rng);
  const auto counts = sys.structural_counts();
  CHECK(counts[1] == 6);              // int(3481 * 0.002) primary radicals
  CHECK(counts[0] == 3341 - 6);       // int(3481 * 0.96) monovinyl units minus the initiators
  CHECK(counts[3] == 3481 - 3341);    // cross-linkers
  CHECK(counts[2] + counts[4] + counts[5] + counts[6] + counts[7] + counts[8] + counts[9] == 0);
}

// ------------------------------------------------------------------------------------------------
void test_config_and_restart() {
  namespace fs = std::filesystem;
  const fs::path dir = fs::temp_directory_path() / "bfm_test_cfg";
  fs::create_directories(dir);
  {
    std::ofstream f(dir / "a.ini");
    f << "[system]\nlattice = 32 32 32\nn_units = 3481   # comment\n[scheme]\nmechanism = FRP\nk_term=0.9\n"
         "[run]\nsamples = 2\nextra_snapshots = 0.969, 0.5\n";
  }
  Config cfg = Config::from_file(dir / "a.ini");
  CHECK(cfg.lattice_x == 32 && cfg.n_units == 3481 && cfg.mechanism == "frp");
  CHECK_NEAR(cfg.rates.k_term, 0.9, 1e-12);
  cfg.apply_overrides({"--seed=42", "run.samples=3", "stop_conversion=0.96"});
  CHECK(cfg.seed == 42 && cfg.samples == 3);
  cfg.validate();
  const auto levels = cfg.snapshot_levels();
  CHECK(levels.size() == 94);  // 0.03..0.96 step 0.01 (94 values); 0.5 duplicates a regular level
  CHECK_NEAR(levels.front(), 0.03, 1e-12);
  CHECK_NEAR(levels.back(), 0.96, 1e-12);
  bool threw = false;
  try {
    cfg.apply_overrides({"nonsense=1"});
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw);

  // restart round trip through the legacy format
  const BondVectorSet set;
  Rng rng(99);
  BfmSystem sys(PeriodicLattice(16, 16, 16), set);
  const auto frp = make_scheme("frp", RateConstants{1e-3, 0.5, 0.25, 0.95}, CompetitionRule::Legacy);
  sys.place_randomly(frp->initial_composition(400, 0.9, 0.01), rng);
  for (int mcs = 0; mcs < 50; ++mcs) {
    for (Index i = 0; i < sys.size(); ++i) sys.attempt_move(static_cast<Index>(rng.below(400)), rng.random_direction());
  }
  int bonded = 0;
  for (Index id = 0; id < sys.size() && bonded < 40; ++id) {
    for (int slot = 0; slot < kNumReactionOffsets; ++slot) {
      const Index p = sys.partner_at(id, slot);
      if (p == kNone || sys.monomer(id).nbonds >= 1 || sys.monomer(p).nbonds >= 1) continue;
      sys.bond_through_offset(id, p, slot);
      ++bonded;
      break;
    }
  }
  write_legacy_configuration(dir / "restart.dat", sys);
  const Vec3i box{16, 16, 16};
  const auto loaded = read_legacy_configuration(dir / "restart.dat", box, 400);
  CHECK(loaded.units.size() == 400 && loaded.unwrapped.size() == 400);
  bool same = true;
  bool any_image = false;
  for (Index i = 0; i < sys.size(); ++i) {
    const Monomer& a = sys.monomer(i);
    const Monomer& b = loaded.units[static_cast<std::size_t>(i)];
    same = same && a.anchor() == b.anchor() && sys.unwrapped(i) == loaded.unwrapped[static_cast<std::size_t>(i)] &&
           a.nbonds == b.nbonds && a.state == b.state;
    any_image = any_image || a.anchor() != sys.unwrapped(i);
    for (int k = 0; k < a.nbonds; ++k) same = same && a.bond(k).partner == b.bond(k).partner && a.bond(k).vector == b.bond(k).vector;
  }
  CHECK(same);
  CHECK(any_image);  // some units crossed the boundary in 50 MCS, so images are exercised
  BfmSystem again(PeriodicLattice(16, 16, 16), set);
  again.load(loaded.units, loaded.unwrapped);  // validates occupancy and bond geometry
  bool same_images = true;
  for (Index i = 0; i < sys.size(); ++i) same_images = same_images && again.unwrapped(i) == sys.unwrapped(i);
  CHECK(same_images);
  fs::remove_all(dir);
}


// ------------------------------------------------------------------------------------------------
// Reaction attempts restricted to reactive units must reproduce the legacy attempt distribution.
void test_reactive_sampling() {
  // IndexSet against std::set under random insertions and removals
  {
    IndexSet set(200);
    std::set<Index> ref;
    Rng rng(3);
    for (int i = 0; i < 20000; ++i) {
      const auto id = static_cast<Index>(rng.below(200));
      const bool add = rng.bernoulli(0.5);
      set.set(id, add);
      if (add) {
        ref.insert(id);
      } else {
        ref.erase(id);
      }
    }
    CHECK(set.size() == ref.size());
    std::set<Index> members;
    for (std::size_t k = 0; k < set.size(); ++k) members.insert(set[k]);
    CHECK(members == ref);
    bool consistent = true;
    for (Index id = 0; id < 200; ++id) consistent = consistent && set.contains(id) == ref.contains(id);
    CHECK(consistent);
  }

  // the offset in the opposite slot is the negated offset
  const BondVectorSet bonds;
  const auto offsets = bonds.reaction_offsets();
  bool opposite = true;
  for (int s = 0; s < kNumReactionOffsets; ++s) {
    opposite = opposite && offsets[static_cast<std::size_t>(BondVectorSet::opposite_offset(s))].offset ==
                               -offsets[static_cast<std::size_t>(s)].offset;
  }
  CHECK(opposite);

  // Scheme contract: a pair of non-reactive units never reacts, a non-reactive unit never
  // undergoes a unimolecular step, even with all rate constants at 1.
  const RateConstants ones{1.0, 1.0, 1.0, 1.0};
  Rng rng(17);
  for (const char* mech : {"frp", "atrp"}) {
    const auto scheme = make_scheme(mech, ones, CompetitionRule::Legacy);
    std::vector<UnitView> views;
    for (int code = 0; code < kNumLegacyTypeCodes; ++code) views.push_back(view_for(code));
    views.push_back(view_for(1, true));
    for (const UnitView& a : views) {
      bool uni = false;
      for (int t = 0; t < 20; ++t) uni = uni || scheme->propose_unimolecular(0, a, rng).has_value();
      CHECK(!uni || scheme->is_reactive(a));
      for (const UnitView& b : views) {
        if (scheme->is_reactive(a) || scheme->is_reactive(b)) continue;
        bool reacted = false;
        for (int t = 0; t < 20; ++t) reacted = reacted || scheme->propose_pair(0, a, 1, b, rng).has_value();
        CHECK(!reacted);
      }
    }
  }
  const auto frp = make_scheme("frp", ones, CompetitionRule::Legacy);
  const auto atrp = make_scheme("atrp", ones, CompetitionRule::Legacy);
  CHECK(frp->is_reactive(view_for(9)) && !atrp->is_reactive(view_for(9)));  // dormant initiator
  CHECK(frp->is_reactive(view_for(1)) && atrp->is_reactive(view_for(7)) && !frp->is_reactive(view_for(5)));

  // map_reaction_draw: enumerating every draw (x, slot) yields every legacy attempt that
  // involves a reactive unit exactly once, and no other attempt.
  BfmSystem sys(PeriodicLattice(16, 16, 16), bonds);
  sys.place_randomly(std::vector<UnitState>(300, UnitState::monovinyl_monomer()), rng);
  for (int mcs = 0; mcs < 20; ++mcs) {
    for (Index i = 0; i < sys.size(); ++i) sys.attempt_move(static_cast<Index>(rng.below(300)), rng.random_direction());
  }
  const auto n = static_cast<std::uint32_t>(sys.size());
  const auto enumerate = [&](const IndexSet& reactive) {
    std::vector<ReactionAttempt> got;
    for (std::uint32_t x = 0; x < n; ++x) {
      for (int s = 0; s < kNumReactionOffsets; ++s) {
        if (const auto a = map_reaction_draw(sys, reactive, x, s)) got.push_back(*a);
      }
    }
    return got;
  };
  IndexSet reactive(300);
  for (Index id = 0; id < 300; id += 7) reactive.insert(id);
  auto got = enumerate(reactive);
  std::set<std::pair<Index, int>> unique;
  for (const auto& a : got) unique.insert({a.center, a.slot});
  CHECK(unique.size() == got.size());  // no attempt produced twice
  std::set<std::pair<Index, int>> expected;
  for (Index c = 0; c < sys.size(); ++c) {
    for (int s = 0; s < kNumReactionOffsets; ++s) {
      const Index p = sys.partner_at(c, s);
      if (reactive.contains(c) || (p != kNone && reactive.contains(p))) expected.insert({c, s});
    }
  }
  CHECK(unique == expected);
  CHECK(expected.size() > 43 * 54);  // reverse attempts (partner reactive, center not) occur

  CHECK(enumerate(IndexSet(300)).empty());  // nothing reactive: nothing attempted
  IndexSet dense(300);
  for (Index id = 0; id < 160; ++id) dense.insert(id);  // 2|A| > n: plain legacy sampling
  CHECK(enumerate(dense).size() == static_cast<std::size_t>(n) * kNumReactionOffsets);
}

}  // namespace

int main() {
  test_random();
  test_bond_set();
  test_lattice_and_moves();
  test_unit_codes();
  test_scheme_fidelity();
  test_categories();
  test_composition();
  test_config_and_restart();
  test_reactive_sampling();
  std::cout << g_checks - g_failures << " / " << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
