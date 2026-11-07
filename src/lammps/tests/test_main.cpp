// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Self-contained unit tests (no framework, no LAMMPS) of the reaction files and the rule engine
// of fix poly/react. The most important one checks the radical preset against the transition
// table of the original 3DBFM code for every pair of structural units.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "polyreact/engine.hpp"
#include "polyreact/presets.hpp"
#include "polyreact/reaction_file.hpp"

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

#define CHECK_NEAR(a, b, tol) CHECK(std::fabs((a) - (b)) <= (tol))

using namespace polyreact;

struct TestRng {
  std::mt19937_64 engine;
  explicit TestRng(std::uint64_t seed) : engine(seed) {}
  double uniform() { return std::uniform_real_distribution<double>(0.0, 1.0)(engine); }
};

template <class F>
bool throws(F&& f) {
  try {
    f();
  } catch (const std::exception&) {
    return true;
  }
  return false;
}

ReactionFile parse(const std::string& text) {
  std::istringstream in(text);
  return parse_reaction_file(in, "test");
}

// ------------------------------------------------------------------------------------------------
void test_equations() {
  const Equation a = parse_equation("FA(3, 4) + MA(1,2) -> FA(1):DA(2)");
  CHECK(a.kind == Kind::Association && a.reactants.size() == 2);
  CHECK((a.reactants[0] == Reactant{"FA", 3, 4}) && (a.reactants[1] == Reactant{"MA", 1, 2}));
  CHECK((a.products[0] == Product{"FA", 1}) && (a.products[1] == Product{"DA", 2}));
  const Equation d = parse_equation("C(0, 1):E(0, 1) -> A(1) + E(-1)");
  CHECK(d.kind == Kind::Dissociation && (d.products[1] == Product{"E", -1}));
  const Equation x = parse_equation("C(0, 1):E(1, 2) + W(0, 2) -> A(1):Z(1) + D(1)");
  CHECK(x.kind == Kind::Exchange && x.reactants.size() == 3);
  // products in the order of the reactants: C -> A, E -> D, W -> Z
  CHECK(x.products[0].type == "A" && x.products[1].type == "D" && x.products[2].type == "Z");
  const Equation u = parse_equation(" INI -> PRI ");
  CHECK(u.kind == Kind::Unimolecular && u.reactants[0].min == INT_MIN && u.products[0].delta == 0);
  CHECK(u.reactants[0].accepts(-5) && u.reactants[0].accepts(7));
  CHECK(throws([] { (void)parse_equation("A + B -> C:D -> E"); }));
  CHECK(throws([] { (void)parse_equation("A(2,2) + B -> A:B"); }));
  CHECK(throws([] { (void)parse_equation("A(1) + B -> A:B"); }));
  CHECK(throws([] { (void)parse_equation("A + B -> A + B"); }));
  CHECK(throws([] { (void)parse_equation("1A + B -> A:B"); }));
  CHECK(throws([] { (void)parse_equation("A + B"); }));
}

const std::string kChemLabStyle = R"(# A file in the format of ChemLab
[general]
interval: 200
nearest=0

[ext_change_type]
ext_type=ChangeNeighboursProperty
type_transfers=MA:2->PA,ML:1->PL(state=1)

[reaction_b]
reaction: DA(3, 4) + MA(1, 2) -> RA(1):DA(2)
cutoff: 1.2
rate: 0.8
intramolecular: 1
intraresidual: 0
active: True
group: reaction_1

[group_reaction_1]
potential=Harmonic
potential_options=K=30.0,r0=0.97
extensions=change_type
bond_type = 2   # LAMMPS bond type of the new bonds

[reaction_a]
reaction = FA(3,4) + MA(1,2) -> FA(1):DA(2)
cutoff = 1.0
min_cutoff = 0.2
rate = 0.5
group = reaction_1

[reaction_c]
reaction: I(0,5):I(0,5) -> I(1) + I(1)
diss_rate: 0.000001
cutoff: 0.5
rate: 0.0
alpha: 0.00001
group: reaction_1
)";

void test_reaction_file() {
  const ReactionFile f = parse(kChemLabStyle);
  CHECK(f.general.interval == 200 && !f.general.nearest && f.general.competition == Competition::Independent);
  CHECK(f.reactions.size() == 3 && f.groups.size() == 1 && f.extensions.size() == 1);
  CHECK(f.extensions[0].type == "ChangeNeighboursProperty" && f.extensions[0].options.at("type_transfers") == "MA:2->PA,ML:1->PL(state=1)");
  const Group& g = f.groups[0];
  CHECK(g.name == "reaction_1" && g.bond_type == 2 && g.extensions == std::vector<std::string>{"change_type"});
  CHECK(g.options.at("potential_options") == "K=30.0,r0=0.97");
  const Reaction& b = f.reactions[0];
  CHECK(b.name == "b" && b.rate == 0.8 && b.cutoff == 1.2 && b.intramolecular && !b.intraresidual && b.active);
  const Reaction& a = f.reactions[1];
  CHECK(!a.intramolecular && !a.intraresidual && !a.virtual_bond);  // ChemLab defaults
  CHECK(a.min_cutoff == 0.2 && a.equation.reactants[0].type == "FA");
  CHECK(f.reactions[2].equation.kind == Kind::Dissociation && f.reactions[2].diss_rate == 1e-6);
  CHECK(f.notes.size() == 3);  // potential, potential_options, alpha
  const auto names = f.type_names();
  CHECK((names == std::vector<std::string>{"DA", "MA", "RA", "FA", "I"}));

  CHECK(throws([] { (void)parse("[general]\ninterval: 0\n"); }));
  CHECK(throws([] { (void)parse("[general]\nintervall: 10\n"); }));
  CHECK(throws([] { (void)parse("[reactions]\n"); }));
  CHECK(throws([] { (void)parse("[group_g]\n[reaction_a]\nreaction: A + B -> A:B\nrate: 0.1\ngroup: g\n"); }));  // no cutoff
  CHECK(throws([] { (void)parse("[group_g]\n[reaction_a]\nreaction: A + B -> A:B\nrate: 1.5\ncutoff: 1\ngroup: g\n"); }));
  CHECK(throws([] { (void)parse("[reaction_a]\nreaction: A + B -> A:B\nrate: 0.1\ncutoff: 1\ngroup: g\n"); }));  // no group g
  CHECK(throws([] { (void)parse("[group_g]\nextensions = missing\n"); }));
  CHECK(throws([] { (void)parse("[group_g]\n[group_g]\n"); }));
  bool located = false;
  try {
    (void)parse("[general]\ninterval: 10\n[group_g]\n[reaction_a]\nreaction: A + B -> C\n");
  } catch (const std::runtime_error& e) {
    located = std::string(e.what()).starts_with("test:5:");
  }
  CHECK(located);
}

// ------------------------------------------------------------------------------------------------
void test_engine_matching() {
  const ReactionFile f = parse(R"([general]
competition = partition
[group_g]
[reaction_bond]
reaction = A(1,2) + B(1,3) -> C(-1):D(-1)
rate = 1
cutoff = 1.5
min_cutoff = 0.5
intramolecular = no
group = g
[reaction_cut]
reaction = C(0,1):D(0,2) -> A(1) + B(1)
rate = 1
cutoff = 2
group = g
[reaction_mark]
reaction = A(1,2) + E -> A:F
virtual = true
intramolecular = yes
rate = 1
cutoff = 1
group = g
[reaction_decay]
reaction = E(0,1) -> F(3)
rate = 1
group = g
)");
  const Engine e(f);
  const int A = e.type_index("A"), B = e.type_index("B"), C = e.type_index("C"), D = e.type_index("D"), E = e.type_index("E"),
            F = e.type_index("F");
  CHECK(e.num_types() == 6 && e.max_cutoff() == 2 && e.has_unimolecular());
  CHECK(e.reactive(A) && e.reactive(E) && !e.reactive(F));
  CHECK(e.seeks(A) && e.seeks(C) && !e.seeks(B) && !e.seeks(E) && !e.seeks(F));
  CHECK(e.may_react(A, B) && e.may_react(B, A) && e.may_react(E, A) && !e.may_react(A, A) && !e.may_react(B, E) && !e.may_react(-1, A));
  TestRng rng(1);
  const PairContext near{1.0, false, false, false};
  auto r = e.pair({A, 1}, {B, 2}, near, rng);
  CHECK(r && r->creates_bond && !r->swapped && r->center == (Bead{C, 0}) && r->partner == (Bead{D, 1}) && r->reaction == 0);
  r = e.pair({B, 1}, {A, 1}, near, rng);  // the partner plays the first bead
  CHECK(r && r->swapped && r->center == (Bead{D, 0}) && r->partner == (Bead{C, 0}));
  CHECK(!e.pair({A, 2}, {B, 2}, near, rng));                                     // functionality outside [1, 2)
  CHECK(!e.pair({A, 1}, {B, 3}, near, rng));                                     // outside [1, 3)
  CHECK(!e.pair({A, 1}, {B, 1}, {1.5, false, false, false}, rng));               // at the cutoff
  CHECK(!e.pair({A, 1}, {B, 1}, {0.4, false, false, false}, rng));               // below min_cutoff
  CHECK(!e.pair({A, 1}, {B, 1}, {1.0, true, false, false}, rng));                // already bonded
  CHECK(!e.pair({A, 1}, {B, 1}, {1.0, false, true, false}, rng));                // same molecule
  CHECK(!e.pair({A, 1}, {B, 1}, {1.0, false, false, true}, rng));                // same residue
  r = e.pair({C, 0}, {D, 1}, {1.8, true, true, true}, rng);                       // dissociation needs the bond
  CHECK(r && r->breaks_bond && r->center == (Bead{A, 1}) && r->partner == (Bead{B, 2}));
  CHECK(!e.pair({C, 0}, {D, 1}, {1.8, false, false, false}, rng));
  r = e.pair({E, 0}, {A, 1}, {0.9, true, true, false}, rng);                      // virtual: no bond, bonded allowed
  CHECK(r && !r->creates_bond && r->center == (Bead{F, 0}) && r->partner == (Bead{A, 1}));
  r = e.single({E, 0}, rng);
  CHECK(r && r->kind == Kind::Unimolecular && r->center == (Bead{F, 3}) && r->reaction == 3);
  CHECK(!e.single({E, 1}, rng) && !e.single({A, 1}, rng) && !e.pair({7, 0}, {A, 1}, near, rng));

  // Types given by the caller, e.g. LAMMPS types; unknown names are an error.
  const std::map<std::string, int, std::less<>> lammps{{"A", 0}, {"B", 2}, {"C", 3}, {"D", 4}, {"E", 5}, {"F", 6}};
  const Engine mapped(f, lammps, 8);
  CHECK(mapped.num_types() == 8 && mapped.reactive(0) && !mapped.reactive(1) && !mapped.reactive(7));
  CHECK(throws([&] { (void)Engine(f, {{"A", 0}}, 2); }));
  CHECK(throws([] { (void)Engine(parse("[group_g]\n[reaction_x]\nreaction = A:B + C -> A:C + B\nrate = 0.1\ncutoff = 1\ngroup = g\n")); }));
}

/// Frequencies of the outcomes of reaction 0 and reaction 1 for one pair.
std::pair<double, double> frequencies(const Engine& e, Bead a, Bead b, int trials, TestRng& rng) {
  int n0 = 0, n1 = 0;
  for (int t = 0; t < trials; ++t) {
    const auto r = e.pair(a, b, {1.0, false, false, false}, rng);
    if (r) (r->reaction == 0 ? n0 : n1) += 1;
  }
  return {static_cast<double>(n0) / trials, static_cast<double>(n1) / trials};
}

void test_competition() {
  const int trials = 200000;
  const double tol = 0.006;
  TestRng rng(2);
  const auto file = [](std::string_view competition, double a, double b) {
    return parse(std::format("[general]\ncompetition = {}\n[group_g]\n"
                             "[reaction_a]\nreaction = X + Y -> P:Q\nrate = {}\ncutoff = 2\ngroup = g\n"
                             "[reaction_b]\nreaction = X + Y -> R:S\nrate = {}\ncutoff = 2\ngroup = g\n",
                             competition, a, b));
  };
  {
    const Engine e(file("independent", 0.3, 0.4));
    const auto [pa, pb] = frequencies(e, {0, 0}, {1, 0}, trials, rng);
    CHECK_NEAR(pa, 0.3 * 0.6 + 0.3 * 0.4 / 2, tol);  // a passes alone, or both pass and a is drawn
    CHECK_NEAR(pb, 0.4 * 0.7 + 0.3 * 0.4 / 2, tol);
  }
  {
    const Engine e(file("partition", 0.3, 0.4));
    const auto [pa, pb] = frequencies(e, {0, 0}, {1, 0}, trials, rng);
    CHECK_NEAR(pa, 0.3, tol);
    CHECK_NEAR(pb, 0.4, tol);
    const Engine full(file("partition", 0.8, 0.4));  // rescaled: one of them always happens
    const auto [fa, fb] = frequencies(full, {0, 0}, {1, 0}, trials, rng);
    CHECK_NEAR(fa, 0.8 / 1.2, tol);
    CHECK_NEAR(fb, 0.4 / 1.2, tol);
  }
  {
    const Engine e(file("legacy", 0.3, 0.4));
    const auto [pa, pb] = frequencies(e, {0, 0}, {1, 0}, trials, rng);
    CHECK_NEAR(pa, 0.3, tol);
    CHECK_NEAR(pb, 0.7 * 0.7 * 0.7, tol);
    const Engine full(file("legacy", 0.6, 0.7));
    const double p = 0.6 / 1.3;
    const auto [fa, fb] = frequencies(full, {0, 0}, {1, 0}, trials, rng);
    CHECK_NEAR(fa, p, tol);
    CHECK_NEAR(fb, (1 - p) * (1 - p), tol);
  }
}

// ------------------------------------------------------------------------------------------------
// The transition table of the original 3DBFM code: legacy structural units of the center unit a
// and its partner b, the products of the propagation (or of the termination for Term rows), the
// products of the termination branch for Competition rows, and how the pair reacts.
enum class Cls { First, Second, Term, Competition };

struct LegacyRow {
  int a, b;
  int prod_a, prod_b;
  int term_a, term_b;
  Cls cls;
};

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

/// Legacy structural unit of a preset type (PRI and ONE are both unit 1).
int legacy_code(std::string_view type) {
  if (type == "PRI") return 1;
  const auto it = std::ranges::find(kRadicalTypes, type);
  const int i = static_cast<int>(it - kRadicalTypes.begin());
  return i == 0 ? 0 : i - 1;  // MON 0, ONE 1, TWO 2, ..., INI 9
}

/// Preset type of a legacy unit; unit 1 is the chain radical unless `primary`.
std::string_view preset_type(int code, bool primary = false) {
  if (code == 1) return primary ? "PRI" : "ONE";
  return kRadicalTypes[static_cast<std::size_t>(code == 0 ? 0 : code + 1)];
}

struct PairStats {
  int propagations{0}, terminations{0};
  int prod_a{-1}, prod_b{-1}, term_a{-1}, term_b{-1};
};

PairStats sample(const Engine& e, const ReactionFile& f, std::string_view a, std::string_view b, int trials, TestRng& rng) {
  PairStats s;
  const Bead ba{e.type_index(a), 0}, bb{e.type_index(b), 0};
  for (int t = 0; t < trials; ++t) {
    const auto r = e.pair(ba, bb, {1.0, false, false, false}, rng);
    if (!r) continue;
    const bool term = f.reactions[r->reaction].name.starts_with("term:");
    const int pa = legacy_code(e.type_names()[static_cast<std::size_t>(r->center.type)]);
    const int pb = legacy_code(e.type_names()[static_cast<std::size_t>(r->partner.type)]);
    if (term) {
      ++s.terminations;
      if (s.term_a < 0) s.term_a = pa, s.term_b = pb;
    } else {
      ++s.propagations;
      if (s.prod_a < 0) s.prod_a = pa, s.prod_b = pb;
    }
  }
  return s;
}

void test_radical_preset_fidelity() {
  const RadicalRates k{0.0, 0.7, 0.3, 0.5};  // k_p2 + k_term = 0.8 < 1: legacy p = k_p2
  const ReactionFile frp = parse(radical_reaction_file(Mechanism::Frp, k, {}));
  const ReactionFile atrp = parse(radical_reaction_file(Mechanism::Atrp, k, {}));
  CHECK(frp.reactions.size() == 36 && atrp.reactions.size() == 20 && frp.general.competition == Competition::Legacy);
  const Engine ef(frp), ea(atrp);
  TestRng rng(3);
  const int trials = 40000;
  const double tol = 0.015;
  std::set<std::pair<int, int>> in_table;
  for (const LegacyRow& row : kLegacyTable) {
    in_table.insert({row.a, row.b});
    // the legacy unit 1 is a chain radical here; a primary radical gives the same products
    for (const bool primary_a : {false, true}) {
      for (const bool primary_b : {false, true}) {
        if ((primary_a && row.a != 1) || (primary_b && row.b != 1)) continue;
        const auto ta = preset_type(row.a, primary_a), tb = preset_type(row.b, primary_b);
        const PairStats f = sample(ef, frp, ta, tb, trials, rng);
        const PairStats a = sample(ea, atrp, ta, tb, trials, rng);
        const double pf = static_cast<double>(f.propagations) / trials, tf = static_cast<double>(f.terminations) / trials;
        const double pa = static_cast<double>(a.propagations) / trials;
        switch (row.cls) {
          case Cls::First:
          case Cls::Second: {
            const double rate = row.cls == Cls::First ? k.k_p1 : k.k_p2;
            CHECK_NEAR(pf, rate, tol);
            CHECK_NEAR(pa, rate, tol);
            CHECK(f.terminations == 0 && a.terminations == 0);
            CHECK(f.prod_a == row.prod_a && f.prod_b == row.prod_b && a.prod_a == row.prod_a && a.prod_b == row.prod_b);
            break;
          }
          case Cls::Term:
            CHECK(f.propagations == 0 && a.propagations == 0 && a.terminations == 0);
            CHECK_NEAR(tf, k.k_term, tol);
            CHECK(f.term_a == row.prod_a && f.term_b == row.prod_b);
            break;
          case Cls::Competition: {
            const double p = k.k_p2;
            CHECK_NEAR(pf, p, tol);
            CHECK_NEAR(tf, (1 - p) * (1 - p) * std::min(1.0, k.k_term + p), tol);
            CHECK(f.prod_a == row.prod_a && f.prod_b == row.prod_b && f.term_a == row.term_a && f.term_b == row.term_b);
            CHECK_NEAR(pa, p, tol);  // ATRP: a plain pendant propagation
            CHECK(a.terminations == 0 && a.prod_a == row.prod_a && a.prod_b == row.prod_b);
            break;
          }
        }
      }
    }
  }
  // every other pair of units is inert
  for (int ca = 0; ca <= 9; ++ca) {
    for (int cb = 0; cb <= 9; ++cb) {
      if (in_table.contains({ca, cb})) continue;
      const PairStats f = sample(ef, frp, preset_type(ca), preset_type(cb), 2000, rng);
      CHECK(f.propagations == 0 && f.terminations == 0);
    }
  }
}

void test_radical_preset_categories() {
  // Reaction n of the file is the legacy reaction category n (Macromolecules 2025 numbering).
  RadicalRates k{0.0, 1.0, 1.0, 0.0};
  const ReactionFile prop = parse(radical_reaction_file(Mechanism::Frp, k, {}));
  k = {0.0, 0.0, 0.0, 1.0};
  const ReactionFile term = parse(radical_reaction_file(Mechanism::Frp, k, {}));
  const Engine ep(prop), et(term);
  TestRng rng(4);
  const auto category = [&](const Engine& e, std::string_view a, std::string_view b) {
    const auto r = e.pair({e.type_index(a), 0}, {e.type_index(b), 0}, {1.0, false, false, false}, rng);
    return r ? static_cast<int>(r->reaction) + 1 : 0;
  };
  CHECK(category(ep, "PRI", "MON") == 1 && category(ep, "ONE", "MON") == 2 && category(ep, "FOU", "MON") == 3);
  CHECK(category(ep, "SIX", "MON") == 4 && category(ep, "SEV", "MON") == 5 && category(ep, "PRI", "THR") == 6);
  CHECK(category(ep, "ONE", "THR") == 7 && category(ep, "ONE", "FOU") == 9 && category(ep, "PRI", "FIV") == 10);
  CHECK(category(ep, "FOU", "THR") == 12 && category(ep, "FOU", "FOU") == 15 && category(ep, "SIX", "FOU") == 17);
  CHECK(category(ep, "SEV", "FIV") == 20 && category(ep, "MON", "SEV") == 5);
  CHECK(category(et, "PRI", "PRI") == 21 && category(et, "PRI", "ONE") == 22 && category(et, "ONE", "ONE") == 23);
  CHECK(category(et, "FOU", "PRI") == 24 && category(et, "ONE", "FOU") == 25 && category(et, "SIX", "ONE") == 27);
  CHECK(category(et, "SEV", "PRI") == 28 && category(et, "FOU", "FOU") == 30 && category(et, "SIX", "FOU") == 31);
  CHECK(category(et, "SEV", "FOU") == 32 && category(et, "SIX", "SIX") == 33 && category(et, "SEV", "SIX") == 34);
  CHECK(category(et, "SEV", "SEV") == 35);
}

void test_radical_preset_options() {
  RadicalRates k{0.3, 0.2, 0.1, 0.5};
  k.k_px = 0.6;
  PresetOptions o;
  o.cyclization = false;
  o.cutoff = 1.3;
  o.interval = 50;
  o.angle_type = 2;
  const ReactionFile f = parse(radical_reaction_file(Mechanism::Frp, k, o));
  CHECK(f.general.interval == 50 && f.groups[0].bond_type == 1 && f.groups[0].angle_type == 2);
  CHECK(std::ranges::all_of(f.reactions, [](const Reaction& r) { return r.equation.kind == Kind::Unimolecular || (!r.intramolecular && r.cutoff == 1.3); }));
  const Engine e(f);
  TestRng rng(5);
  int x = 0, m = 0, init = 0;
  const int trials = 100000;
  for (int t = 0; t < trials; ++t) {
    x += e.pair({e.type_index("ONE"), 0}, {e.type_index("THR"), 0}, {1.0, false, false, false}, rng) ? 1 : 0;
    m += e.pair({e.type_index("ONE"), 0}, {e.type_index("MON"), 0}, {1.0, false, false, false}, rng) ? 1 : 0;
    init += e.single({e.type_index("INI"), 0}, rng) ? 1 : 0;
  }
  CHECK_NEAR(static_cast<double>(x) / trials, 0.6, 0.006);  // k_px for an unreacted cross-linker
  CHECK_NEAR(static_cast<double>(m) / trials, 0.2, 0.006);
  CHECK_NEAR(static_cast<double>(init) / trials, 0.3, 0.006);
  CHECK(!e.pair({e.type_index("ONE"), 0}, {e.type_index("THR"), 0}, {1.0, false, true, false}, rng));  // no cyclization
  CHECK(throws([] { (void)radical_reaction_file(Mechanism::Atrp, {0, 1.5, 0, 0}, {}); }));
  CHECK(parse_mechanism("FRP") == Mechanism::Frp && !parse_mechanism("raft"));
}

}  // namespace

int main() {
  test_equations();
  test_reaction_file();
  test_engine_matching();
  test_competition();
  test_radical_preset_fidelity();
  test_radical_preset_categories();
  test_radical_preset_options();
  std::cout << g_checks - g_failures << " / " << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
