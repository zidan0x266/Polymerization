// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Self-contained unit tests (no framework dependency). The most important one compares every
// time step with a brute-force implementation of the dynamic lattice liquid rules.

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "dll/io/config.hpp"
#include "dll/model/initial.hpp"
#include "dll/model/system.hpp"
#include "dll/sim/motion.hpp"

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
// Brute force: a site is on a loop if following the attempted moves leads back to it after at
// least three steps; loops that would break a bond are removed in rounds until none does.
struct Outcome {
  std::vector<Bead> occupant;
  std::vector<Vec3i> unwrapped;
};

Outcome reference_step(const DllSystem& sys, const std::vector<std::uint8_t>& field) {
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
  const auto new_pos = [&](Bead b) {
    return moving(b) ? sys.unwrapped(b) + FccGeometry::kNeighbors[field[sys.site(b)]] : sys.unwrapped(b);
  };
  for (;;) {
    std::set<Site> violators;
    for (Site s = 0; s < n; ++s) {
      if (loop_id[s] == n || !accepted.contains(loop_id[s])) continue;
      const Bead b = sys.occupant(s);
      for (const Bead p : sys.bonds(b)) {
        if (lat.minimum_image(new_pos(p) - new_pos(b)).norm2() != 2) violators.insert(loop_id[s]);
      }
    }
    if (violators.empty()) break;
    for (const Site v : violators) accepted.erase(v);
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
  };
  const Case cases[] = {{InitialConfiguration::RandomWalk, 2, 4, 3},  {InitialConfiguration::RandomWalk, 3, 6, 9},
                        {InitialConfiguration::Aligned, 2, 4, 8},     {InitialConfiguration::Aligned, 3, 3, 36},
                        {InitialConfiguration::Aligned, 3, 2, 30},    {InitialConfiguration::RandomWalk, 4, 8, 0}};
  poly::Rng rng(31);
  std::size_t moved = 0, rejected = 0;
  for (const Case& c : cases) {
    DllSystem sys = make_system(c.kind, c.cells, c.chain_length, c.n_chains, 5);
    Motion motion(sys);
    bool same = true;
    for (int step = 0; step < 3000; ++step) {
      const auto field = random_field(sys.size(), rng);
      const Outcome expected = reference_step(sys, field);
      const StepStats stats = motion.step(sys, field);
      moved += stats.moved;
      rejected += stats.rejected;
      for (Site s = 0; s < sys.size(); ++s) same = same && sys.occupant(s) == expected.occupant[s];
      for (Bead b = 0; b < sys.size(); ++b) same = same && sys.unwrapped(b) == expected.unwrapped[b];
      if (!same) break;
    }
    CHECK(same);
    CHECK(sys.validate().has_value());
  }
  CHECK(moved > 0 && rejected > 0);  // both outcomes were exercised
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
  cfg.apply_overrides({"cells=6", "--chain_length=12", "n_chains=40", "run.seed=3", "dump_every=10"});
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
    bad.cells = 2;  // 32 sites for 512 x 16 beads
    bad.validate();
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);
}

}  // namespace

int main() {
  test_initial_configurations();
  test_uniform_field();
  test_against_reference();
  test_invariants();
  test_config();
  std::cout << g_checks - g_failures << " / " << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
