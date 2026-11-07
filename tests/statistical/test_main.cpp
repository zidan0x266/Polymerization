// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
//
// Self-contained unit tests (no framework dependency) of the statistical model: exact kinetics
// of the Gillespie algorithm, the choice of events and partners, and the bookkeeping.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "poly/chem/progress.hpp"
#include "poly/chem/scheme.hpp"
#include "statistical/io/config.hpp"
#include "statistical/model/network.hpp"
#include "statistical/sim/kinetics.hpp"
#include "statistical/sim/run.hpp"

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

using namespace statistical;
using poly::chem::RadicalChainScheme;
using poly::chem::RateConstants;
using poly::chem::ReactionProgress;
using poly::chem::Species;

struct Setup {
  Network network;
  ReactionProgress progress;
  Kinetics kinetics;

  Setup(std::string_view mechanism, std::size_t n, double monovinyl, double radical, RateConstants k, bool cyclization = true)
      : network(poly::chem::initial_composition(n, monovinyl, radical, poly::chem::find_mechanism(mechanism).features.explicit_initiation)),
        kinetics(network, poly::chem::find_mechanism(mechanism).features, k, cyclization) {
    progress.reset(poly::chem::category_labels(poly::chem::find_mechanism(mechanism).features.termination).size());
    for (const UnitState s : network.states()) progress.count_unit(s);
  }
  /// Run until the conversion reaches `conversion`; returns false if nothing (allowed) can react
  /// before.
  bool run_to(double conversion, poly::Rng& rng) {
    long idle = 0;
    while (progress.conversion() < conversion) {
      const Event e = kinetics.step(rng, progress);
      if (e == Event::Exhausted || (e == Event::Idle && ++idle > 1'000'000)) return false;
      if (e == Event::Reaction) idle = 0;
    }
    return true;
  }
};

void test_network() {
  const UnitState mono = UnitState::monovinyl_monomer(), di = UnitState::divinyl_monomer(), rad = UnitState::primary_radical();
  Network net({rad, mono, di, di, mono});
  CHECK(net.validate().has_value() && net.molecule_radicals(0) == 1 && net.molecule_pendants(2) == 0);
  net.set_state(0, {Species::Initiator, 0, 0});
  net.set_state(2, {Species::Divinyl, 1, 1});
  net.bond(0, 2);  // initiator fragment adds to a cross-linker, which keeps a pendant vinyl
  CHECK(net.molecule_radicals(0) == 1 && net.molecule_pendants(0) == 1 && net.molecule(0) == net.molecule(2));
  net.bond(3, 4);
  net.bond(2, 3);
  net.bond(4, 0);  // closes a ring
  const ConnectivityTally& t = net.connectivity();
  CHECK(t.monomer_monomer == 2 && t.polymer_polymer == 1 && t.intramolecular == 1 && t.polymer_monomer == 0);
  CHECK(net.bonded(4, 0) && !net.bonded(1, 0) && net.molecule(1) == 1);
  CHECK(net.validate().has_value());
  bool threw = false;
  try {
    net.bond(0, 2);
  } catch (const std::logic_error&) {
    threw = true;
  }
  CHECK(threw);
}

void test_exact_kinetics() {
  // Living polymerization of a monovinyl monomer: R radicals stay constant, so the monomer
  // decays as exp(-k_p1 R t / N).
  {
    poly::Rng rng(1);
    Setup s("atrp", 200000, 1.0, 0.01, RateConstants{0, 0.3, 0, 0});
    const double m0 = static_cast<double>(s.kinetics.monomers());
    CHECK(s.run_to(0.6, rng));
    const double rate = std::log(m0 / static_cast<double>(s.kinetics.monomers())) / s.kinetics.time();
    CHECK(std::fabs(rate / (0.3 * 2000.0 / 200000.0) - 1.0) < 0.01);
  }
  // Initiator decomposition alone: first order with k_init.
  {
    poly::Rng rng(2);
    Setup s("frp", 100000, 1.0, 0.2, RateConstants{0.05, 0, 0, 0});
    const std::size_t d0 = s.kinetics.dormant();
    while (2 * s.kinetics.dormant() > d0) s.kinetics.step(rng, s.progress);
    CHECK(std::fabs(s.kinetics.time() / (std::log(2.0) / 0.05) - 1.0) < 0.02);
    CHECK(s.progress.initiations == d0 - s.kinetics.dormant() && s.kinetics.radicals() == s.progress.initiations);
  }
}

void test_event_choice() {
  // Early in a copolymerization, cross-linker and monomer additions occur in the ratio
  // k_px X / (k_p1 M); k_px unset means k_p1.
  for (const bool set : {true, false}) {
    poly::Rng rng(set ? 3 : 4);
    RateConstants k{0, 0.2, 0.1, 0};
    if (set) k.k_px = 0.6;
    Setup s("atrp", 400000, 0.8, 0.005, k);
    const double m = static_cast<double>(s.kinetics.monomers()), x = static_cast<double>(s.kinetics.crosslinkers());
    CHECK(s.run_to(0.03, rng));
    const auto& c = s.progress.category_counts;  // one-based legacy categories
    const double onto_monomer = static_cast<double>(c[1] + c[2] + c[3] + c[4] + c[5]);
    const double onto_crosslinker = static_cast<double>(c[6] + c[7] + c[12] + c[13] + c[14]);
    const double expected = (set ? 0.6 : 0.2) * x / (0.2 * m);
    CHECK(std::fabs(onto_crosslinker / onto_monomer / expected - 1.0) < 0.06);
  }
}

void test_living_chains() {
  // Ideal living polymerization: the number of monomers added to a chain is Poisson distributed.
  poly::Rng rng(5);
  Setup s("atrp", 400000, 1.0, 0.005, RateConstants{0, 1.0, 0, 0});
  CHECK(s.run_to(0.5, rng));
  std::map<Unit, double> size;
  for (Unit u = 0; u < s.network.size(); ++u) {
    if (!s.network.bonds(u).empty()) ++size[s.network.molecule(u)];
  }
  double sum = 0, sum2 = 0;
  for (const auto& [root, n] : size) {
    sum += n - 1;
    sum2 += (n - 1) * (n - 1);
  }
  const double count = static_cast<double>(size.size()), mean = sum / count;
  CHECK(size.size() == 2000);
  CHECK(std::fabs((sum2 / count - mean * mean) / mean - 1.0) < 0.1);
}

void test_cyclization() {
  // Without cyclization no bond forms within a molecule, so every molecule stays a tree.
  for (const bool cyclization : {false, true}) {
    poly::Rng rng(cyclization ? 6 : 7);
    RateConstants k{0, 0.1, 0.1, 0};
    k.k_px = 0.2;
    Setup s("atrp", 100000, 0.1, 0.1, k, cyclization);
    s.run_to(0.6, rng);
    std::size_t bonds = 0, molecules = 0;
    for (Unit u = 0; u < s.network.size(); ++u) {
      bonds += s.network.bonds(u).size();
      molecules += s.network.molecule(u) == u ? 1 : 0;
    }
    bonds /= 2;
    if (cyclization) {
      CHECK(s.network.connectivity().intramolecular > 0 && bonds > s.network.size() - molecules);
    } else {
      CHECK(s.network.connectivity().intramolecular == 0 && bonds == s.network.size() - molecules);
    }
    CHECK(s.network.validate().has_value() && s.kinetics.validate().has_value());
  }
}

void test_bookkeeping() {
  for (const bool cyclization : {true, false}) {
    poly::Rng rng(cyclization ? 8 : 9);
    RateConstants k{0.002, 0.5, 0.25, 0.9};
    k.k_px = 0.8;
    Setup s("frp", 20000, 0.9, 0.02, k, cyclization);
    bool ok = true;
    for (int e = 0; e < 40000 && ok; ++e) {
      if (s.kinetics.step(rng, s.progress) == Event::Exhausted) break;
      if (e % 2000 != 0) continue;
      ok = s.network.validate().has_value() && s.kinetics.validate().has_value();
      std::uint64_t radicals = 0, consumed = 0, bonds = 0;
      for (Unit u = 0; u < s.network.size(); ++u) {
        const UnitState st = s.network.state(u);
        radicals += st.radicals;
        consumed += static_cast<std::uint64_t>(poly::chem::max_vinyls(st.species) - st.vinyls) + (st.species == Species::Initiator ? 1 : 0);
        bonds += s.network.bonds(u).size();
      }
      ok = ok && radicals == s.progress.live_radicals && consumed == s.progress.consumed_vinyls &&
           bonds == 2 * (s.progress.propagations + s.progress.terminations) && radicals + 2 * s.progress.terminations == s.progress.initiations;
    }
    CHECK(ok);
    CHECK(s.progress.terminations > 0 && s.progress.conversion() > 0.3);
  }
}

void test_config_and_driver() {
  Config cfg;
  CHECK(!cfg.set("chain_length", "4").has_value());
  cfg.apply_overrides({"n_units=20000", "mechanism=FRP", "k_init=0.001", "k_p1=5", "k_px=10", "k_p2=2.5", "k_term=50",
                       "radical_fraction=0.01", "cyclization=no", "stop_conversion=0.4", "snapshot_start=0.1", "snapshot_step=0.1",
                       "samples=2", "seed=3"});
  cfg.validate();  // rate coefficients may exceed one
  CHECK(cfg.mechanism == "frp" && cfg.rates.k_px == 10.0 && !cfg.cyclization);
  const auto dir = std::filesystem::temp_directory_path() / "stat_test_driver";
  std::filesystem::remove_all(dir);
  cfg.output_dir = dir;
  {
    std::ofstream out(dir.string() + ".ini");
    out << cfg.to_ini();
  }
  CHECK(Config::from_file(dir.string() + ".ini").to_ini() == cfg.to_ini());
  std::ostringstream log;
  const Network last = run(cfg, log);
  CHECK(last.validate().has_value() && log.str().contains("stop_conversion reached"));
  for (const char* name : {"run_config.ini", "structural_2.dat", "propagation_1.dat", "reactype_2.dat", "config_1_0.200000.psf",
                           "config_2_0.400000.psf"}) {
    CHECK(std::filesystem::exists(dir / name));
  }
  std::filesystem::remove_all(dir);
  std::filesystem::remove(dir.string() + ".ini");
  bool threw = false;
  try {
    Config bad;
    bad.rates.k_p2 = -1;
    bad.validate();
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);
}

}  // namespace

int main() {
  test_network();
  test_exact_kinetics();
  test_event_choice();
  test_living_chains();
  test_cyclization();
  test_bookkeeping();
  test_config_and_driver();
  std::cout << g_checks - g_failures << " / " << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
