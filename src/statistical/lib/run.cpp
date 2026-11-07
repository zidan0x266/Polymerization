// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "statistical/sim/run.hpp"

#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <optional>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <span>
#include <stdexcept>
#include <string>

#include "poly/chem/psf.hpp"
#include "poly/chem/tables.hpp"
#include "statistical/sim/kinetics.hpp"

namespace statistical {

namespace {

/// Idle events in a row after which no allowed reaction is assumed to remain (every radical
/// left has only partners in its own molecule, or none at all).
constexpr long kMaxIdleEvents = 10'000'000;

class Logger {
 public:
  explicit Logger(std::ostream& out) : out_(out), started_(std::chrono::steady_clock::now()) {}
  void operator()(const std::string& message) const {
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started_).count();
    std::println(out_, "[{:9.1f} s] {}", elapsed, message);
    out_.flush();
  }

 private:
  std::ostream& out_;
  std::chrono::steady_clock::time_point started_;
};

Network run_sample(const Config& config, const poly::chem::Mechanism& mechanism, poly::Rng& rng, const Logger& say, int sample) {
  Network network(poly::chem::initial_composition(config.n_units, config.monovinyl_fraction, config.radical_fraction,
                                                  mechanism.features.explicit_initiation));
  const auto labels = poly::chem::category_labels(mechanism.features.termination);
  poly::chem::ReactionProgress progress;
  progress.reset(labels.size());
  std::size_t initiators = 0, monovinyl = 0, divinyl = 0;
  for (const UnitState s : network.states()) {
    progress.count_unit(s);
    initiators += s.species == poly::chem::Species::Initiator;
    monovinyl += s.species == poly::chem::Species::Monovinyl;
    divinyl += s.species == poly::chem::Species::Divinyl;
  }
  say(std::format("sample {}: {} units ({} initiator-derived, {} monovinyl, {} divinyl); {} vinyl equivalents, conversion {:.6f}",
                  sample, network.size(), initiators, monovinyl, divinyl, progress.total_vinyls, progress.conversion()));

  Kinetics kinetics(network, mechanism.features, config.rates, config.cyclization);
  const std::vector<double> levels = config.snapshot_levels();
  std::size_t next_level = 0;
  while (next_level < levels.size() && progress.conversion() >= levels[next_level]) ++next_level;

  const auto snapshot = [&](double level) {
    poly::chem::SampleTables tables(config.output_dir, sample, labels, "time");
    tables.append(kinetics.time(), progress.conversion(), network.structural_counts(), network.connectivity(),
                  std::span<const std::uint64_t>(progress.category_counts).subspan(1));
    if (config.write_psf) {
      poly::chem::write_psf(
          config.output_dir / std::format("config_{}_{}.psf", sample, poly::chem::legacy_conversion_tag(level)), network.size(),
          [&](std::size_t u) { return poly::chem::legacy_residue_name(network.view(static_cast<Unit>(u)).legacy_code()); },
          [&](std::size_t u) { return network.bonds(static_cast<Unit>(u)); });
    }
    const ConnectivityTally& t = network.connectivity();
    say(std::format("sample {} snapshot at conversion {:.3f} (actual {:.4f}) time={:.6g} radicals={} bonds mm/pm/pp/intra={}/{}/{}/{}",
                    sample, level, progress.conversion(), kinetics.time(), progress.live_radicals, t.monomer_monomer,
                    t.polymer_monomer, t.polymer_polymer, t.intramolecular));
  };
  const auto report = [&](long reactions) {
    say(std::format("sample {} reactions={} time={:.6g} conversion={:.4f} propagations={} terminations={} initiations={} radicals={}",
                    sample, reactions, kinetics.time(), progress.conversion(), progress.propagations, progress.terminations,
                    progress.initiations, progress.live_radicals));
  };

  long reactions = 0, idle = 0;
  std::string ending;
  if (progress.conversion() >= config.stop_conversion - 1e-12) ending = "stop_conversion reached";
  while (ending.empty()) {
    const Event event = kinetics.step(rng, progress);
    if (event == Event::Exhausted) {
      ending = "no reaction possible any more";
      break;
    }
    if (event == Event::Idle) {
      if (++idle >= kMaxIdleEvents) ending = "no allowed reaction left";
      continue;
    }
    idle = 0;
    ++reactions;
    while (next_level < levels.size() && progress.conversion() >= levels[next_level] - 1e-12) snapshot(levels[next_level++]);
    if (progress.conversion() >= config.stop_conversion - 1e-12) ending = "stop_conversion reached";
    if (reactions % config.log_every == 0) report(reactions);
  }
  report(reactions);
  say(std::format("sample {}: {} at conversion {:.6f}", sample, ending, progress.conversion()));
  return network;
}

}  // namespace

Network run(Config config, std::ostream& log) {
  config.validate();
  const Logger say(log);
  config.seed = config.seed ? config.seed : poly::Rng::entropy_seed();  // the dumped configuration then reproduces the run
  std::filesystem::create_directories(config.output_dir);
  {
    std::ofstream dump(config.output_dir / "run_config.ini");
    std::print(dump, "# effective configuration of this run (seed resolved)\n{}", config.to_ini());
  }
  const poly::chem::Mechanism mechanism = poly::chem::find_mechanism(config.mechanism);
  const auto& k = config.rates;
  say(std::format("statistical model, mechanism {}, seed {}", mechanism.name, config.seed));
  say(std::format("rate coefficients: k_init={} k_p1={} k_px={} k_p2={} k_term={}{}", k.k_init, k.k_p1,
                  k.k_px ? std::format("{}", *k.k_px) : std::string("k_p1"), k.k_p2, k.k_term,
                  config.cyclization ? "" : ", no cyclization"));
  poly::Rng rng(config.seed);
  std::optional<Network> last;
  for (int s = 0; s < config.samples; ++s) {
    last.reset();
    last.emplace(run_sample(config, mechanism, rng, say, config.first_sample + s));
  }
  say("all samples finished");
  return std::move(*last);
}

}  // namespace statistical
