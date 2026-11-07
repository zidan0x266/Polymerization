// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/sim/run.hpp"

#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <memory>
#include <optional>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <span>
#include <stdexcept>
#include <string>

#include "dll/io/writers.hpp"
#include "dll/sim/motion.hpp"
#include "dll/sim/reaction.hpp"
#include "poly/chem/tables.hpp"

namespace dll {

namespace {

class Logger {
 public:
  explicit Logger(std::ostream& out) : out_(out), started_(std::chrono::steady_clock::now()) {}
  [[nodiscard]] double elapsed() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - started_).count();
  }
  void operator()(const std::string& message) const {
    std::println(out_, "[{:9.1f} s] {}", elapsed(), message);
    out_.flush();
  }

 private:
  std::ostream& out_;
  std::chrono::steady_clock::time_point started_;
};

/// Motion statistics accumulated between two log lines.
class MotionReport {
 public:
  void add(const StepStats& s) noexcept {
    stats_.loops += s.loops;
    stats_.rejected += s.rejected;
    stats_.blocked += s.blocked;
    stats_.moved += s.moved;
    ++steps_;
  }
  /// Summary since the last call, then start a new interval.
  std::string take(const DllSystem& system, double now) {
    const double steps = steps_ ? static_cast<double>(steps_) : 1.0;
    const double tried = static_cast<double>(stats_.loops + stats_.rejected + stats_.blocked);
    const auto share = [tried](std::size_t k) { return tried > 0 ? 100.0 * static_cast<double>(k) / tried : 0.0; };
    std::string text = std::format(
        "{:.4f} of the beads moved per step in {:.1f} loops; {:.1f} % of the loops broke a bond, {:.2f} % would cross a "
        "chain; {:.1f} ns per site and step",
        static_cast<double>(stats_.moved) / steps / static_cast<double>(system.size()), static_cast<double>(stats_.loops) / steps,
        share(stats_.rejected), share(stats_.blocked), 1e9 * (now - since_) / steps / static_cast<double>(system.size()));
    stats_ = {};
    steps_ = 0;
    since_ = now;
    return text;
  }
  void restart(double now) noexcept {
    stats_ = {};
    steps_ = 0;
    since_ = now;
  }

 private:
  StepStats stats_{};
  long steps_{0};
  double since_{0.0};
};

class SampleRun {
 public:
  SampleRun(const Config& config, const poly::chem::ReactionScheme& scheme, poly::Rng& rng, const Logger& say, int sample)
      : config_(config), scheme_(scheme), rng_(rng), say_(say), sample_(sample), system_(prepare()) {}

  DllSystem run() {
    equilibrate();
    react();
    if (auto ok = system_.validate(); !ok) throw std::logic_error("final configuration check failed: " + ok.error());
    if (config_.write_restart) write_configuration(config_.output_dir / std::format("final_{}.dat", sample_), system_);
    return std::move(system_);
  }

 private:
  DllSystem prepare() {
    if (!config_.initial_configuration.empty()) {
      say_(std::format("sample {}: loading configuration {}", sample_, config_.initial_configuration.string()));
      return read_configuration(config_.initial_configuration);
    }
    const auto chains = initial_chains(config_.initial, config_.cells, config_.chain_length, config_.n_chains, rng_);
    const auto units = scheme_.initial_composition(config_.n_units, config_.monovinyl_fraction, config_.radical_fraction);
    const auto points = random_free_points(config_.cells, chains, units.size(), rng_);
    loaded_ = false;
    return DllSystem(config_.cells, config_.chain_length, chains, points, units);
  }

  void describe() {
    if (auto ok = system_.validate(); !ok) throw std::logic_error("initial configuration check failed: " + ok.error());
    std::size_t initiators = 0, monovinyl = 0, divinyl = 0;
    progress_.reset(scheme_.category_labels().size());
    for (const UnitState s : system_.states()) {
      progress_.count_unit(s);
      initiators += s.species == poly::chem::Species::Initiator;
      monovinyl += s.species == poly::chem::Species::Monovinyl;
      divinyl += s.species == poly::chem::Species::Divinyl;
    }
    const std::size_t solvent = system_.size() - system_.solute_beads();
    say_(std::format("sample {}: {} sites ({}^3 FCC cells): {} units ({} initiator-derived, {} monovinyl, {} divinyl), "
                     "{} chains of {} beads, {} solvent; {} vinyl equivalents, conversion {:.6f}",
                     sample_, system_.size(), system_.lattice().cells(), system_.n_units(), initiators, monovinyl, divinyl,
                     system_.n_chains(), system_.chain_length(), solvent, progress_.total_vinyls, progress_.conversion()));
  }

  void equilibrate() {
    describe();
    if (config_.equilibration_steps > 0) {
      report_.restart(say_.elapsed());
      for (long step = 1; step <= config_.equilibration_steps; ++step) {
        report_.add(motion_.step(system_, rng_));
        if (step % config_.log_every == 0 || step == config_.equilibration_steps) {
          say_(std::format("sample {} equilibration step {}: {}", sample_, step, report_.take(system_, say_.elapsed())));
        }
      }
    }
    if (config_.write_restart && !loaded_) {
      write_configuration(config_.output_dir / std::format("equilibrated_{}.dat", sample_), system_);
    }
  }

  void react() {
    const std::vector<double> levels = config_.snapshot_levels();
    std::size_t next_level = 0;
    while (next_level < levels.size() && progress_.conversion() >= levels[next_level]) ++next_level;
    system_.reset_connectivity_tally();
    Reactor reactor(system_, scheme_, config_.reaction_partner, config_.cyclization);
    const bool reactive = system_.n_units() > 0;
    if (reactive) {
      say_(std::format("sample {}: reaction phase starts at conversion {:.6f} with {} reactive units", sample_,
                       progress_.conversion(), reactor.reactive_count()));
    }

    std::ofstream trajectory;
    if (config_.dump_every > 0) {
      trajectory.open(config_.output_dir / std::format("trajectory_{}.lammpstrj", sample_), std::ios::binary);
      if (!trajectory) throw std::runtime_error("cannot write the trajectory file");
      write_dump_frame(trajectory, system_, 0, config_.dump_solvent);
    }

    report_.restart(say_.elapsed());
    long step = 0;
    std::string ending;
    if (reactive && progress_.conversion() >= config_.stop_conversion - 1e-12) {
      ending = "stop_conversion reached";
    } else if (reactive && reactor.reactive_count() == 0) {
      ending = "no reactive units";
    }
    while (ending.empty() && step < config_.steps) {
      ++step;
      report_.add(motion_.step(system_, rng_));
      reactor.step(system_, rng_, progress_);
      if (config_.dump_every > 0 && step % config_.dump_every == 0) write_dump_frame(trajectory, system_, step, config_.dump_solvent);
      while (next_level < levels.size() && progress_.conversion() >= levels[next_level] - 1e-12) {
        snapshot(step, levels[next_level]);
        ++next_level;
      }
      if (reactive && progress_.conversion() >= config_.stop_conversion - 1e-12) {
        ending = "stop_conversion reached";
      } else if (reactive && reactor.reactive_count() == 0) {
        ending = "no reactive units left";
      }
      if (step % config_.log_every == 0 || !ending.empty() || step == config_.steps) log_progress(step);
    }
    if (ending.empty()) ending = "step limit reached";
    say_(std::format("sample {}: {} after {} steps at conversion {:.6f}", sample_, ending, step, progress_.conversion()));
  }

  void log_progress(long step) {
    const std::string motion = report_.take(system_, say_.elapsed());
    if (system_.n_units() == 0) {
      say_(std::format("sample {} step {}: {}", sample_, step, motion));
      return;
    }
    say_(std::format("sample {} step {}: conversion={:.4f} propagations={} terminations={} initiations={} radicals={}; {}", sample_,
                     step, progress_.conversion(), progress_.propagations, progress_.terminations, progress_.initiations,
                     progress_.live_radicals, motion));
  }

  void snapshot(long step, double level) {
    poly::chem::SampleTables tables(config_.output_dir, sample_, scheme_.category_labels(), "step");
    tables.append(step, progress_.conversion(), system_.structural_counts(), system_.connectivity(),
                  std::span<const std::uint64_t>(progress_.category_counts).subspan(1));
    const std::string tag = poly::chem::legacy_conversion_tag(level);
    if (config_.write_pdb) write_pdb(config_.output_dir / std::format("config_{}_{}.pdb", sample_, tag), system_);
    if (config_.write_psf) write_psf(config_.output_dir / std::format("config_{}_{}.psf", sample_, tag), system_);
    const ConnectivityTally& t = system_.connectivity();
    say_(std::format("sample {} snapshot at conversion {:.3f} (actual {:.4f}) step={} radicals={} bonds mm/pm/pp/intra={}/{}/{}/{}",
                     sample_, level, progress_.conversion(), step, progress_.live_radicals, t.monomer_monomer, t.polymer_monomer,
                     t.polymer_polymer, t.intramolecular));
  }

  const Config& config_;
  const poly::chem::ReactionScheme& scheme_;
  poly::Rng& rng_;
  const Logger& say_;
  int sample_;
  bool loaded_{true};
  DllSystem system_;
  Motion motion_{system_, config_.allow_crossing};
  poly::chem::ReactionProgress progress_{};
  MotionReport report_{};
};

}  // namespace

DllSystem run(Config config, std::ostream& log) {
  config.validate();
  const Logger say(log);
  config.seed = config.seed ? config.seed : poly::Rng::entropy_seed();  // the dumped configuration then reproduces the run
  std::filesystem::create_directories(config.output_dir);
  {
    std::ofstream dump(config.output_dir / "run_config.ini");
    std::print(dump, "# effective configuration of this run (seed resolved)\n{}", config.to_ini());
  }
  const auto scheme = poly::chem::make_scheme(config.mechanism, config.rates, config.competition_rule);
  const auto& k = scheme->rate_constants();
  say(std::format("dynamic lattice liquid, mechanism {}, seed {}", scheme->name(), config.seed));
  say(std::format("rate constants: k_init={} k_p1={} k_p2={} k_term={}, competition rule {}, reaction partner {}{}", k.k_init,
                  k.k_p1, k.k_p2, k.k_term, poly::chem::to_string(config.competition_rule), to_string(config.reaction_partner),
                  config.cyclization ? "" : ", no cyclization"));

  poly::Rng rng(config.seed);
  std::optional<DllSystem> last;
  for (int s = 0; s < config.samples; ++s) {
    last.reset();  // free the previous sample before building the next
    last.emplace(SampleRun(config, *scheme, rng, say, config.first_sample + s).run());
  }
  say("all samples finished");
  return std::move(*last);
}

}  // namespace dll
