// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "bfm/sim/simulation.hpp"

#include <format>
#include <fstream>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>

#include "bfm/io/writers.hpp"

namespace bfm {

Simulation::Simulation(Config config, std::unique_ptr<ReactionScheme> scheme, std::ostream& log)
    : config_(std::move(config)), scheme_(std::move(scheme)), log_(log) {
  config_.validate();
  const std::uint64_t seed = config_.seed ? config_.seed : Rng::entropy_seed();
  config_.seed = seed;  // the dumped configuration then reproduces the run
  rng_ = std::make_unique<Rng>(seed);
  started_ = std::chrono::steady_clock::now();
}

void Simulation::log_line(const std::string& message) const {
  const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started_).count();
  std::println(log_, "[{:9.1f} s] {}", elapsed, message);
  log_.flush();
}

void Simulation::run() {
  std::filesystem::create_directories(config_.output_dir);
  {
    std::ofstream dump(config_.output_dir / "run_config.ini");
    std::print(dump, "# effective configuration of this run (seed resolved)\n{}", config_.to_ini());
  }
  log_line(std::format("3DBFM reactive Monte Carlo, mechanism {}, seed {}", scheme_->name(), config_.seed));
  const auto& k = scheme_->rate_constants();
  log_line(std::format("rate constants: k_init={} k_p1={} k_p2={} k_term={}, competition rule {}", k.k_init, k.k_p1, k.k_p2,
                       k.k_term, to_string(config_.competition_rule)));
  for (int s = 0; s < config_.samples; ++s) run_sample(config_.first_sample + s);
  log_line("all samples finished");
}

void Simulation::prepare_system(int sample) {
  system_ = std::make_unique<BfmSystem>(PeriodicLattice(config_.lattice_x, config_.lattice_y, config_.lattice_z), bonds_);
  const auto composition = scheme_->initial_composition(config_.n_units, config_.monovinyl_fraction, config_.radical_fraction);

  progress_.reset(scheme_->category_labels().size());
  std::size_t n_init = 0;
  for (const UnitState& s : composition) {
    if (s.species == Species::Initiator) ++n_init;
  }

  if (!config_.initial_configuration.empty()) {
    log_line(std::format("sample {}: loading configuration {}", sample, config_.initial_configuration.string()));
    auto restart = read_legacy_configuration(config_.initial_configuration, system_->lattice().dimensions(), config_.n_units);
    system_->load(std::move(restart.units), restart.unwrapped);
    for (const Monomer& m : system_->monomers()) progress_.count_unit(m.state);
    return;
  }
  for (const UnitState& s : composition) progress_.count_unit(s);

  system_->place_randomly(composition, *rng_);
  const auto n_mono_total = static_cast<std::size_t>(static_cast<double>(config_.n_units) * config_.monovinyl_fraction);
  log_line(std::format("sample {}: placed {} units ({} initiator-derived, {} monovinyl, {} divinyl) on a {}x{}x{} lattice; "
                       "{} vinyl equivalents",
                       sample, config_.n_units, n_init, n_mono_total - n_init, config_.n_units - n_mono_total, config_.lattice_x,
                       config_.lattice_y, config_.lattice_z, progress_.total_vinyls));
  relax(config_.mcs_initial_relaxation, "initial relaxation");
  relax(config_.mcs_equilibration, "equilibration");
  if (config_.write_restart) {
    write_legacy_configuration(config_.output_dir / std::format("equilibrated_{}.dat", sample), *system_);
  }
}

void Simulation::relax(long n_mcs, const char* label) {
  const auto n = static_cast<std::uint32_t>(system_->size());
  for (long mcs = 0; mcs < n_mcs; ++mcs) {
    for (std::uint32_t i = 0; i < n; ++i) {
      const auto [unit, dir] = rng_->below2(n, kNumDirections);
      system_->attempt_move(static_cast<Index>(unit), static_cast<Direction>(dir));
    }
  }
  log_line(std::format("{}: {} MCS of athermal moves done", label, n_mcs));
}

std::optional<ReactionAttempt> map_reaction_draw(const BfmSystem& system, const IndexSet& reactive, std::uint32_t x,
                                                 int slot) noexcept {
  const auto n = static_cast<std::uint32_t>(system.size());
  const auto na = static_cast<std::uint32_t>(reactive.size());
  if (2 * static_cast<std::uint64_t>(na) > n) return ReactionAttempt{static_cast<Index>(x), slot};
  if (x < na) return ReactionAttempt{reactive[x], slot};
  if (x < 2 * na) {
    const Index center = system.partner_at(reactive[x - na], slot);
    if (center == kNone || reactive.contains(center)) return std::nullopt;
    return ReactionAttempt{center, BondVectorSet::opposite_offset(slot)};
  }
  return std::nullopt;
}

void Simulation::rebuild_reactive_set() {
  reactive_.reset(static_cast<std::size_t>(system_->size()));
  for (Index id = 0; id < system_->size(); ++id) update_reactive(id);
}

void Simulation::update_reactive(Index id) {
  const Monomer& m = system_->monomer(id);
  reactive_.set(id, scheme_->is_reactive(UnitView{m.state, m.nbonds}));
}

bool Simulation::attempt_reaction(Index center, int slot) {
  const Index partner = system_->partner_at(center, slot);
  const Monomer& mc = system_->monomer(center);
  const UnitView vc{mc.state, mc.nbonds};

  if (partner == kNone) {
    if (scheme_->has_unimolecular_step() && !config_.initiation_requires_neighbor) {
      if (const auto r = scheme_->propose_unimolecular(vc, *rng_)) {
        apply(*r, center, kNone, slot);
        return true;
      }
    }
    return false;
  }

  if (mc.bonded_to(partner)) return false;  // no double bonds between the same pair

  const Monomer& mp = system_->monomer(partner);
  const UnitView vp{mp.state, mp.nbonds};
  if (const auto r = scheme_->propose_pair(vc, vp, *rng_)) {
    apply(*r, center, partner, slot);
    return true;
  }
  if (scheme_->has_unimolecular_step()) {
    // Legacy quirk: the initiator branch is only reached when a partner was found.
    if (const auto r = scheme_->propose_unimolecular(vc, *rng_)) {
      apply(*r, center, kNone, slot);
      return true;
    }
  }
  return false;
}

void Simulation::apply(const Reaction& reaction, Index center, Index partner, int slot) {
  if (reaction.kind == EventKind::Initiation) {
    poly::chem::apply(reaction, system_->monomer(center).state, progress_);
  } else {
    poly::chem::apply(reaction, system_->monomer(center).state, system_->monomer(partner).state, progress_);
    system_->bond_through_offset(center, partner, slot);
  }
  // The attacking unit first: the order of the reactive set decides which unit a draw picks.
  const bool partner_first = reaction.kind == EventKind::Propagation && reaction.second_attacks;
  update_reactive(partner_first ? partner : center);
  if (partner != kNone) update_reactive(partner_first ? center : partner);
}

void Simulation::snapshot(int sample, long mcs, double level) {
  SampleTables tables(config_.output_dir, sample, scheme_->category_labels());
  tables.append(mcs, level, system_->structural_counts(), system_->connectivity(),
                std::span<const std::uint64_t>(progress_.category_counts).subspan(1));
  const std::string tag = legacy_conversion_tag(level);
  if (config_.write_pdb) write_pdb(config_.output_dir / std::format("config_{}_{}.pdb", sample, tag), *system_);
  if (config_.write_psf) write_psf(config_.output_dir / std::format("config_{}_{}.psf", sample, tag), *system_);
  const ConnectivityTally& t = system_->connectivity();
  log_line(std::format("sample {} snapshot at conversion {:.3f} (actual {:.3f}) mcs={} radicals={} bonds mm/pm/pp/intra={}/{}/{}/{}",
                       sample, level, progress_.conversion(), mcs, progress_.live_radicals, t.monomer_monomer, t.polymer_monomer,
                       t.polymer_polymer, t.intramolecular));
}

void Simulation::react(int sample) {
  const auto n = static_cast<std::uint32_t>(system_->size());
  const std::vector<double> levels = config_.snapshot_levels();
  std::size_t next_level = 0;
  while (next_level < levels.size() && progress_.conversion() >= levels[next_level]) ++next_level;
  system_->reset_connectivity_tally();
  rebuild_reactive_set();

  log_line(std::format("sample {}: reaction phase starts at conversion {:.6f} with {} reactive units", sample,
                       progress_.conversion(), reactive_.size()));
  long mcs = 0;
  bool done = false;
  bool exhausted = reactive_.empty();
  for (; mcs < config_.max_mcs && !done && !exhausted; ++mcs) {
    for (std::uint32_t i = 0; i < n; ++i) {
      const auto [unit, dir] = rng_->below2(n, kNumDirections);
      system_->attempt_move(static_cast<Index>(unit), static_cast<Direction>(dir));
      const auto [x, slot] = rng_->below2(n, kNumReactionOffsets);
      const auto attempt = map_reaction_draw(*system_, reactive_, x, static_cast<int>(slot));
      if (attempt && attempt_reaction(attempt->center, attempt->slot)) {
        while (next_level < levels.size() && progress_.conversion() >= levels[next_level] - 1e-12) {
          snapshot(sample, mcs, levels[next_level]);
          ++next_level;
        }
        if (progress_.conversion() >= config_.stop_conversion - 1e-12) {
          done = true;
          break;
        }
        if (reactive_.empty()) {  // no radical or dormant initiator left: nothing can react any more
          exhausted = true;
          break;
        }
      }
    }
    if (!done && !exhausted && (mcs + 1) % config_.log_every_mcs == 0) {
      log_line(std::format("sample {} mcs={} conversion={:.4f} propagations={} terminations={} initiations={} radicals={}", sample,
                           mcs + 1, progress_.conversion(), progress_.propagations, progress_.terminations, progress_.initiations,
                           progress_.live_radicals));
    }
  }
  if (exhausted) {
    log_line(std::format("sample {}: no reactive units left, stopping before stop_conversion", sample));
  } else if (!done) {
    log_line(std::format("sample {}: max_mcs reached before stop_conversion", sample));
  }
  if (config_.write_restart) {
    write_legacy_configuration(config_.output_dir / std::format("final_{}.dat", sample), *system_);
  }
  log_line(std::format("sample {}: finished after {} MCS at conversion {:.6f}", sample, mcs, progress_.conversion()));
}

void Simulation::run_sample(int sample) {
  std::filesystem::create_directories(config_.output_dir);
  prepare_system(sample);
  react(sample);
}

}  // namespace bfm
