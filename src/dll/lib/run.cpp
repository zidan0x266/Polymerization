// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "dll/sim/run.hpp"

#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <stdexcept>

#include "dll/io/writers.hpp"
#include "dll/sim/motion.hpp"

namespace dll {

DllSystem run(Config config, std::ostream& log) {
  config.validate();
  const auto started = std::chrono::steady_clock::now();
  const auto elapsed = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count(); };
  const auto say = [&](const std::string& message) {
    std::println(log, "[{:9.1f} s] {}", elapsed(), message);
    log.flush();
  };

  config.seed = config.seed ? config.seed : poly::Rng::entropy_seed();  // the dumped configuration then reproduces the run
  std::filesystem::create_directories(config.output_dir);
  {
    std::ofstream dump(config.output_dir / "run_config.ini");
    std::print(dump, "# effective configuration of this run (seed resolved)\n{}", config.to_ini());
  }

  poly::Rng rng(config.seed);
  DllSystem system(config.cells, config.chain_length,
                   initial_chains(config.initial, config.cells, config.chain_length, config.n_chains, rng));
  if (auto ok = system.validate(); !ok) throw std::logic_error("initial configuration check failed: " + ok.error());
  say(std::format("dynamic lattice liquid: {} chains of {} beads in {} sites ({}^3 FCC cells), polymer fraction {:.4f}, seed {}",
                  system.n_chains(), system.chain_length(), system.size(), config.cells,
                  static_cast<double>(system.chain_beads()) / static_cast<double>(system.size()), config.seed));
  say(std::format("starting configuration: {}", to_string(config.initial)));

  std::ofstream trajectory;
  if (config.dump_every > 0) {
    trajectory.open(config.output_dir / "trajectory.lammpstrj", std::ios::binary);
    if (!trajectory) throw std::runtime_error("cannot write the trajectory file");
    write_dump_frame(trajectory, system, 0, config.dump_solvent);
  }

  Motion motion(system);
  StepStats interval;
  double interval_started = elapsed();
  for (long step = 1; step <= config.steps; ++step) {
    const StepStats s = motion.step(system, rng);
    interval.loops += s.loops;
    interval.rejected += s.rejected;
    interval.moved += s.moved;
    if (config.dump_every > 0 && step % config.dump_every == 0) write_dump_frame(trajectory, system, step, config.dump_solvent);
    if (step % config.log_every == 0 || step == config.steps) {
      const long done = step % config.log_every == 0 ? config.log_every : step % config.log_every;
      const double per_step = static_cast<double>(done);
      const double seconds = elapsed() - interval_started;
      say(std::format("step {}: {:.4f} of the beads moved per step in {:.1f} loops; {:.1f} % of the loops broke a bond; "
                      "{:.1f} ns per site and step",
                      step, static_cast<double>(interval.moved) / per_step / static_cast<double>(system.size()),
                      static_cast<double>(interval.loops) / per_step,
                      interval.loops + interval.rejected
                          ? 100.0 * static_cast<double>(interval.rejected) / static_cast<double>(interval.loops + interval.rejected)
                          : 0.0,
                      1e9 * seconds / per_step / static_cast<double>(system.size())));
      interval = {};
      interval_started = elapsed();
    }
  }

  if (auto ok = system.validate(); !ok) throw std::logic_error("final configuration check failed: " + ok.error());
  say(std::format("finished {} steps; final configuration checked", config.steps));
  return system;
}

}  // namespace dll
