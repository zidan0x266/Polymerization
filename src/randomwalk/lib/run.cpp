// SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
// Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>

#include "rw/sim/run.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <ostream>  // std::print / std::println overloads for std::ostream
#include <print>
#include <stdexcept>

#include "rw/io/writers.hpp"
#include "rw/sim/packer.hpp"

namespace rw {

Packing run(Config config, std::ostream& log) {
  config.validate();
  const auto started = std::chrono::steady_clock::now();
  const auto say = [&](const std::string& message) {
    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    std::println(log, "[{:9.3f} s] {}", elapsed, message);
    log.flush();
  };

  config.seed = config.seed ? config.seed : poly::Rng::entropy_seed();  // the dumped configuration then reproduces the run
  std::filesystem::create_directories(config.output_dir);
  {
    std::ofstream dump(config.output_dir / "run_config.ini");
    std::print(dump, "# effective configuration of this run (seed resolved)\n{}", config.to_ini());
  }

  const int cells = config.box_cells();
  const std::size_t n_beads = config.n_chains * config.chain_length + config.n_crosslinkers;
  const double sites = 4.0 * cells * cells * cells;
  const double box = std::sqrt(2.0) * cells * config.bond_length;
  say(std::format("random-walk packing: {} chains of {} beads and {} cross-linkers, seed {}", config.n_chains,
                  config.chain_length, config.n_crosslinkers, config.seed));
  say(std::format("FCC lattice of {}^3 cells, box edge {:.4f}, {:.0f} sites; occupancy {:.4f}, number density {:.4f}", cells,
                  box, sites, static_cast<double>(n_beads) / sites, static_cast<double>(n_beads) / (box * box * box)));

  poly::Rng rng(config.seed);
  Packer packer(cells, config.walk, rng);
  Packing packing = packer.pack(config.n_chains, config.chain_length, config.n_crosslinkers, config.bond_length);
  if (auto ok = validate(packing, config.walk.angles); !ok) throw std::logic_error("packing check failed: " + ok.error());
  say(std::format("packed and checked; {} trapped chain attempts were grown again", packer.dead_ends()));

  const std::filesystem::path data = config.output_dir / config.data_file;
  write_lammps_data(data, packing, config.mass,
                    std::format("LAMMPS data file: {} chains of {} beads and {} cross-linkers, FCC random walks "
                                "(Polymerization randomwalk, seed {})",
                                config.n_chains, config.chain_length, config.n_crosslinkers, config.seed));
  say(std::format("wrote {}", data.string()));
  if (config.write_xyz) {
    const std::filesystem::path xyz = std::filesystem::path(data).replace_extension(".xyz");
    write_xyz(xyz, packing);
    say(std::format("wrote {}", xyz.string()));
  }
  return packing;
}

}  // namespace rw
