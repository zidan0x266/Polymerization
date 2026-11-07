# 3dbfm — reactive Monte Carlo on the 3D bond fluctuation model

Monomers, cross-linkers and initiators move on a cubic lattice under excluded volume and bond
constraints and react with neighboring units according to the chosen polymerization
mechanism: free radical (FRP) or atom transfer radical (ATRP) polymerization. Everything is set
at run time from a small configuration file. Build instructions are in the [repository README](../../README.md).

## Run

```bash
build/src/3dbfm/bfm_react examples/example_1/atrp_paper.ini
build/src/3dbfm/bfm_react examples/example_1/frp_paper.ini seed=12345 samples=1   # override keys
build/src/3dbfm/bfm_react --print-config                                          # all keys with defaults
```

Each run writes `run_config.ini`, the full configuration including the random seed, so any run
can be repeated exactly. Example inputs are described in
[`examples/example_1`](../../examples/example_1/README.md).

## Configuration

Keys can be grouped under optional `[system]`, `[scheme]`, `[run]` and `[output]` headers.

| | Key | Default | Meaning |
|:--|:--|:--|:--|
| **System** | `lattice` | `100 100 100` | box size `lx ly lz` (even integers from 4 to 32766) |
| | `n_units` | `106250` | total number of units (initiators + monomers + cross-linkers) |
| | `monovinyl_fraction` | `0.96` | fraction of mono-vinyl units, initiators included; the rest are cross-linkers |
| | `radical_fraction` | `0.002` | fraction of initiator units |
| **Scheme** | `mechanism` | `atrp` | `frp` or `atrp` |
| | `k_init` | `0` | initiation probability |
| | `k_p1`, `k_p2` | `0.005`, `0.0025` | propagation probability onto a free / pendant vinyl |
| | `k_term` | `0` | termination probability |
| | `competition_rule` | `legacy` | see *Model notes* |
| | `initiation_requires_neighbor` | `true` | see *Model notes* |
| **Run** | `samples`, `first_sample` | `5`, `1` | number of replicas and number of the first one |
| | `seed` | `0` | random seed; `0` picks a fresh one |
| | `mcs_initial_relaxation`, `mcs_equilibration` | `1000`, `10000` | Monte Carlo steps before the reaction |
| | `max_mcs` | `990000000` | upper limit on reaction steps |
| | `stop_conversion` | `0.96` | vinyl conversion at which the run ends |
| | `snapshot_start`, `snapshot_step` | `0.03`, `0.01` | conversions at which output is written |
| | `extra_snapshots` | — | additional conversions, comma separated |
| | `log_every_mcs` | `1000` | progress report interval |
| | `initial_configuration` | — | restart file to start from |
| **Output** | `directory` | `.` | output directory |
| | `write_pdb`, `write_psf`, `write_restart` | `true` | which files to write |

## Output (per sample `s`)

- `structural_s.dat`: number of units of each structural type
- `propagation_s.dat`: new bonds classified as monomer–monomer, polymer–monomer, polymer–polymer or intramolecular
- `reactype_s.dat`: counts of each reaction type
- `config_s_<conversion>.pdb/.psf`: configuration and bonds for visualization and network analysis
- `equilibrated_s.dat`, `final_s.dat`: restart files

## Model notes

The defaults reproduce the implementation used in the publications listed in the repository
README. Two choices can be changed:

- `competition_rule`: how a pair that could either propagate or terminate is resolved. `legacy`
  (default) keeps the published behavior; `partition` uses a single draw with probabilities
  proportional to `k_p2` and `k_term`.
- `initiation_requires_neighbor`: with `true` (default) an initiator decomposes only when a
  neighboring unit is sampled; `false` makes initiation independent of neighbors.

Reaction attempts in which neither unit is reactive (no radical, no dormant initiator) cannot
succeed, so they are skipped. All other attempts keep their original probability, and results
are statistically unchanged. A sample also stops early once no reactive unit is left.

## Performance

Time per unit and Monte Carlo step for the paper-size system (100³ lattice, 106,250 units,
rate constants of the 2025 paper), on one core of an Apple M4 with Apple Clang 21. The original
program used in the publications was compiled with the same optimization flags.

| phase | original program | this code | with profile-guided optimization |
|---|---|---|---|
| relaxation (moves only) | 16.2 ns | 11.6 ns | 10.7 ns |
| reaction, ATRP | 26.8 ns | 14.7 ns | 13.4 ns |
| reaction, FRP | 27.2 ns | 15.1 ns | 13.3 ns |

This code is 1.4 times faster than the original program for moves and 1.8 times faster for
reactions. One reaction step of this system takes about 1.6 ms, so a million steps take about
26 minutes. Output writing is not included.

## Code structure

```
include/bfm/
  core/      basic types, random numbers and index sets
  lattice/   bond-vector set and periodic lattice
  model/     unit states and the simulated system
  reaction/  reaction schemes (FRP, ATRP)
  io/        configuration and output
  sim/       simulation driver
lib/         implementation
apps/        bfm_react
```

Unit tests are in [`tests/3dbfm`](../../tests/3dbfm). A new polymerization mechanism is added as a
new `ReactionScheme`, which also declares which units are reactive. The lattice and system
classes are independent of the chemistry.
