# statistical — coordinate-free kinetic Monte Carlo

The statistical model follows a radical polymerization reaction by reaction without any spatial
information: every unit can meet every other one, as in a perfectly mixed system. It is the
kinetic Monte Carlo (Gillespie) counterpart of the lattice models, with the same chemistry: free
radical (FRP) and atom transfer radical (ATRP) polymerization of monovinyl monomers and divinyl
cross-linkers, the same unit states, reaction categories and output tables, and the network
topology in PSF files. It is the model of the "w.c." and "wo.c." simulations of refs. 2 and 3.

## Run

```bash
build/src/statistical/stat_react examples/example_4/divinyl_2018.ini
build/src/statistical/stat_react examples/example_4/divinyl_2018.ini cyclization=false   # "wo.c." model
build/src/statistical/stat_react --print-config                                        # all keys with defaults
```

Each run writes `run_config.ini`, the full configuration including the random seed, so any run
can be repeated exactly.

## Configuration

Keys can be grouped under optional `[system]`, `[scheme]`, `[run]` and `[output]` headers. The
composition and snapshot keys mean the same as in the [3dbfm module](../3dbfm/README.md).

| | Key | Default | Meaning |
|:--|:--|:--|:--|
| **System** | `n_units` | `1000000` | number of units (2 to 4 × 10⁹) |
| | `monovinyl_fraction` | `0.96` | fraction of the units that are initiators or monovinyl monomers; the rest are divinyl cross-linkers |
| | `radical_fraction` | `0.002` | fraction of the units that are initiators |
| **Scheme** | `mechanism` | `atrp` | `frp` or `atrp` |
| | `k_init` | `0` | initiator decomposition rate (FRP) |
| | `k_p1`, `k_p2` | `0.005`, `0.0025` | propagation onto an unreacted monomer, onto a pendant vinyl |
| | `k_px` | not set | propagation onto an unreacted cross-linker; if not set, `k_p1` |
| | `k_term` | `0` | termination by recombination (FRP) |
| | `cyclization` | `true` | reactions within a molecule (rings, intramolecular cross-links) take place like any other; `false` forbids them |
| **Run** | `samples`, `first_sample` | `1`, `1` | number of independent samples and the number of the first |
| | `seed` | `0` | random seed; `0` picks a fresh one |
| | `stop_conversion` | `0.96` | end of the run |
| | `snapshot_start`, `snapshot_step` | `0.03`, `0.01` | conversions at which snapshots are written |
| | `extra_snapshots` | — | further snapshot conversions, comma separated |
| | `log_every` | `1000000` | progress report interval in reactions |
| **Output** | `directory` | `.` | output directory |
| | `write_psf` | `true` | PSF snapshots |

## Output

For sample *s*:

- `structural_s.dat`, `propagation_s.dat`, `reactype_s.dat`: the tables of the lattice modules
  (structural units, connectivity of new bonds, reaction categories), one line per snapshot.
  `structural` and `reactype` start with the time and the conversion, `propagation` with the
  conversion.
- `config_s_<conversion>.psf`: the network topology at each snapshot conversion, in the layout of
  the lattice modules (legacy residue names, every bond listed from both ends). There are no
  coordinates, so no PDB files are written.
- `run_config.ini`: the configuration of the run.
- The log reports the time, conversion, reaction counts and radicals.

## Model

The system holds N units in a volume of N, so the concentration of a species is its number divided
by N. Each elementary reaction has a rate:

| Reaction | Rate |
|:--|:--|
| initiator decomposition (FRP) | `k_init` × dormant initiators |
| propagation onto a monomer | `k_p1` × radicals × unreacted monomers / N |
| propagation onto a cross-linker | `k_px` × radicals × unreacted cross-linkers / N |
| propagation onto a pendant vinyl | `k_p2` × radicals × pendant vinyls / N |
| termination by recombination (FRP) | `k_term` × radicals × (radicals − 1) / 2N |

Each step of Gillespie's direct method [1] advances the time by an exponentially distributed
interval with the total rate as its mean frequency, chooses a reaction with probability
proportional to its rate, and then draws its radical, and its partner, at random. A cross-linker
carries two vinyl groups, so equal reactivity per vinyl group, the convention of refs. 2 and 3, is
`k_px = 2 k_p1` with `k_p2 = k_p1`. Only the ratios of the rate constants affect the structure; their
size sets the unit of time.

As in ref. 3, ATRP is treated as an ideal living polymerization: all initiators carry a radical
from the start, there is no termination, and the fast exchange between active and dormant chain
ends, which changes only the time scale, is not followed explicitly.

A radical may react with a pendant vinyl of its own molecule, with the probability that the drawn
pendant vinyl lies on that molecule ("w.c." model of ref. 3). With `cyclization = false`, partners
from the radical's own molecule are drawn again, so every reaction joins two molecules and the
network stays free of cycles, as Flory–Stockmayer theory assumes ("wo.c." model). In both cases a
unit never reacts with itself or with a unit it is already bonded to. If no allowed partner is left
for any radical, the run ends.

Checks of the implementation:

- Unit tests verify the exact kinetics (first-order monomer consumption in living
  polymerization and initiator decomposition), the choice between monomer and cross-linker
  additions, Poisson chain lengths in living polymerization, the absence of rings without
  cyclization, and the bookkeeping of radicals, vinyl groups, bonds and molecules.
- Gel points of ATRP homopolymerization of a divinyl monomer, from the maximum of the reduced
  degree of polymerization, 10⁶ units, equal reactivity per vinyl group (`divinyl_2018.ini`),
  against the simulations of ref. 3 and Flory–Stockmayer theory, √([I]₀/2[X]₀). Vinyl conversion
  at the gel point, cross-linker conversion in parentheses:

  | [X]₀/[I]₀ | with cyclization | ref. 3 (w.c.) | without cyclization | ref. 3 (wo.c.) | Flory–Stockmayer |
  |:--|:--|:--|:--|:--|:--|
  | 100/1 | 0.073 (0.141) | 0.07 (0.14) | 0.075 (0.145) | 0.07 (0.14) | 0.071 |
  | 100/10 | 0.228 (0.405) | 0.24 (0.42) | 0.223 (0.396) | 0.23 (0.41) | 0.224 |

  Without cyclization the gel point is the Flory–Stockmayer value, as it must be, and all values
  agree with ref. 3 to within 0.01 (between seeds they vary by about ±0.003). Ref. 3 drew a
  radical's partner by primary chain (a chain first, then a pendant vinyl on it); drawing the
  pendant vinyl directly, as here, gives the same gel points.
- Free radical polymerization of a monovinyl monomer (10⁷ units, `k_init = 10⁻⁵`, `k_p1 = 1`,
  `k_term = 1000`, 1 % initiator): the number of radicals follows the steady state
  √(`k_init` × initiators × N / `k_term`) within 1–3 %, every termination makes one dead chain,
  and the dead chains have the number-average length 2(ν + 1) of recombination, with ν the
  kinetic chain length, and a dispersity of 1.46 (1.5 for long chains).

## Performance

Every reaction takes a constant time: 10⁶ units of an ATRP copolymerization reach conversion 0.96
in 0.3 s on one core of an Apple M4. Writing the PSF snapshots takes most of the time of a run.

## Code structure

```
include/statistical/
  model/   the network: units, bonds and molecules
  sim/     the kinetic Monte Carlo and the driver
  io/      configuration and snapshots
lib/       implementation
apps/      stat_react
```

Unit tests are in [`tests/statistical`](../../tests/statistical). The chemistry (unit states,
mechanisms, reaction categories, bookkeeping, tables and PSF output) and the random number
generator are shared with the lattice modules (`src/common`).

## References

1. D. T. Gillespie. Exact stochastic simulation of coupled chemical reactions. *J. Phys. Chem.*
   **1977**, 81, 2340–2361. [10.1021/j100540a008](https://doi.org/10.1021/j100540a008)
2. H. Gao, P. Polanowski, K. Matyjaszewski. Gelation in living copolymerization of monomer and
   divinyl cross-linker: comparison of ATRP experiments with Monte Carlo simulations.
   *Macromolecules* **2009**, 42, 5925–5932. [10.1021/ma901005d](https://doi.org/10.1021/ma901005d)
3. J. Lyu, Y. Gao, Z. Zhang, U. Greiser, P. Polanowski, J. K. Jeszka, K. Matyjaszewski, H. Tai,
   W. Wang. Monte Carlo Simulations of Atom Transfer Radical (Homo)polymerization of Divinyl
   Monomers: Applicability of Flory–Stockmayer Theory. *Macromolecules* **2018**, 51, 6673–6681.
   [10.1021/acs.macromol.8b01630](https://doi.org/10.1021/acs.macromol.8b01630)
