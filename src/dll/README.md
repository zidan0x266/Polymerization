# dll — dynamic lattice liquid

The dynamic lattice liquid (DLL) model of Pakula fills every site of a lattice with a bead:
initiators, monomers, cross-linkers, polymer units and solvent. A bead can move only together
with its neighbors, along closed loops in which each bead takes the place of the next.
Polanowski, Jeszka and co-workers (Lodz University of Technology, Poland) used it to simulate
polymerization and gelation [1–4]. This module implements the model from the published
descriptions; no original code was available. The chemistry, free radical (FRP) and atom
transfer radical (ATRP) polymerization of monomers and divinyl cross-linkers, is shared with the
3dbfm module.

**Status — on-going.** Phase 1, the motion of chains and solvent, and phase 2, the reactions, are
done. The published descriptions leave some details open (see [Model](#model)); they are being
checked against the published DLL results, in particular the gel points of ref. 4. Faster,
parallel time steps are planned.

## Run

```bash
build/src/dll/dll_react examples/example_3/atrp_gel.ini                  # ATRP gelation in bulk
build/src/dll/dll_react examples/example_3/atrp_gel.ini seed=12345 k_p1=0.01  # override keys
build/src/dll/dll_react examples/example_3/solution.ini                  # chains in solvent, no reactions
build/src/dll/dll_react --print-config                                    # all keys with defaults
```

Each run writes `run_config.ini`, the full configuration including the random seed, so any run
can be repeated exactly.

## Configuration

Keys can be grouped under optional `[system]`, `[scheme]`, `[run]` and `[output]` headers. The
composition, chemistry and snapshot keys mean the same as in the
[3dbfm module](../3dbfm/README.md).

| | Key | Default | Meaning |
|:--|:--|:--|:--|
| **System** | `cells` | `32` | cubic FCC cells per box edge (2 to 1000); the box has 4 × cells³ sites |
| | `n_units` | `131072` | reactive units: initiators, monomers and cross-linkers; sites left over hold solvent |
| | `monovinyl_fraction` | `0.96` | fraction of the units that are initiators or monovinyl monomers; the rest are divinyl cross-linkers |
| | `radical_fraction` | `0.002` | fraction of the units that are initiators |
| | `n_chains` | `0` | inert precursor chains, placed before the units |
| | `chain_length` | `16` | beads per precursor chain |
| | `initial` | `random_walk` | precursor chains: `random_walk` (self-avoiding walks, up to a chain fraction of about 0.7) or `aligned` (straight chains along lattice rows, up to a melt) |
| | `allow_crossing` | `false` | let chains pass through each other, for comparison only (see [Model](#model)) |
| **Scheme** | `mechanism` | `atrp` | `frp` or `atrp` |
| | `k_init` | `0` | initiator decomposition per time step (FRP) |
| | `k_p1`, `k_p2` | `0.005`, `0.0025` | propagation onto a monomer or cross-linker vinyl, onto a pendant vinyl |
| | `k_term` | `0` | termination by recombination (FRP) |
| | `competition_rule` | `legacy` | propagation against termination for one pair, as in 3dbfm |
| | `reaction_partner` | `random_neighbor` | which neighbors are tried in a time step: `random_neighbor` or `all_neighbors` (see [Model](#model)) |
| **Run** | `samples`, `first_sample` | `1`, `1` | number of independent samples and the number of the first |
| | `seed` | `0` | random seed; `0` picks a fresh one |
| | `equilibration_steps` | `0` | time steps of motion before the reactions start |
| | `steps` | `1000000000` | most time steps of the reaction phase; a system without units moves for this many steps |
| | `stop_conversion` | `0.96` | end of the run |
| | `snapshot_start`, `snapshot_step` | `0.03`, `0.01` | conversions at which snapshots are written |
| | `extra_snapshots` | — | further snapshot conversions, comma separated |
| | `log_every` | `1000` | progress report interval in time steps |
| | `initial_configuration` | — | configuration file to start from (`equilibrated_*.dat`, `final_*.dat`); the system keys are then ignored |
| **Output** | `directory` | `.` | output directory |
| | `write_pdb`, `write_psf` | `true` | PDB and PSF snapshots |
| | `write_restart` | `true` | configuration files after the equilibration and at the end |
| | `dump_every` | `0` | write a trajectory frame every this many steps; `0` writes none |
| | `dump_solvent` | `false` | include the solvent in the trajectory |

As in 3dbfm, the initiators are counted among the monovinyl units: a composition
[I]₀ : [X]₀ : [M]₀ = 1 : 2 : 100 is `radical_fraction = 1/103` and `monovinyl_fraction = 101/103`.
The conversion is that of 3dbfm too: vinyl groups consumed, with every initiator counted as one
consumed vinyl.

## Output

For sample *s*:

- `structural_s.dat`, `propagation_s.dat`, `reactype_s.dat`: the tables of the 3dbfm module
  (structural units, connectivity of new bonds, reaction categories), one line per snapshot. The
  first column is the time step, the second the conversion reached at that step.
- `config_s_<conversion>.pdb` and `.psf`: snapshots in the 3dbfm layout, with the units and
  precursor-chain beads (residue `CHN`) but not the solvent; lengths in bond lengths. A snapshot
  is written at the end of the time step in which the conversion passes its level.
- `equilibrated_s.dat`, `final_s.dat`: configuration files (unwrapped positions, unit states and
  bonds) from which a run can continue with `initial_configuration`.
- `trajectory_s.lammpstrj`: frames in the LAMMPS dump format (read by OVITO, VMD and most analysis
  tools) with unwrapped positions in bond lengths. `mol` identifies the molecule (0 for solvent);
  `type` is 1 for precursor chains, 2 for solvent and 3 plus the legacy structural code (0–9) for
  units.
- `run_config.ini`: the configuration of the run.
- The log reports the conversion and reaction counts, the fraction of beads that moved per step,
  the shares of loops removed because they would break a bond or let chains cross, and the time
  per site and step.

## Model

The lattice is face-centered cubic with twelve neighbors per site; bonds join neighboring sites,
so their length is fixed. One time step consists of:

1. Every bead picks one of its twelve neighbors at random as the direction of an attempted move.
2. Because all sites are occupied, a bead can move only if the bead ahead moves out of the way.
   Following the attempted moves, the successful ones are the closed loops of at least three
   beads; two beads trying to swap places form no loop.
3. A loop is removed if it would leave a chain bond longer than one lattice step, taking into
   account which neighbors move as well. Removing a loop can make a bond in another loop break,
   so removal is repeated in rounds until no remaining loop breaks a bond.
4. A loop is removed if it would let chains pass through each other (below); since this changes
   the moves around it, steps 3 and 4 are repeated until no loop fails.
5. All remaining loops move at once, each bead into the site of the next.
6. The reaction stage follows (below). Motion and reaction stages alternate, so their order
   within a step does not matter.

Every site holds one bead before and after a step, so the excluded volume is always kept. The
cooperative moves can still carry a bond through another bond within one step, without any
two beads meeting, which would let a chain cut through another chain or through itself and
change its entanglements. To decide this, every moving bead is taken to travel in a straight
line to its new site, all at the same speed, and a bond is the straight segment between its
beads. A loop is removed in step 4 if, at some instant, two bonds would touch, or a bead without
bonds (solvent or a free monomer) would touch a bond; all loops taking part are removed. The test
is exact (`sim/crossing.hpp`): segments can only meet when their four ends are coplanar, which
happens at the real roots of a cubic in time. Only bonds whose two beads both move can take part
in a contact, so only those are checked, against the moving bonds and beads around them. In
solutions of 16- and 64-bead chains, the rule removes about 2 × 10⁻⁵ of the loops attempted; it
hardly changes the mobility, but over long runs it keeps the chain topology.
`allow_crossing = true` switches it off, for comparison only.

In the reaction stage, units on neighboring sites react with the probabilities given by the rate
constants, by the rules of the 3dbfm module: a radical adds across the vinyl group of a monomer
(`k_p1`) or a pendant vinyl (`k_p2`), two radicals recombine (`k_term`, FRP), and a dormant
initiator decomposes (`k_init` per step, FRP). A new bond always joins neighbors, and the motion
keeps it intact, so branched molecules and networks move without further rules. Only units that
carry a radical, or dormant initiators, can start a reaction; those reactive at the start of
the stage take part, in random order, and a radical passed on during the stage acts from the
next step. `reaction_partner` sets which pairs are tried:

- `random_neighbor`: every unit is the center of one attempt per step, with one of its twelve
  neighbors drawn at random. The rate constants are then probabilities per attempt, as in 3dbfm,
  and a given pair of neighbors is tried twice per step with probability 1/12 each.
- `all_neighbors`: every pair of neighbors is tried once per step, so the rate constants are
  probabilities per pair and time step.

The published descriptions leave some details open, settled here as follows:

- Bond check: loops are removed in simultaneous rounds (step 3), so the outcome does not depend
  on the order in which loops are examined.
- Chain crossing: later papers of the Lodz group state that loops that would make chains cross
  are immobilized [5], but not how this is decided in 3D; the straight-line rule above is ours.
  Loops that break a bond are removed first, so that a move that does not happen cannot block
  another.
- Reaction attempts: the papers give a reaction probability per time step [5] but not whether an
  active unit tries one neighbor or all of them; both rules are available.

Checks of the implementation, motion:

- A chain in solvent samples its conformations uniformly, as it should without interactions. For
  chains of 4 and 6 beads, the distributions of end-to-end distance and bond angle agree with
  exact enumeration of all self-avoiding walks within statistical error, and so does the mean
  square end-to-end distance of dilute 8-bead chains (9.49 ± 0.09 against 9.556 b²). Longer
  chains grow with the Flory exponent of a good solvent, ν ≈ 0.59.
- Solvent beads diffuse normally: the mean square displacement grows by 0.055 b² per step from
  10 to 10,000 steps, equal to the fraction of beads that move per step, so successive moves are
  uncorrelated.
- Dense systems move much more slowly than solutions: about 6 % of the beads move per step in
  pure solvent, 1 % at a polymer fraction of 0.5, and 0.1 % in a relaxed melt of tetramers. A
  melt started from straight chains of 8 or more beads hardly moves at all, so dense systems
  are best started from a random-walk configuration where possible.

- Crossing: the contact test agrees with a brute-force search over time and both segments for
  439,299 random lattice configurations, and no contact occurs unless a bond moves at both ends.
  Every time step is compared with a brute-force implementation that tests every pair of bonds
  in every periodic image. With crossing forbidden, dilute 8-bead chains still give the exact
  mean square end-to-end distance (9.43 ± 0.09 against 9.556 b²).

Checks of the reactions:

- Unit tests compare the probability that a pair reacts in one step with its exact value for
  both partner rules (radical and monomer, two radicals, initiator decomposition), and check the
  bookkeeping at every step of moving, reacting ATRP and FRP systems: radicals, conversion,
  bonds and reaction counts. A living polymerization of a monovinyl monomer gives one linear
  chain per initiator.
- Living polymerization in bulk, [M]₀/[I]₀ = 49, `k_p1 = 0.002`, 32,000 sites: at conversion 0.5
  the chain lengths have a variance-to-mean ratio of 1.15, close to the Poisson value 1 of ideal
  living polymerization; it rises to 1.6 at conversion 0.8, when the polymer slows down the
  monomers.
- Gel points of ATRP copolymerization in bulk, [I]₀ : [X]₀ : [M]₀ = 1 : x : 100,
  `k_p1 = k_p2 = 0.02`, 32,000 sites, from the maximum of the reduced degree of polymerization:
  vinyl conversion 0.54 for x = 5 and 0.84 for x = 2, and no gel up to 0.90 for x = 1, against
  0.32, 0.50 and 0.71 from Flory–Stockmayer theory. Cross-links within the same molecule delay
  gelation, the effect that brought the DLL gel points of ref. 2 close to experiment.

## Performance

One time step costs about 10 ns per lattice site on one core of an Apple M4 (Apple Clang 21),
independent of the composition: 0.17 ms for 16,384 sites, 1.4 ms for 131,072 sites and 13 ms for
about a million sites. The crossing check adds about 4 % and the reaction stage little (12 ns per
site and step in bulk ATRP with 1 % initiator): `atrp_gel.ini` reaches conversion 0.9 in about
14,000 steps and 24 s. The motion
is local and the same in every part of the box, so it parallelizes naturally; this is planned for
a later phase.

## Code structure

```
include/dll/
  model/   the fully occupied lattice with its units and molecules, and the starting configurations
  sim/     the motion and reaction stages of a time step, the crossing test, and the driver
  io/      configuration, snapshots, trajectories and configuration files
lib/       implementation
apps/      dll_react
```

Unit tests are in [`tests/dll`](../../tests/dll); every time step of the motion is compared there
with a brute-force implementation of the rules above. The chemistry (unit states, reaction
schemes, bookkeeping and tables), the FCC geometry, the random number generator and the
configuration reader are shared with the other modules (`src/common`), and the random-walk
starting configurations come from the randomwalk module.

## References

1. P. Polanowski, T. Pakula. Studies of polymer conformation and dynamics in two dimensions using
   simulations based on the Dynamic Lattice Liquid (DLL) model. *J. Chem. Phys.* **2002**, 117,
   4022–4029. [10.1063/1.1495836](https://doi.org/10.1063/1.1495836)
2. P. Polanowski, J. K. Jeszka, K. Matyjaszewski. Modeling of branching and gelation in living
   copolymerization of monomer and divinyl cross-linker using dynamic lattice liquid model (DLL)
   and Flory–Stockmayer model. *Polymer* **2010**, 51, 6084–6092.
   [10.1016/j.polymer.2010.10.006](https://doi.org/10.1016/j.polymer.2010.10.006)
3. P. Polanowski, J. K. Jeszka, W. Li, K. Matyjaszewski. Effect of dilution on branching and
   gelation in living copolymerization of monomer and divinyl cross-linker: Modeling using
   dynamic lattice liquid model (DLL) and Flory–Stockmayer (FS) model. *Polymer* **2011**, 52,
   5092–5101. [10.1016/j.polymer.2011.08.055](https://doi.org/10.1016/j.polymer.2011.08.055)
4. J. Lyu, Y. Gao, Z. Zhang, U. Greiser, P. Polanowski, J. K. Jeszka, K. Matyjaszewski, H. Tai,
   W. Wang. Monte Carlo Simulations of Atom Transfer Radical (Homo)polymerization of Divinyl
   Monomers: Applicability of Flory–Stockmayer Theory. *Macromolecules* **2018**, 51, 6673–6681.
   [10.1021/acs.macromol.8b01630](https://doi.org/10.1021/acs.macromol.8b01630)
5. K. Hałagan, M. Banaszak, J. Jung, P. Polanowski, A. Sikorski. Computer simulation of the
   dynamic behavior of double polymer brush–solvent systems. *Comput. Methods Sci. Technol.*
   **2021**, 27, 141–149. [10.12921/cmst.2021.0000030](https://doi.org/10.12921/cmst.2021.0000030)
