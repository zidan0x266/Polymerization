# dll — dynamic lattice liquid

The dynamic lattice liquid (DLL) model of Pakula fills every site of a lattice with a bead:
chain units, solvent and, later, monomers. A bead can move only together with its neighbors,
along closed loops in which each bead takes the place of the next. Polanowski, Jeszka and
co-workers (Lodz University of Technology, Poland) used it to simulate polymerization and
gelation [1–4]. This module implements the model from the published descriptions; no original
code was available.

**Status — on-going.** Phase 1, the motion of chains and solvent, is done. Reactions (free
radical and atom transfer radical polymerization, sharing the chemistry of the 3dbfm module) are
in preparation.

## Run

```bash
build/src/dll/dll_react examples/example_3/solution.ini
build/src/dll/dll_react examples/example_3/solution.ini steps=20000 seed=12345   # override keys
build/src/dll/dll_react --print-config                                          # all keys with defaults
```

Each run writes `run_config.ini`, the full configuration including the random seed, so any run
can be repeated exactly.

## Configuration

Keys can be grouped under optional `[system]`, `[run]` and `[output]` headers.

| | Key | Default | Meaning |
|:--|:--|:--|:--|
| **System** | `cells` | `16` | cubic FCC cells per box edge (2 to 1000); the box has 4 × cells³ sites |
| | `chain_length` | `16` | beads per chain |
| | `n_chains` | `512` | number of chains; all other sites hold solvent |
| | `initial` | `random_walk` | starting chains: `random_walk` (self-avoiding walks, up to a polymer fraction of about 0.7) or `aligned` (straight chains along lattice rows, up to a melt) |
| **Run** | `seed` | `0` | random seed; `0` picks a fresh one |
| | `steps` | `10000` | number of time steps |
| | `log_every` | `1000` | progress report interval |
| **Output** | `directory` | `.` | output directory |
| | `dump_every` | `0` | write a trajectory frame every this many steps; `0` writes none |
| | `dump_solvent` | `false` | include the solvent in the trajectory |

## Output

- `trajectory.lammpstrj`: frames in the LAMMPS dump format (read by OVITO, VMD and most analysis
  tools) with unwrapped positions in units of the bond length; chains have their molecule
  number and type 1, solvent has molecule 0 and type 2.
- `run_config.ini`: the configuration of the run.
- The log reports the fraction of beads that moved per step, the number of loops, and the share
  of loops removed because they would break a bond.

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
4. All remaining loops move at once, each bead into the site of the next.

The published descriptions leave two details open, settled here as follows:

- Bond check: loops are removed in simultaneous rounds (step 3), so the outcome does not depend
  on the order in which loops are examined.
- Bonds passing through each other during a cooperative move are not prevented. This matters
  only for entangled dynamics; the 2D rule of ref. 1 (solvent jumping through a bond) has no
  direct 3D counterpart on the FCC lattice.

Checks of the implementation:

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

## Performance

One time step costs about 10 ns per lattice site on one core of an Apple M4 (Apple Clang 21):
0.17 ms for the 16,384 sites of `solution.ini`, 1.4 ms for 131,072 sites and 13 ms for about a
million sites, independent of the polymer fraction. The motion is local and the same in every
part of the box, so it parallelizes naturally; this is planned for a later phase.

## Code structure

```
include/dll/
  model/   the fully occupied lattice and the starting configurations
  sim/     the DLL time step and the driver
  io/      configuration and trajectory output
lib/       implementation
apps/      dll_react
```

Unit tests are in [`tests/dll`](../../tests/dll); every time step is compared there with a
brute-force implementation of the rules above. The FCC geometry, random number generator and
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
