# randomwalk — chain packing for end-linked networks

Precursor chains are grown as self-avoiding random walks on a face-centered cubic (FCC) lattice
and single-bead cross-linkers are put on random empty sites. The result is written as a LAMMPS
data file, the starting point of coarse-grained end-linking simulations with
`fix bond/create`. The original program, the *Amorphous Polymer Generator*, was obtained from
Mark A. Tschopp and Don K. Ward, Center for Advanced Vehicular Systems, Mississippi State
University, who used it in [1, 2]; the end-linking protocol follows Riggleman and co-workers
[3–5]. This module built the monodisperse, end-linked networks of the 2026 fracture paper listed
in the [repository README](../../README.md).

## Run

```bash
build/src/randomwalk/rw_pack examples/example_2/endlink_paper.ini
build/src/randomwalk/rw_pack examples/example_2/endlink_paper.ini chain_length=8 n_chains=12416 n_crosslinkers=6208
build/src/randomwalk/rw_pack --print-config        # all keys with defaults
```

Each run writes `run_config.ini`, the full configuration including the random seed, so any run
can be repeated exactly. The example inputs and the LAMMPS scripts that relax, compress and
end-link the packed system are described in [`examples/example_2`](../../examples/example_2/README.md).

## Configuration

Keys can be grouped under optional `[system]`, `[walk]`, `[run]` and `[output]` headers.

| | Key | Default | Meaning |
|:--|:--|:--|:--|
| **System** | `chain_length` | `18` | beads per precursor chain (at least 2); both end beads are reactive |
| | `n_chains` | `5680` | number of precursor chains |
| | `n_crosslinkers` | `2840` | number of single-bead cross-linkers |
| | `occupancy` | `0.05` | sets the box: the smallest number of cells in which the chain beads fill less than this fraction of the lattice sites |
| | `cells` | `0` | cubic cells per box edge (2 to 1000); `0` derives it from `occupancy` |
| | `bond_length` | `1.0` | bond length b, the nearest-neighbor distance of the lattice |
| | `mass` | `1.0` | mass of every atom type |
| **Walk** | `bond_angles` | `60, 90, 120` | bond angles a chain may form, from 60, 90, 120 and 180 degrees |
| | `site_weighting` | `legacy` | `legacy`, `free_neighbors` or `uniform`; see *Model notes* |
| | `max_attempts` | `1000` | growth attempts per chain before the run stops with an error |
| **Run** | `seed` | `0` | random seed; `0` picks a fresh one |
| **Output** | `directory` | `.` | output directory |
| | `data_file` | `config.lmp` | name of the LAMMPS data file |
| | `write_xyz` | `true` | also write an XYZ file for visualization |

## Output

- `config.lmp`: LAMMPS data file, atom style `molecular`, with molecule IDs and image flags.
  Chains come first (chain *c* holds atoms *c*·*L* + 1 to (*c* + 1)·*L*), then the cross-linkers,
  each a molecule of its own. Atom types carry the labels `bead` (1, chain interior), `end`
  (2, reactive chain end), `xlink` (3, cross-linker) and `xlink_full` (4, cross-linker with all
  four sites reacted, assigned during end-linking). Bonds have type 1.
- `config.xyz`: extended XYZ file with the atom type and the wrapped position of every atom.
- `run_config.ini`: the configuration of the run.

## Model

Sites of the FCC lattice are one bond length apart, so every bond has length b and no two beads
are closer than b. The box edge is √2·b times the number of cells. Each chain starts on a random
empty site; every further bead goes to an empty neighbor of the last one that forms an allowed
bond angle with the previous bond, chosen according to `site_weighting`. A chain whose last bead
has no such neighbor is removed and grown again from a new random site. The cross-linkers take
random empty sites after all chains are in place.

The defaults reproduce the program used in the publications, including two of its quirks:

- `bond_angles`: the original program meant to allow 90 and 120 degrees only, but its
  floating-point angle test also lets 60-degree folds through, so the default is `60, 90, 120`.
- `site_weighting`: the original program prefers sites in less crowded regions. It weights each
  candidate by 12 minus the number of neighboring beads from which a chain has already grown,
  but it takes the weights, in its scan order, from the list of all empty neighbors and applies
  them to the shorter list of those with an allowed angle. `legacy` (default) reproduces this;
  `free_neighbors` weights each candidate by its own number of empty neighbors, as intended;
  `uniform` gives every candidate the same probability.

With the defaults, every candidate list and weight matches the original program step by step
(checked on traced runs), and bond-angle fractions, end-to-end distances and bead contacts agree
within statistical error. Three of its faults are not reproduced:

- Image flags are now written, so LAMMPS no longer warns about inconsistent image flags.
- A trapped chain is grown again elsewhere. The original program backed up one bead but kept
  growing next to the removed bead, which can leave a bond of twice the bond length, and stopped
  when this happened twice. In 40 test runs with the paper's settings, 3 contained such a bond
  and 1 stopped.
- A cross-linker can no longer land on an occupied site.

## Performance

Wall time and peak memory of a complete run, including the output files, on one core of an
Apple M4 with Apple Clang 21; the original program was compiled with the same optimization
flags and its screen output discarded.

| system | original program | this code |
|---|---|---|
| fracture paper, No = 20 (105,080 beads) | 0.15 s, 100 MB | 0.024 s, 9 MB |
| ten times larger (1,050,800 beads) | 1.46 s, 470 MB | 0.26 s, 78 MB |

The original program is limited to 300 cells per edge and stops when a chain is trapped a
second time, which happened in every test at occupancies of 0.2 and above or with chains of 1000
beads. This code handles these cases and boxes of up to 1000 cells per edge.

## Code structure

```
include/rw/
  model/   the FCC lattice with its bond angles, and the packed system with its consistency check
  sim/     chain growth and cross-linker placement, and the driver
  io/      configuration and output files
lib/       implementation
apps/      rw_pack
```

Unit tests are in [`tests/randomwalk`](../../tests/randomwalk). The random number generator and
the configuration reader are shared with the other modules (`src/common`).

## References

1. D. Hossain, M. A. Tschopp, D. K. Ward, J. L. Bouvard, P. Wang, M. F. Horstemeyer. Molecular
   dynamics simulations of deformation mechanisms of amorphous polyethylene. *Polymer* **2010**,
   51, 6071–6083. [10.1016/j.polymer.2010.10.009](https://doi.org/10.1016/j.polymer.2010.10.009)
2. M. A. Tschopp, J. L. Bouvard, D. K. Ward, M. F. Horstemeyer. Atomic scale deformation
   mechanisms of amorphous polyethylene under tensile loading. In *Supplemental Proceedings:
   Materials Fabrication, Properties, Characterization, and Modeling* (TMS 2011); Wiley, **2011**;
   pp 789–794. [10.1002/9781118062142.ch95](https://doi.org/10.1002/9781118062142.ch95)
3. Z. Ye, R. A. Riggleman. Molecular View of Cavitation in Model-Solvated Polymer Networks.
   *Macromolecules* **2020**, 53, 7825–7834.
   [10.1021/acs.macromol.0c00994](https://doi.org/10.1021/acs.macromol.0c00994)
4. C. W. Barney, Z. Ye, I. Sacligil, K. R. McLeod, H. Zhang, G. N. Tew, R. A. Riggleman,
   A. J. Crosby. Fracture of model end-linked networks. *Proc. Natl. Acad. Sci. U.S.A.* **2022**,
   119, e2112389119. [10.1073/pnas.2112389119](https://doi.org/10.1073/pnas.2112389119)
5. H. Zhang, R. A. Riggleman. Predicting failure locations in model end-linked polymer networks.
   *Phys. Rev. Mater.* **2024**, 8, 035604.
   [10.1103/PhysRevMaterials.8.035604](https://doi.org/10.1103/PhysRevMaterials.8.035604)
