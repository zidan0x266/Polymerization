# common — code shared by the modules

Everything here is used by at least two modules and knows nothing about any one of them. Code
moves here when a second module needs it; until then it stays in its module. The LAMMPS fix in
[`src/lammps`](../lammps) is licensed separately and uses nothing from this directory.

| Folder | Namespace | Contents | Used by |
|:--|:--|:--|:--|
| `poly/` | `poly`, `poly::ini` | random numbers (xoshiro256++), configuration-file reader, index set | all modules |
| `poly/geometry/` | `poly::geometry` | integer lattice vector, periodic face-centered cubic (FCC) lattice | 3dbfm, dll, randomwalk |
| `poly/topology/` | `poly::topology` | molecules as disjoint sets of units (union-find), bond lists of fixed capacity | 3dbfm, dll, statistical |
| `poly/chem/` | `poly::chem` | unit states, reaction schemes (FRP, ATRP), reaction bookkeeping, per-sample output tables | 3dbfm, dll, statistical |
| `poly/io/` | `poly::io` | output files that report failure, PDB and PSF writers | all modules |

The chemistry is the static library `poly::chemistry` (`lib/`); the rest is header-only and comes
with `poly::common`.

The three directions the toolkit grows in map onto this layout: a new simulation method is a new
module under `src/`, a new chemistry is a new reaction scheme in `poly/chem/`, and a new polymer
system (linear, branched, network, copolymer) is a new example, not new code.
