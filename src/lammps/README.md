# lammps — fix poly/react

`fix poly/react` brings the chemistry of this toolkit into coarse-grained molecular dynamics with
[LAMMPS](https://www.lammps.org): beads react with neighbors according to a reaction file, bonds
form and break, and bead types change, while LAMMPS moves the beads. One fix serves different
mechanisms (radical chain growth, step growth, end-linking, reversible reactions), written as
sets of reactions in the format of ChemLab (refs. 1, 2).

**Status — on-going.** Phases 1 and 2 are done: the reaction files and the rule engine that
decides the reactions, the reaction files of free radical (FRP) and atom transfer radical
polymerization (ATRP) with the chemistry of the [3dbfm module](../3dbfm/README.md), and the fix
itself as a LAMMPS plugin, for association and unimolecular reactions on any number of processors.
Phase 3 is under way: the fix reports the bond tally and the molecule statistics of the lattice
modules, and [example 5](../../examples/example_5/README.md) forms the FRP and ATRP networks of
ref. 3 in a Kremer–Grest melt and compares them with the 3dbfm and statistical modules; the
comparison runs on a cluster. Next:

| Phase | Content |
|:--|:--|
| 3 | FRP and ATRP against the 3dbfm and statistical modules (Macromolecules 2025 cases): results of the comparison |
| 4 | step growth: AB₂ hyperbranched polymers (ref. 2), A₂ + B₃, end-linking |
| 5 | reversible reactions and byproducts, explicit ATRP activation, performance |

## License

This directory is licensed under the **GNU General Public License, version 2 only**
([LICENSE](LICENSE)), the license of LAMMPS, so that the fix can be built into LAMMPS and
distributed with it. The rest of the repository keeps its PolyForm Noncommercial license; this
directory uses none of its code.

## Reaction files

An INI file with four kinds of sections; keys are written `key = value` or `key: value`, and lines
starting with `#` or `;` are comments.

```ini
[general]
interval = 200            # MD steps between two checks for reactions
competition = independent # how one of several possible reactions for a pair is chosen

[group_growth]
bond_type = 1             # LAMMPS bond type of the bonds created by the reactions of this group
angle_type = 1            # optional: LAMMPS angle type of the angles a new bond closes

[reaction_propagation]
reaction = DA(3, 4) + MA(1, 2) -> RA(1):DA(2)
rate = 0.8                # probability per check
cutoff = 1.2              # capture radius
intramolecular = true     # allow bonds within a molecule (rings)
group = growth
```

A bead is a type (a LAMMPS type label) with an integer functionality. On the left side of an
equation, `T(min, max)` accepts a bead of type `T` with functionality in [min, max); on the right
side, `T'(d)` gives the new type and the change of functionality. A type without parentheses
accepts any functionality and keeps it.

| Equation | Reaction |
|:--|:--|
| `A(…) + B(…) -> A'(…):B'(…)` | association: a bond forms between two beads within the cutoff |
| `A(…):B(…) -> A'(…) + B'(…)` | dissociation: the bond between two bonded beads breaks (parsed; carried out in phase 5) |
| `A(…):B(…) + C(…) -> A'(…):C'(…) + B'(…)` | exchange (parsed; carried out in phase 5) |
| `A(…) -> A'(…)` | unimolecular, e.g. initiator decomposition (not in ChemLab) |

Keys of a reaction: `reaction`, `rate`, `cutoff`, `min_cutoff`, `group`, `intramolecular`
(bonds within a molecule; default `false`), `intraresidual` (bonds within a residue, the initial
molecule; default `false`), `virtual` (change the beads without a bond), `active`, `diss_rate`
(dissociation probability, if not given by `rate`), `exclude_extensions`. The defaults are those of
ChemLab, so its reaction files can be used as they are; keys that the fix does not use yet
(`potential`, `potential_options`, `sigma`, `alpha`, …) are accepted and reported. The bond
potential itself is set in the LAMMPS input script.

When more than one reaction is possible for a pair, `competition` decides:

| Value | Rule |
|:--|:--|
| `independent` | each is tried with its rate; one of those that pass is taken at random (ChemLab) |
| `partition` | one draw divides [0, 1) by the rates, rescaled if they sum to more than one |
| `legacy` | for two reactions a and b (in the order of the file), the rule of the original 3DBFM code: a with p = r_a, or r_a/(r_a + r_b) if r_a + r_b ≥ 1; otherwise b with probability (1 − p) min(1, r_b + p) |

## Build the plugin

The fix is a [LAMMPS plugin](https://docs.lammps.org/plugin.html): LAMMPS must be built with the
PLUGIN and MOLECULE packages (`-D PKG_PLUGIN=on -D PKG_MOLECULE=on`). From the repository root,
with the same compiler as LAMMPS:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPOLYREACT_LAMMPS_SOURCE=/path/to/lammps/src \
      -DPOLYREACT_LAMMPS_EXECUTABLE=/path/to/lammps/build/lmp   # optional, for the tests
cmake --build build --target polyreactplugin                    # build/src/lammps/polyreactplugin.so
```

| Option | Meaning |
|:--|:--|
| `POLYREACT_LAMMPS_SOURCE` | the `src` directory of the LAMMPS the plugin is loaded into |
| `POLYREACT_LAMMPS_CXX_STANDARD` | the C++ standard of that LAMMPS build (default `17`, the LAMMPS default); it must match, since LAMMPS formats its messages with its own copy of fmt in C++17 and with `std::format` from C++20 on |
| `POLYREACT_LAMMPS_SIZES` | `LAMMPS_SIZES` of that build (default `smallbig`) |
| `POLYREACT_LAMMPS_MPI` | `ON` (default) for LAMMPS with MPI, `OFF` for a serial LAMMPS with the MPI stubs |
| `POLYREACT_LAMMPS_EXECUTABLE` | that `lmp`, to run the plugin tests with `ctest` |

Developed and tested with LAMMPS 30 Sep 2026 (Update 1).

## Use in LAMMPS

```
plugin load /path/to/polyreactplugin.so
...
labelmap atom 1 MON 2 PRI 3 ONE 4 TWO 5 THR 6 FOU 7 FIV 8 SIX 9 SEV 10 EIG 11 INI
fix react all poly/react frp.cfg 2025 every 20 log reactions.txt
thermo_style custom step f_react[1] f_react[3] f_react[7] f_react[10] f_react[12] bonds angles
```

```
fix ID group-ID poly/react file seed keyword value ...
```

`file` is the reaction file and `seed` a positive integer for the random numbers. Keywords:

| Keyword | Meaning |
|:--|:--|
| `every N` | MD steps between two checks (default: `interval` of the file) |
| `partner all` | every pair of beads within the cutoff is tried once per check (default) |
| `partner random` | every bead that can take part in a pair reaction tries one neighbor in the group within the cutoff, chosen at random: as in the lattice models, a pair of such beads (a radical and a monomer, or two radicals) has two chances per check, one from each bead |
| `functionality name` | the integer per-atom vector `i_name` (`fix property/atom`) that holds the functionalities; needed when the reactions use them, otherwise all are 0 |
| `max N` | at most N reactions per check (default: `max_per_interval` of the file, unlimited) |
| `log file` | write every reaction: step, reaction, center and partner atom IDs (partner 0 for a unimolecular reaction), distance, and 1 if the bond closes a ring |

The bead types of the reaction file are atom type labels (`labelmap atom`) or atom type numbers.
Only beads in the group react, and only with beads in the group. The fix needs a molecular atom
style with bonds (`bond`, `angle`, `molecular`, `full`); room for the new bonds, angles and special
neighbors (`extra/bond/per/atom`, `extra/angle/per/atom`, `extra/special/per/atom` of `read_data` or
`create_box`); and ghost atoms out to the largest reaction cutoff plus the neighbor skin
(`comm_modify cutoff`). A residue is an initial molecule, given by the molecule ID; ID 0 means no
residue.

At every check, right after the positions are updated:

1. each processor tries its own beads, the unimolecular reactions and the pairs with their
   neighbors, with the rules of the file; the reactions that pass are proposals;
2. the proposals of all processors are gathered on every processor;
3. every processor goes through them in the same order, random or, with `nearest`, by distance,
   and accepts a proposal if neither bead has reacted in this check yet: each bead reacts at most
   once per check, and the outcome does not depend on the number of processors except through the
   random numbers. A bond within a molecule is refused if the reaction is not `intramolecular`,
   also when bonds accepted earlier in the same check have joined the two molecules;
4. each processor carries out the accepted reactions on its beads: types, functionalities, the new
   bond, the special lists as LAMMPS builds them for the `special_bonds` weights (only 1-2 neighbors
   when the 1-3 and 1-4 weights are 1, as with FENE), and the angles the bond closes, with the angle
   type of its group (if an angle style is defined), for `newton_bond` on and off. LAMMPS then
   reneighbors at once.

Molecules are followed with a union-find over atom IDs that every processor keeps, together with
the number of bonds of every atom. The fix provides a global vector for `thermo` and other
commands:

| Element | Meaning |
|:--|:--|
| 1 | reactions in all |
| 2 | reactions in the last check |
| 3 | bonds created, of which: |
| 4 | between two beads without bonds (monomer–monomer) |
| 5 | between a bead with and a bead without bonds (polymer–monomer) |
| 6 | between beads with bonds in different molecules (polymer–polymer) |
| 7 | within a molecule (rings) |
| 8 | molecules of two or more atoms |
| 9, 10 | their number- and weight-average size Mn and Mw, in atoms |
| 11 | the largest molecule |
| 12 | Mw without the largest molecule, the reduced degree of polymerization, whose maximum marks the gel point |
| 13, 14, … | reactions of each `[reaction_…]` section, in the order of the file |

Elements 4–7 classify new bonds as the lattice modules do. Elements 8–12 count all atoms,
are computed only when asked for, and take a pass over all atoms each time.

The counts (elements 1–7 and 13 on) are kept in restart files; after `read_restart`, define the
fix again with the same ID.
As with other LAMMPS fixes, the state of the random numbers is not kept.

Not yet supported: dissociation and exchange reactions (phase 5) and `run_style respa`. Bonds
changed by other fixes (`fix bond/create`, `bond/break`, …) are not seen by the molecule tracking.
Since every processor holds all proposals of a check and the molecules of all atoms, very large
systems with many reactions per check will need the distributed version planned for phase 5. In a
Kremer–Grest melt of 32,000 beads on 4 processes with a check every 20 steps, the run takes 42 %
longer than without the fix; about half of that is the reneighboring after checks with reactions,
as with `fix bond/create`.

## Radical polymerization

`polyreact_preset` writes the reaction file for radical copolymerization of monovinyl monomers
and divinyl cross-linkers with the chemistry of the 3dbfm module:

```bash
build/src/lammps/polyreact_preset frp k_init=1e-5 k_p1=0.5 k_p2=0.25 k_term=0.95 > frp.cfg
build/src/lammps/polyreact_preset atrp k_p1=0.005 k_p2=0.0025 cutoff=1.2 interval=100 > atrp.cfg
```

Keys: the rate constants `k_init`, `k_p1`, `k_px` (unreacted cross-linker; default `k_p1`),
`k_p2`, `k_term`, as probabilities per check; `interval`, `cutoff`, `bond_type`, `angle_type`,
`competition` (default `legacy`) and `cyclization` (default `true`). Each bead is a monomer unit,
and its type is its structural unit in the original code:

| Type | Unit | Type | Unit |
|:--|:--|:--|:--|
| `MON` | unreacted monovinyl monomer | `FOU` | cross-linker with a radical and a pendant vinyl |
| `PRI` | initiator fragment with its radical | `FIV` | cross-linker with a pendant vinyl |
| `ONE` | monovinyl unit with a radical | `SIX` | cross-linker with a radical, both vinyls reacted |
| `TWO` | saturated unit or used initiator | `SEV` | cross-linker with two radicals |
| `THR` | unreacted divinyl cross-linker | `EIG` | fully reacted cross-linker |
| `INI` | dormant initiator (FRP) | | |

The reactions are listed in the order of the reaction categories of the 3dbfm module (propagation
1–20, termination 21–35, then initiation), so the counts per reaction compare column by column.

[Example 5](../../examples/example_5/README.md) runs these reactions in a Kremer–Grest melt (FENE bonds,
WCA repulsion, cosine bending with k = 1.5, ρ = 0.85, Langevin thermostat, Δt = 0.005, as in ref. 3),
with the composition of the 3dbfm examples: 0.2 % initiator and a cross-linker fraction `xl`. The
beads start on an FCC lattice and are melted; every 20 steps, each bead that can react (radicals,
monomers, cross-linkers with a vinyl group, initiators) tries one random partner within 1.2 σ. The vinyl conversion is computed from the bead types as in the lattice
modules, the run stops at the conversion `stop`, and the thermo output gives the conversion, the
radicals, the bond tally and the molecule statistics:

```bash
cd examples/example_5
mpirun -np 8 lmp -in in.radical -var plugin /path/to/polyreactplugin.so -var mechanism atrp -var xl 0.08
```

`atrp.cfg` and `frp.cfg` were written by `polyreact_preset`, with the rate constants of ref. 3 for
FRP (`k_init` = 10⁻⁵, `k_p1` = 0.5, `k_p2` = 0.25, `k_term` = 0.95) and ten times those of ref. 3
for ATRP (`k_p1` = 0.05, `k_p2` = 0.025), per check, which keeps the runs affordable;
`atrp_paper.cfg` has the ATRP rate constants of ref. 3 (`-var rf atrp_paper.cfg`).

## Checks

The unit tests ([`tests`](tests)) check the parsing of equations and files, the matching of
types, functionalities, distances and bond conditions, the three competition rules against their
exact probabilities, and the radical reaction files against the transition table of the original
3DBFM code: for every pair of structural units, the same products with the same probabilities, for
FRP (including the legacy competition between propagation and termination) and ATRP, and the same
reaction categories. All twelve reaction files of the ChemLab examples are read without change.

The plugin tests ([`tests/plugin`](tests/plugin), `ctest -R lammps.plugin`) run LAMMPS on 1 and 4
processors (also checked on 8) with systems whose outcome is known:

| Test | Check |
|:--|:--|
| pairs | 1000 isolated A–B pairs, also across periodic and processor boundaries: bonds only between partners, with the new types |
| probabilities | the probability per check of the pair reaction and of the competing decay B → D of the same bead, from 5 runs, both partner modes (with `partner random` both beads try the pair), within 4σ of the exact value |
| triplets | B–A–B: A takes one bond per check up to its functionality window and closes the angle B–A–B; the first bonds are counted as monomer–monomer, the second as polymer–monomer, and the molecule statistics are exact |
| squares | two bonded dimers side by side: the second bond closes a ring and is refused unless the reaction is intramolecular, also when both bonds are proposed in the same check; the first bond is counted as polymer–polymer and, with `intramolecular = true`, the second as a ring |
| kg | a Kremer–Grest melt polymerizing with the FRP preset for 3000 steps: the numbers of bonds and angles, and the pair, bond and angle energies, equal (to 1e-10) those of the final configuration read again by LAMMPS, which builds the special lists from scratch; with FENE special bonds, with 1-3 and 1-4 weights, and with `newton_bond` on and off. The bond tally adds up to the bonds created, and the molecule statistics equal those of the bonds in the data file |

## Code structure

```
include/polyreact/
  reaction_file.hpp   reaction files and equations
  engine.hpp          the decision on reactions for a pair or a single bead
  presets.hpp         reaction files of radical polymerization
lib/                  implementation
apps/                 polyreact_preset
plugin/
  fix_poly_react.*    the fix
  polyreact_plugin.cpp  registration with LAMMPS
tests/                unit tests
tests/plugin/         LAMMPS inputs, reaction files and the test driver run_tests.py
```

The plugin code compiles as C++17 (the engine headers it includes are C++17) and links the rule
engine, which is C++23.

## Credits and references

The structure of `fix poly/react` borrows the ideas of ChemLab, the reactive coarse-grained
molecular dynamics code of Jakub Krajniak, Zidan Zhang, Erik Nies and Giovanni Samaey, KU Leuven
(see the [acknowledgements](../../README.md#acknowledgements) of the repository): the reaction format, and reactions as rules on bead types and functionalities,
checked every few MD steps within a capture radius, one reaction per bead, partners drawn at
random. The code here is new. The updates of bond storage, special lists and angles in the fix
follow `fix bond/create` of LAMMPS (MC package, GPL-2.0); `rebuild_special_one`, `create_angles`
and `dedup` are adapted from it.

1. J. Krajniak, Z. Zhang, S. Pandiyan, E. Nies, G. Samaey. Coarse-grained molecular dynamics
   simulations of polymerization with forward and backward reactions. *J. Comput. Chem.* **2018**,
   39, 1764–1778. [10.1002/jcc.25348](https://doi.org/10.1002/jcc.25348)
2. Z. Zhang, J. Krajniak, G. Samaey, E. Nies. A parallel multiscale simulation framework for
   complex polymerization: AB₂-type monomer hyperbranched polymerization as an example. *Adv.
   Theory Simul.* **2019**, 2, 1800102. [10.1002/adts.201800102](https://doi.org/10.1002/adts.201800102)
3. Z. Zhang, J. Krajniak, A. Z. Dookhith, Y. Tian, H. S. Sachar, N. Marioni, T. J. Duncan, J. Liu,
   G. E. Sanoja, V. Ganesan. Topology and Mechanical Properties of Polymer Networks Formed under
   Free Radical and Atom Transfer Radical Polymerizations. *Macromolecules* **2025**, 58,
   3168–3187. [10.1021/acs.macromol.4c03138](https://doi.org/10.1021/acs.macromol.4c03138)
