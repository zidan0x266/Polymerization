# lammps — fix poly/react

`fix poly/react` brings the chemistry of this toolkit into coarse-grained molecular dynamics with
[LAMMPS](https://www.lammps.org): beads react with neighbors according to a reaction file, bonds
form and break, and bead types change, while LAMMPS moves the beads. One fix serves different
mechanisms (radical chain growth, step growth, end-linking, reversible reactions), written as
sets of reactions in the format of ChemLab (refs. 1, 2).

**Status — on-going.** Phase 1 is done: the reaction files, the rule engine that decides the
reactions, and the reaction files of free radical (FRP) and atom transfer radical polymerization
(ATRP) with the chemistry of the [3dbfm module](../3dbfm/README.md), tested without LAMMPS. Next:

| Phase | Content |
|:--|:--|
| 2 | the fix as a LAMMPS plugin: neighbor search, one reaction per bead per check across processors, bonds, angles and type changes, molecule tracking, output |
| 3 | FRP and ATRP against the 3dbfm and statistical modules (Macromolecules 2025 cases) |
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
| `A(…):B(…) -> A'(…) + B'(…)` | dissociation: the bond between two bonded beads breaks |
| `A(…):B(…) + C(…) -> A'(…):C'(…) + B'(…)` | exchange (parsed; carried out in a later phase) |
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

## Checks

The unit tests ([`tests`](tests)) check the parsing of equations and files, the matching of
types, functionalities, distances and bond conditions, the three competition rules against their
exact probabilities, and the radical reaction files against the transition table of the original
3DBFM code: for every pair of structural units, the same products with the same probabilities, for
FRP (including the legacy competition between propagation and termination) and ATRP, and the same
reaction categories. All twelve reaction files of the ChemLab examples are read without change.

## Code structure

```
include/polyreact/
  reaction_file.hpp   reaction files and equations
  engine.hpp          the decision on reactions for a pair or a single bead
  presets.hpp         reaction files of radical polymerization
lib/                  implementation
apps/                 polyreact_preset
tests/                unit tests
```

The fix will follow in `fix_poly_react.{h,cpp}`, built as a LAMMPS plugin.

## Credits and references

The reaction format and the approach — reactions as rules on bead types and functionalities,
checked every few MD steps within a capture radius, one reaction per bead, partners drawn at random
— follow ChemLab by Jakub Krajniak and the papers below; the code here is new.

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
