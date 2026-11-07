# Example 3 — dll

| file | system |
|---|---|
| `atrp_gel.ini` | ATRP copolymerization of a monomer and a cross-linker in bulk, [I]₀ : [X]₀ : [M]₀ = 1 : 2 : 100 (131,072 sites; about 25 s) |
| `divinyl_2018.ini` | ATRP homopolymerization of a divinyl monomer in bulk, reproducing the DLL gel points of Lyu et al., *Macromolecules* 2018 (about 10 s) |
| `frp_solution.ini` | FRP copolymerization of a monomer and a cross-linker at 50 % solvent (55,296 sites; about 35 s) |
| `solution.ini` | chains of 16 beads in solvent at polymer fraction 0.5, no reactions |
| `smoke_react.ini` | small reactive system for the smoke tests (ATRP; FRP with `mechanism=frp`) |
| `smoke_small.ini` | small system of chains in solvent for the smoke tests |

From the repository root:

```bash
build/src/dll/dll_react examples/example_3/atrp_gel.ini
build/src/dll/dll_react examples/example_3/atrp_gel.ini reaction_partner=all_neighbors k_p1=0.005 k_p2=0.005
build/src/dll/dll_react examples/example_3/frp_solution.ini samples=5
build/src/dll/dll_react examples/example_3/solution.ini initial=aligned n_chains=1024   # a melt, no reactions
```

The gel point can be found from the PSF snapshots, for example as the maximum of the reduced
degree of polymerization. All keys are described in [`src/dll/README.md`](../../src/dll/README.md).
