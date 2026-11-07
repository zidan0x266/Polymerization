# Example 4 — statistical

| file | system |
|---|---|
| `divinyl_2018.ini` | ATRP homopolymerization of a divinyl monomer, the "w.c." model of Lyu et al., *Macromolecules* 2018 (200,000 units; a few seconds, about 0.4 GB of PSF snapshots) |
| `frp_copolymer.ini` | FRP copolymerization of a monomer and 4 % cross-linker, the composition of the 3dbfm FRP example (106,250 units; a few seconds) |
| `smoke_small.ini` | small system for the smoke tests (ATRP; FRP with `mechanism=frp`) |

From the repository root:

```bash
build/src/statistical/stat_react examples/example_4/divinyl_2018.ini
build/src/statistical/stat_react examples/example_4/divinyl_2018.ini cyclization=false   # "wo.c." model
build/src/statistical/stat_react examples/example_4/frp_copolymer.ini samples=5
```

The gel point can be found from the PSF snapshots, for example as the maximum of the reduced
degree of polymerization. All keys are described in
[`src/statistical/README.md`](../../src/statistical/README.md).
