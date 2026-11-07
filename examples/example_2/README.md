# Example 2 — randomwalk

| file | purpose |
|---|---|
| `endlink_paper.ini` | precursors of the end-linked network with No = 20 (Zhang et al., *Macromolecules* 2026, Table 1) |
| `smoke_small.ini` | small system that packs in milliseconds |
| `in.prepare` | LAMMPS: relaxation, equilibration and compression to number density 0.85 |
| `in.react` | LAMMPS: end-linking with `fix bond/create` |

From the repository root:

```bash
build/src/randomwalk/rw_pack examples/example_2/endlink_paper.ini
cd endlink_No20
lmp -in ../examples/example_2/in.prepare -var data config.lmp
lmp -in ../examples/example_2/in.react
```

The other systems of Table 1 differ in three keys, given on the command line:

| No | `chain_length` | `n_chains` | `n_crosslinkers` |
|--:|--:|--:|--:|
| 20 | 18 | 5680 | 2840 |
| 10 | 8 | 12416 | 6208 |
| 7 | 6 | 19520 | 9760 |
| 6 | 5 | 27154 | 13577 |
| 5 | 4 | 35020 | 17510 |

```bash
build/src/randomwalk/rw_pack examples/example_2/endlink_paper.ini chain_length=8 n_chains=12416 n_crosslinkers=6208 directory=endlink_No10
```

The LAMMPS scripts need LAMMPS 30 Sep 2026 or newer, built with the MOLECULE and MC packages.
They use the step counts of the paper (610,000 steps of preparation and up to 120 million steps
of end-linking, after which more than 96 % of the cross-linker sites have reacted); shorter
runs are set with `-var pre_steps …`, `-var compress_steps …` or `-var react_steps …`, and the
seed with `-var seed …`. All keys of `rw_pack` are described in
[`src/randomwalk/README.md`](../../src/randomwalk/README.md).
