# Example 3 — dll

| file | system |
|---|---|
| `solution.ini` | chains of 16 beads in solvent at polymer fraction 0.5 |
| `smoke_small.ini` | small system that runs in a fraction of a second |

From the repository root:

```bash
build/src/dll/dll_react examples/example_3/solution.ini
build/src/dll/dll_react examples/example_3/solution.ini initial=aligned n_chains=1024   # a melt
```

The module is in phase 1: chains and solvent move, reactions are in preparation. All keys are
described in [`src/dll/README.md`](../../src/dll/README.md).
