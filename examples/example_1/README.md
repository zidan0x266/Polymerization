# Example 1 — 3dbfm

| file | system |
|---|---|
| `atrp_paper.ini` | ATRP network, 4 % cross-linker (Zhang et al., *Macromolecules* 2025) |
| `frp_paper.ini` | FRP network, 4 % cross-linker (same paper) |
| `smoke_small.ini` | small 32³ system that finishes in seconds |

From the repository root:

```bash
build/src/3dbfm/bfm_react examples/example_1/atrp_paper.ini
```

Other cross-linker contents are set with `monovinyl_fraction` (e.g. 0.92 for 8 %). All keys are
described in [`src/3dbfm/README.md`](../../src/3dbfm/README.md).
