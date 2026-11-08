# Example 5 — lammps

The FRP and ATRP networks of Zhang et al., *Macromolecules* **2025**, 58, 3168, formed in
coarse-grained molecular dynamics with `fix poly/react` ([`src/lammps`](../../src/lammps/README.md)),
and compared with the 3dbfm and statistical modules.

| File | Purpose |
|:--|:--|
| `in.radical` | LAMMPS input: a Kremer–Grest melt of monomers, cross-linkers and initiators that polymerizes up to a given vinyl conversion |
| `frp.cfg`, `atrp.cfg`, `atrp_paper.cfg` | reaction files written by `polyreact_preset`: FRP with the rate constants of the paper, ATRP with ten times and with the rate constants of the paper |
| `cases.py` | writes the job lists of the comparison: one bash command per run |
| `run_job.sh` | runs one line of a job list (the line number, or `$SLURM_ARRAY_TASK_ID`) |
| `slurm_lattice.sbatch`, `slurm_md.sbatch` | SLURM job arrays for the lattice runs (one core each) and the LAMMPS runs (one node each) |
| `analyze.py` | turns the PSF snapshots of a lattice run into cluster tables, and compares all runs |

## The LAMMPS example

With the plugin built (see [`src/lammps/README.md`](../../src/lammps/README.md)), from this
directory:

```bash
mpirun -np 8 lmp -in in.radical -var plugin /path/to/polyreactplugin.so -var mechanism atrp -var xl 0.08
```

Variables (`-var name value`): `mechanism` (`atrp` or `frp`), `rf` (reaction file, default
`<mechanism>.cfg`), `cells` (FCC cells per box edge, 4 × cells³ beads; default 20, i.e. 32,000),
`xl` (cross-linker fraction, default 0.04), `init` (initiator fraction, default 0.002), `seed`,
`stop` (vinyl conversion at which the run ends, default 0.96), `out` (steps between output lines)
and `tag` (name of the output files `log_<tag>.lammps` and `final_<tag>.data`). The thermo output
gives the vinyl conversion, the radicals, the bond tally and the molecule statistics of the fix.

## The comparison

The networks are formed in three models with the same chemistry and compared at equal vinyl
conversion:

| Model | Module | Space and motion | Partner of a reaction attempt |
|:--|:--|:--|:--|
| 3dbfm | [3dbfm](../../src/3dbfm/README.md) | cubic lattice, bond fluctuation moves (the model of the paper) | one of the 54 nearby lattice positions, at random |
| statistical | [statistical](../../src/statistical/README.md) | none: every unit can meet every other | any unit, with the rates of the elementary reactions |
| LAMMPS | [lammps](../../src/lammps/README.md), `in.radical` | Kremer–Grest melt, Langevin dynamics | one neighbor within 1.2 σ, at random, every 20 steps (`partner random`) |

In 3dbfm and LAMMPS every unit that can react draws a partner, so a radical and a monomer, or two
radicals, have two chances per step or check, one from each side; this keeps the ratios of the rate
constants those of the statistical model.

All three use the units, structural types and reaction categories of the 3dbfm module (the LAMMPS
reaction files reproduce its transition table), a cross-linker fraction [X] of 4, 8, 12, 16 or
20 %, and the rate constants of the paper: FRP `k_p1` = 0.5, `k_p2` = 0.25, `k_term` = 0.95; ATRP
`k_p1` = 0.005, `k_p2` = 0.0025. They are probabilities per Monte Carlo step, per check or rate
coefficients, so times differ between the models and only quantities at equal conversion are
compared. In LAMMPS, ATRP runs by default with ten times these rate constants to keep runs
affordable; `--rates paper` runs it with the paper's values, to check that the structure does not
depend on this choice.

The comparison runs in two stages:

| Stage | Size | Composition | Where |
|:--|:--|:--|:--|
| building up (`--size mini`) | about 4,000 units: 3dbfm 34³ lattice with 4,176 units, statistical 4,176 units, LAMMPS 4,000 beads | 2 % initiator, FRP `k_init` = 10⁻⁴; 5 samples, up to conversion 0.8 | a workstation, in minutes |
| validation (`--size paper`) | the paper: 3dbfm 100³ lattice with 106,250 units, statistical 106,250 units, LAMMPS 108,000 beads | 0.2 % initiator, FRP `k_init` = 10⁻⁵; 5 samples, up to conversion 0.96, as in the paper | a cluster |

What a small box must hold is many chains, not many units: at the 0.2 % initiator of the paper,
4,000 units hold eight chains of some 500 units, too few for a gel point. The mini systems
therefore have ten times more initiator, and FRP a ten times faster initiation, so that
`k_init`/`k_p1` = 2 × 10⁻⁴, the largest ratio that keeps FRP kinetics realistic according to
Furuya and Koga; all three models get the same composition, so the comparison stays one to one.
As in our study of AB₂ polymerization (Zhang et al., *Adv. Theory Simul.* **2019**, 2, 1800102),
many small samples replace one large system.

### Runs on a workstation

From an empty run directory, with the toolkit built in `build/` and the plugin in `build-plugin/`
(or wherever `POLY_BUILD` points):

```bash
E=/path/to/Polymerization/examples/example_5
python3 $E/cases.py lattice --size mini --x 4 12 20 > jobs_lattice.txt
python3 $E/cases.py md --size mini --x 4 12 20 > jobs_md.txt
seq $(wc -l < jobs_lattice.txt) | xargs -P 4 -I{} $E/run_job.sh jobs_lattice.txt {}
seq $(wc -l < jobs_md.txt) | LMP=/path/to/lmp MPIRUN= xargs -P 4 -I{} $E/run_job.sh jobs_md.txt {}
python3 $E/analyze.py compare . > comparison.md
```

`MPIRUN=` (empty) runs each LAMMPS sample as one process, and `xargs -P 4` runs four side by side,
one per performance core of an Apple M4 Mac mini, which leaves the machine usable. Mini runs stop at
conversion 0.8 for FRP and 0.96 for ATRP, whose gel point comes close to 0.8 at 4 % cross-linker; a
LAMMPS run takes about 30 s, a lattice run seconds to half a minute, and the 30 LAMMPS and 60
lattice runs about seven minutes in all. More samples (`--samples 15`) narrow the error bars.

### Runs on a cluster

For the validation, build the toolkit and the plugin on the cluster (see the
[repository README](../../README.md) and [`src/lammps/README.md`](../../src/lammps/README.md)), then
in an empty run directory:

```bash
export REPO=/path/to/Polymerization            # POLY_BUILD defaults to $REPO/build
python3 $REPO/examples/example_5/cases.py lattice > jobs_lattice.txt   # 3dbfm and statistical, paper size
python3 $REPO/examples/example_5/cases.py md > jobs_md.txt             # LAMMPS
python3 $REPO/examples/example_5/cases.py md --rates paper --mechanisms atrp --x 8 20 --samples 2 > jobs_md_paper.txt
mkdir -p logs
sbatch --array=1-$(wc -l < jobs_lattice.txt) $REPO/examples/example_5/slurm_lattice.sbatch
sbatch --array=1-$(wc -l < jobs_md.txt) --export=ALL,LMP=/path/to/lmp $REPO/examples/example_5/slurm_md.sbatch
```

The templates need the partition and account of the cluster. Where whole nodes are allocated, the
lattice list can run inside one node job instead, e.g.
`seq $(wc -l < jobs_lattice.txt) | xargs -P $SLURM_NTASKS -I{} bash $REPO/examples/example_5/run_job.sh jobs_lattice.txt {}`.
The LAMMPS runs need `LMP`, `PLUGIN` (default `$POLY_BUILD/src/lammps/polyreactplugin.so`) and
`MPIRUN` (`srun`, or `ibrun` on TACC systems).

`cases.py` options: `--size paper` (default), `reduced` (68³ with 33,408 units, 33,408, 32,000
beads; composition of the paper), `mini` or `smoke` (tiny, for a test of the pipeline; composition
of the mini systems); `--init` and `--k-init` (initiator fraction and FRP initiation, defaults as
above); `--samples` (default 5, as in the paper), `--x`, `--mechanisms`, `--stop` (default 0.96, as
in the paper) and, for the lattice list, `--models bfm stat`. The LAMMPS runs write their reaction
file (`react.cfg`) with `polyreact_preset` from these rate constants.

Every run writes into `bfm/`, `stat/`, `md/` or `md_paper/` under `<mechanism>_X<x>_s<sample>`.
The lattice runs turn their PSF snapshots (every 0.01 in conversion) into `clusters_<s>.csv` and
delete them; the LAMMPS runs write the same quantities to their thermo output (`log_*.lammps`) and
the final network to `final_*.data`.

For planning, on one core (lattice) or six cores (LAMMPS) of an Apple M4 at the reduced size: 3dbfm
ATRP reaches a conversion of 0.74 after 675,000 Monte Carlo steps in about ten minutes, and FRP 0.91
in the same time; the statistical model takes seconds; LAMMPS ATRP with 8 % cross-linker reaches 0.8
after 440,000 steps in about twenty minutes. The paper size is about three times larger, and the
last part of the conversion up to 0.96 is the slowest.

### Analysis

```bash
python3 $REPO/examples/example_5/analyze.py compare . > comparison.md
```

For every mechanism and cross-linker fraction, as mean ± standard error over the samples:

- the gel point, the vinyl conversion at the maximum of the reduced degree of polymerization (Mw
  of the molecules of two or more units without the largest), as in the paper;
- the number of molecules at conversion 0.05 (for ATRP, the number of living chains that have
  not joined yet);
- the fraction of new bonds within a molecule (rings and intramolecular cross-links) at conversion
  0.2 and 0.8, and of bonds between two polymers at 0.8; a bond counts as monomer–monomer,
  polymer–monomer, polymer–polymer or intramolecular as in the lattice modules (a unit without
  bonds is a monomer).

A gel point is only reported when RDP has passed its maximum before the end of the run; runs that
end earlier are counted separately.

## Results: mini systems

From the commands under *Runs on a workstation* (5 samples per case; 3dbfm and statistical with
4,176 units, LAMMPS with 4,000 beads; 2 % initiator, so 83 and 80 chains). Mean ± standard error.

Gel point, vinyl conversion at the maximum of RDP:

| [X] | ATRP 3dbfm | ATRP statistical | ATRP LAMMPS | FRP 3dbfm | FRP statistical | FRP LAMMPS |
|:--|--:|--:|--:|--:|--:|--:|
| 4 % | 0.880 ± 0.024 | 0.820 ± 0.016 | 0.884 ± 0.011 | 0.39 ± 0.06 | 0.39 ± 0.04 | 0.40 ± 0.08 |
| 12 % | 0.614 ± 0.016 | 0.498 ± 0.018 | 0.607 ± 0.024 | 0.20 ± 0.03 | 0.26 ± 0.06 | 0.22 ± 0.06 |
| 20 % | 0.488 ± 0.015 | 0.408 ± 0.009 | 0.473 ± 0.028 | 0.14 ± 0.02 | 0.12 ± 0.03 | 0.21 ± 0.03 |

New bonds within a molecule (rings and intramolecular cross-links), as a fraction of all new bonds
up to conversion 0.2 and 0.8:

| [X] | | 3dbfm | statistical | LAMMPS |
|:--|:--|--:|--:|--:|
| 4 % | ATRP, 0.2 | 0.0013 ± 0.0000 | 0.0000 ± 0.0000 | 0.0003 ± 0.0003 |
| | ATRP, 0.8 | 0.0048 ± 0.0005 | 0.0005 ± 0.0002 | 0.0042 ± 0.0003 |
| 12 % | ATRP, 0.2 | 0.0030 ± 0.0009 | 0.0000 ± 0.0000 | 0.0024 ± 0.0009 |
| | ATRP, 0.8 | 0.0314 ± 0.0014 | 0.0293 ± 0.0014 | 0.0332 ± 0.0012 |
| 20 % | ATRP, 0.2 | 0.0059 ± 0.0010 | 0.0002 ± 0.0002 | 0.0062 ± 0.0010 |
| | ATRP, 0.8 | 0.0689 ± 0.0008 | 0.0665 ± 0.0005 | 0.0721 ± 0.0010 |
| 4 % | FRP, 0.2 | 0.0036 ± 0.0009 | 0.0000 ± 0.0000 | 0.0040 ± 0.0008 |
| | FRP, 0.8 | 0.0157 ± 0.0004 | 0.0121 ± 0.0008 | 0.0165 ± 0.0006 |
| 12 % | FRP, 0.2 | 0.0091 ± 0.0017 | 0.0047 ± 0.0006 | 0.0134 ± 0.0024 |
| | FRP, 0.8 | 0.0545 ± 0.0020 | 0.0477 ± 0.0009 | 0.0562 ± 0.0016 |
| 20 % | FRP, 0.2 | 0.0178 ± 0.0035 | 0.0091 ± 0.0013 | 0.0264 ± 0.0027 |
| | FRP, 0.8 | 0.0917 ± 0.0011 | 0.0843 ± 0.0006 | 0.0972 ± 0.0011 |

Molecules of two or more units at conversion 0.05: in ATRP 67, 71 and 72 (3dbfm), 66, 70 and 71
(statistical) and 68, 70 and 72 (LAMMPS) at 4, 12 and 20 % cross-linker; in FRP two or three in
every model. Bonds between two polymers make up 1.2–2.2 % of all new bonds at conversion 0.8 in
ATRP and 0.4–0.6 % in FRP, alike in the three models.

- **ATRP.** LAMMPS gels at the same conversion as 3dbfm at every cross-linker fraction. The
  statistical model, without space, gels earlier by 0.06–0.12: as in the dynamic lattice liquid
  study of Lyu et al. (*Macromolecules* **2018**, 51, 6673), spatial restrictions delay gelation.
- **FRP.** The gel points agree within their large scatter, which comes from the few chains that
  grow at a time in a small box (two or three at conversion 0.05); at 20 % cross-linker LAMMPS gels
  somewhat later.
- **Rings.** Both spatial models form rings from the start, the statistical model hardly before
  the network is large. LAMMPS forms somewhat more rings than 3dbfm early in FRP (0.026 against
  0.018 at 20 % cross-linker and conversion 0.2), a difference of the motion of the beads.
- The models need not agree: they move the units differently, and the ratio of reaction to
  diffusion shapes the networks (Zhang et al., *Macromolecules* **2025**, 58, 3168). What must agree
  is the chemistry, which the reaction files reproduce exactly, and the meaning of the rate
  constants: in all three, a radical and a monomer, or two radicals, have two chances per step or
  check, one from each side. This comparison found that `partner random` first let only the
  radicals draw a partner, which doubled termination relative to propagation; FRP then gelled at
  0.49, 0.29 and 0.22. With every reactive bead drawing a partner, as now, the FRP gel points agree.

These mini systems use 2 % initiator and `k_init` = 10⁻⁴ instead of the paper's 0.2 % and 10⁻⁵, and
ten times the paper's ATRP rate constants in LAMMPS; the comparison at the size and composition of
the paper is to be run on a cluster.
