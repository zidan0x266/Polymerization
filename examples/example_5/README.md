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
one per performance core of an Apple M4 Mac mini, which leaves the machine usable. A LAMMPS run
takes about 30 s (ATRP) or 50 s (FRP) up to conversion 0.8, a lattice run a few seconds, so the 30
LAMMPS and 60 lattice runs take a few minutes. Running on to conversion 0.96 more than doubles the
LAMMPS time (`--stop 0.96`); more samples (`--samples 15`) narrow the error bars.

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

The pipeline has been tested at the smoke size; the results of the full comparison will be added
here.
