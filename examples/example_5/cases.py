#!/usr/bin/env python3
# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
"""Job lists of the model comparison (see README.md): one line per run, a bash command that
run_job.sh executes in the run root.

    cases.py lattice [options] > jobs_lattice.txt    3dbfm and statistical runs (one core each)
    cases.py md [options] > jobs_md.txt              LAMMPS runs (MPI)

Options: --size paper|reduced|smoke, --samples N, --mechanisms atrp frp, --x 4 8 12 16 20 (percent
cross-linker), --stop 0.96, --models bfm stat (lattice only), --rates fast|paper (md only: ATRP with
ten times the rate constants of the paper, the default, or with those of the paper, in md_paper/).
"""
import argparse

SIZES = {  # 3dbfm lattice edge and units (density 0.85), statistical units, LAMMPS FCC cells (4 cells^3 beads)
    "paper": (100, 106250, 106250, 30),
    "reduced": (68, 33408, 33408, 20),
    "smoke": (20, 850, 850, 6),
}
RATES = {  # Zhang et al., Macromolecules 2025, Table 1
    "atrp": "mechanism=atrp k_p1=0.005 k_p2=0.0025",
    "frp": "mechanism=frp k_init=0.00001 k_p1=0.5 k_p2=0.25 k_term=0.95",
}


def seed(model, mechanism, x, sample):
    return 100000 * (1 + ("bfm", "stat", "md").index(model)) + 1000 * x + 100 * ("atrp", "frp").index(mechanism) + sample


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("kind", choices=("lattice", "md"))
    p.add_argument("--size", choices=SIZES, default="paper")
    p.add_argument("--samples", type=int, default=5)
    p.add_argument("--mechanisms", nargs="+", default=["atrp", "frp"])
    p.add_argument("--x", type=int, nargs="+", default=[4, 8, 12, 16, 20])
    p.add_argument("--stop", type=float, default=0.96)
    p.add_argument("--models", nargs="+", default=["bfm", "stat"], help="lattice only")
    p.add_argument("--rates", choices=("fast", "paper"), default="fast", help="md only")
    a = p.parse_args()
    edge, n_bfm, n_stat, cells = SIZES[a.size]
    for m in a.mechanisms:
        for x in a.x:
            for s in range(1, a.samples + 1):
                case = f"{m}_X{x:02d}_s{s}"
                if a.kind == "md":
                    d = "md" if a.rates == "fast" else "md_paper"
                    rf = f"{m}_paper.cfg" if a.rates == "paper" and m == "atrp" else f"{m}.cfg"
                    print(f"mkdir -p {d}/{case} && cd {d}/{case} && $MPIRUN $LMP -in $EXAMPLE/in.radical "
                          f"-var plugin $PLUGIN -var rf $EXAMPLE/{rf} -var mechanism {m} -var cells {cells} "
                          f"-var xl {x / 100:.2f} -var seed {seed('md', m, x, s)} -var stop {a.stop} -var out 1000 "
                          f"-var tag {case} -screen none")
                    continue
                common = (f"radical_fraction=0.002 monovinyl_fraction={1 - x / 100:.2f} samples=1 first_sample={s} "
                          f"stop_conversion={a.stop} snapshot_start=0.01 snapshot_step=0.01 directory=.")
                if "bfm" in a.models:
                    print(f"mkdir -p bfm/{case} && cd bfm/{case} && $POLY_BUILD/src/3dbfm/bfm_react "
                          f"$EXAMPLE/../example_1/{m}_paper.ini 'lattice={edge} {edge} {edge}' n_units={n_bfm} {common} "
                          f"seed={seed('bfm', m, x, s)} write_pdb=false write_restart=false "
                          f"> run.log 2>&1 && python3 $EXAMPLE/analyze.py lattice .")
                if "stat" in a.models:
                    print(f"mkdir -p stat/{case} && cd stat/{case} && $POLY_BUILD/src/statistical/stat_react /dev/null "
                          f"n_units={n_stat} {RATES[m]} cyclization=true {common} seed={seed('stat', m, x, s)} "
                          f"> run.log 2>&1 && python3 $EXAMPLE/analyze.py lattice .")


if __name__ == "__main__":
    main()
