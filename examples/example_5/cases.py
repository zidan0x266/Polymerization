#!/usr/bin/env python3
# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
"""Job lists of the model comparison (see README.md): one line per run, a bash command that
run_job.sh executes in the run root.

    cases.py lattice [options] > jobs_lattice.txt    3dbfm and statistical runs (one core each)
    cases.py md [options] > jobs_md.txt              LAMMPS runs (MPI)

Options: --size paper|reduced|mini|smoke, --samples N, --mechanisms atrp frp, --x 4 8 12 16 20
(percent cross-linker), --init (initiator fraction), --k-init (FRP initiation), --stop,
--models bfm stat (lattice only), --rates fast|paper (md only: ATRP with ten times the rate constants
of the paper, the default, or with those of the paper, in md_paper/).

The paper and reduced sizes have the composition and rate constants of the paper (0.2 % initiator,
k_init = 1e-5) and run up to conversion 0.96. A mini system holds too few chains at 0.2 % initiator,
so mini and smoke default to 2 % initiator and k_init = 1e-4 (k_init / k_p1 = 2e-4). FRP mini runs
stop at conversion 0.8, since their gel points lie far below and the last vinyls take most of the
time; ATRP runs go to 0.96, since with short chains the gel point can come close to 0.8. LAMMPS
writes its output every 100 steps in mini and smoke runs (every 1000 otherwise), so that conversion
advances by less than 0.01 between lines, as between the snapshots of the lattice runs.
"""
import argparse

SIZES = {  # 3dbfm lattice edge and units (density 0.85), statistical units, LAMMPS FCC cells (4 cells^3 beads)
    "paper": (100, 106250, 106250, 30),
    "reduced": (68, 33408, 33408, 20),
    "mini": (34, 4176, 4176, 10),
    "smoke": (20, 850, 850, 6),
}
SMALL = ("mini", "smoke")


def seed(model, mechanism, x, sample):
    return 100000 * (1 + ("bfm", "stat", "md").index(model)) + 1000 * x + 100 * ("atrp", "frp").index(mechanism) + sample


def rates(mechanism, k_init, factor=1.0):
    """Rate constants of the paper (Table 1), ATRP propagation scaled by factor."""
    if mechanism == "atrp":
        return f"k_p1={0.005 * factor:g} k_p2={0.0025 * factor:g}"
    return f"k_init={k_init:g} k_p1=0.5 k_p2=0.25 k_term=0.95"


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("kind", choices=("lattice", "md"))
    p.add_argument("--size", choices=SIZES, default="paper")
    p.add_argument("--samples", type=int, default=5)
    p.add_argument("--mechanisms", nargs="+", default=["atrp", "frp"])
    p.add_argument("--x", type=int, nargs="+", default=[4, 8, 12, 16, 20])
    p.add_argument("--init", type=float, help="initiator fraction (default 0.002, mini and smoke 0.02)")
    p.add_argument("--k-init", type=float, help="FRP initiation (default 1e-5, mini and smoke 1e-4)")
    p.add_argument("--stop", type=float, help="conversion at which runs end (default 0.96; mini and smoke FRP 0.8)")
    p.add_argument("--models", nargs="+", default=["bfm", "stat"], help="lattice only")
    p.add_argument("--rates", choices=("fast", "paper"), default="fast", help="md only")
    a = p.parse_args()
    edge, n_bfm, n_stat, cells = SIZES[a.size]
    init = a.init if a.init is not None else (0.02 if a.size in SMALL else 0.002)
    k_init = a.k_init if a.k_init is not None else (1e-4 if a.size in SMALL else 1e-5)
    for m in a.mechanisms:
        stop = a.stop if a.stop is not None else (0.8 if a.size in SMALL and m == "frp" else 0.96)
        for x in a.x:
            for s in range(1, a.samples + 1):
                case = f"{m}_X{x:02d}_s{s}"
                if a.kind == "md":
                    d = "md" if a.rates == "fast" else "md_paper"
                    k = rates(m, k_init, 10.0 if a.rates == "fast" else 1.0)
                    print(f"mkdir -p {d}/{case} && cd {d}/{case} && $POLY_BUILD/src/lammps/polyreact_preset {m} {k} "
                          f"interval=20 cutoff=1.2 bond_type=1 angle_type=1 > react.cfg && $MPIRUN $LMP "
                          f"-in $EXAMPLE/in.radical -var plugin $PLUGIN -var rf react.cfg -var mechanism {m} "
                          f"-var cells {cells} -var xl {x / 100:.2f} -var init {init:g} -var seed {seed('md', m, x, s)} "
                          f"-var stop {stop} -var out {100 if a.size in SMALL else 1000} -var tag {case} -screen none")
                    continue
                common = (f"radical_fraction={init:g} monovinyl_fraction={1 - x / 100:.2f} samples=1 first_sample={s} "
                          f"stop_conversion={stop} snapshot_start=0.01 snapshot_step=0.01 directory=.")
                if "bfm" in a.models:
                    print(f"mkdir -p bfm/{case} && cd bfm/{case} && $POLY_BUILD/src/3dbfm/bfm_react "
                          f"$EXAMPLE/../example_1/{m}_paper.ini 'lattice={edge} {edge} {edge}' n_units={n_bfm} "
                          f"{rates(m, k_init)} {common} seed={seed('bfm', m, x, s)} write_pdb=false "
                          f"> run.log 2>&1 && python3 $EXAMPLE/analyze.py lattice .")
                if "stat" in a.models:
                    print(f"mkdir -p stat/{case} && cd stat/{case} && $POLY_BUILD/src/statistical/stat_react /dev/null "
                          f"n_units={n_stat} mechanism={m} {rates(m, k_init)} cyclization=true {common} "
                          f"seed={seed('stat', m, x, s)} > run.log 2>&1 && python3 $EXAMPLE/analyze.py lattice .")


if __name__ == "__main__":
    main()
