#!/usr/bin/env python3
# SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
# Copyright (C) 2023-2026 Zidan Zhang <zhangzidan@gmail.com>
"""Analysis of the model comparison (see README.md). Needs numpy and scipy.

    analyze.py lattice DIR      PSF snapshots of a 3dbfm or statistical run -> DIR/clusters_<s>.csv,
                                then deletes the snapshots (--keep keeps them)
    analyze.py compare [ROOT]   tables for all runs below ROOT (bfm/, stat/, md/), in Markdown

For every model the same quantities: molecules of two or more units (number, Mw, largest, and RDP,
Mw without the largest), the gel point as the vinyl conversion at the maximum of RDP, and the
classification of new bonds (monomer-monomer, polymer-monomer, polymer-polymer, intramolecular).
"""
import glob
import os
import re
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

X = (4, 8, 12, 16, 20)
MODELS = (("bfm", "3dbfm"), ("stat", "statistical"), ("md", "LAMMPS"), ("md_paper", "LAMMPS, rates of the paper"))


# ---- lattice runs: PSF snapshots -> cluster tables --------------------------------------------

def psf_sizes(path):
    with open(path) as f:
        lines = f.read().split("\n")
    natom = int(next(l for l in lines if "!NATOM" in l).split()[0])
    j = next(k for k, l in enumerate(lines) if "!NBOND" in l)
    nb = int(lines[j].split()[0])
    nums = np.array(" ".join(lines[j + 1:]).replace("END", " ").split()[: 2 * nb], dtype=np.int64) - 1
    g = coo_matrix((np.ones(nb), (nums[0::2], nums[1::2])), shape=(natom, natom))
    return np.bincount(connected_components(g, directed=False)[1])


def molecule_stats(sizes):
    s = sizes[sizes >= 2].astype(float)
    if len(s) == 0:
        return 0, 0.0, 0.0, 0.0
    s1, s2, top = s.sum(), (s * s).sum(), s.max()
    return len(s), s2 / s1, top, (s2 - top * top) / (s1 - top) if s1 > top else 0.0


def lattice(d, keep=False):
    for prop in sorted(glob.glob(os.path.join(d, "propagation_*.dat"))):
        sample = re.search(r"propagation_(\d+)\.dat$", prop).group(1)
        tally = np.loadtxt(prop, comments="#", ndmin=2)  # conversion (level or actual), mm, pm, pp, intra
        out = os.path.join(d, f"clusters_{sample}.csv")
        rows = {}
        for psf in sorted(glob.glob(os.path.join(d, f"config_{sample}_*.psf"))):
            level = float(re.search(r"config_\d+_([\d.]+)\.psf$", psf).group(1))
            n, mw, top, rdp = molecule_stats(psf_sizes(psf))
            t = tally[int(np.argmin(np.abs(tally[:, 0] - level))), 1:5].astype(int)
            rows[level] = f"{level:.4f},{n},{mw:.4f},{top:.0f},{rdp:.4f},{t[0]},{t[1]},{t[2]},{t[3]}\n"
            if not keep:
                os.remove(psf)
        with open(out, "w") as f:
            f.write("# conversion,molecules,Mw,largest,RDP,mm,pm,pp,intra\n")
            f.writelines(rows[level] for level in sorted(rows))


# ---- reading runs -------------------------------------------------------------------------------

def read_csv(path):
    a = np.loadtxt(path, delimiter=",", comments="#", ndmin=2)
    return dict(zip(("cv", "molecules", "Mw", "largest", "rdp", "mm", "pm", "pp", "intra"), a.T))


def read_log(path):
    """The thermo table of the reaction run of in.radical."""
    lines = open(path).read().splitlines()
    starts = [i for i, l in enumerate(lines) if l.split()[:3] == ["Step", "Time", "v_cv"]]
    if not starts:
        return None
    keys = lines[starts[-1]].split()
    rows = []
    for l in lines[starts[-1] + 1:]:
        w = l.split()
        if len(w) == len(keys) and w[0].isdigit():
            rows.append([float(v) for v in w])
        elif not l.startswith("Fix halt"):
            break
    if not rows:
        return None
    c = dict(zip(keys, np.array(rows).T))
    names = {"cv": "v_cv", "time": "Time", "mm": "f_react[4]", "pm": "f_react[5]", "pp": "f_react[6]",
             "intra": "f_react[7]", "molecules": "f_react[8]", "Mw": "f_react[10]", "largest": "f_react[11]",
             "rdp": "f_react[12]"}
    return {k: c[v] for k, v in names.items()}


def runs(root, model, case):
    if model.startswith("md"):
        found = [read_log(p) for p in sorted(glob.glob(os.path.join(root, model, f"{case}_s*", "log_*.lammps")))]
    else:
        found = [read_csv(p) for p in sorted(glob.glob(os.path.join(root, model, f"{case}_s*", "clusters_*.csv")))]
    return [d for d in found if d is not None and len(d["cv"]) > 1]


def at(d, key, cv):
    """d[key] at the first record with conversion >= cv (None if the run did not get there)."""
    i = int(np.searchsorted(d["cv"], cv - 1e-9))
    return None if i >= len(d["cv"]) else d[key][i]


def fraction(d, key, cv):
    bonds = [at(d, k, cv) for k in ("mm", "pm", "pp", "intra")]
    return None if bonds[0] is None else at(d, key, cv) / max(sum(bonds), 1)


class Censored:
    """A run that ended before the quantity could be found (e.g. RDP still rising at the end)."""


def cell(values, digits):
    censored = sum(isinstance(x, Censored) for x in values)
    v = [x for x in values if x is not None and not isinstance(x, Censored)]
    note = f", {censored} ended before" if censored else ""
    if not v:
        return f"— ({censored} ended before)" if censored else "—"
    if len(v) == 1:
        return f"{v[0]:.{digits}f}" + (f" ({note[2:]})" if note else "")
    return f"{np.mean(v):.{digits}f} ± {np.std(v, ddof=1) / np.sqrt(len(v)):.{digits}f} ({len(v)}{note})"


def gel_point(d):
    """Conversion at the maximum of RDP; Censored if the maximum is at the last record (the run
    ended before RDP went down again), None if the run is too short to tell."""
    if d["cv"][-1] < 0.5:
        return None
    i = int(np.argmax(d["rdp"]))
    return Censored() if i >= len(d["rdp"]) - 1 else d["cv"][i]


def table(root, mechanism, title, value, digits):
    """One row per cross-linker fraction with runs, one column per model with runs."""
    models = [(m, name) for m, name in MODELS if any(runs(root, m, f"{mechanism}_X{x:02d}") for x in X)]
    rows = [x for x in X if any(runs(root, m, f"{mechanism}_X{x:02d}") for m, _ in models)]
    if not rows:
        return
    print(f"\n{mechanism.upper()}: {title}\n")
    print("| [X] | " + " | ".join(name for _, name in models) + " |")
    print("|:--|" + "--:|" * len(models))
    for x in rows:
        cells = [cell([value(d) for d in runs(root, m, f"{mechanism}_X{x:02d}")], digits) for m, _ in models]
        print(f"| {x} % | " + " | ".join(cells) + " |")


def compare(root):
    print("Mean ± standard error (number of runs).")
    for m in ("atrp", "frp"):
        table(root, m, "gel point, vinyl conversion at the maximum of RDP", gel_point, 3)
        table(root, m, "molecules of 2 or more units at conversion 0.05", lambda d: at(d, "molecules", 0.05), 0)
        for cv in (0.2, 0.8):
            table(root, m, f"intramolecular bonds / all new bonds at conversion {cv}", lambda d, cv=cv: fraction(d, "intra", cv), 4)
        table(root, m, "polymer-polymer bonds / all new bonds at conversion 0.8", lambda d: fraction(d, "pp", 0.8), 4)


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "lattice":
        lattice(sys.argv[2], keep="--keep" in sys.argv)
    elif len(sys.argv) >= 2 and sys.argv[1] == "compare":
        compare(sys.argv[2] if len(sys.argv) > 2 else ".")
    else:
        sys.exit(__doc__)
