#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2026 Zidan Zhang <zhangzidan@gmail.com>
"""Tests of fix poly/react in LAMMPS, on each given number of processors.

    run_tests.py --lmp LMP --plugin polyreactplugin.so --preset polyreact_preset
                 [--mpirun mpirun] [--np 1 4] [--workdir DIR]

With --mpirun "" (a serial LAMMPS) only 1 processor is used.

  pairs     isolated A-B pairs: bonds only between partners, with the new types; the probability per
            check of the pair reaction and of the competing decay B -> D (5 seeds, both partner modes)
  triplets  B-A-B: A takes one bond per check up to its functionality and closes the angle B-A-B;
            the first bonds join two monomers, the second a polymer and a monomer
  squares   two bonded dimers: the second bond of a square closes a ring, refused unless the reaction
            is intramolecular, also when both bonds are proposed in the same check
  kg        Kremer-Grest melt with the FRP preset: bonds, angles and pair, bond and angle energies
            equal those of the final configuration read again, with special lists built by LAMMPS
            (1-2 neighbors only, and 1-2 to 1-4 neighbors with newton_bond on and off); the
            molecule statistics of the fix equal those of the bonds in the data file
"""
import argparse
import math
import pathlib
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

HERE = pathlib.Path(__file__).resolve().parent
N = 10  # systems of N^3 pairs, triplets or squares, 3 apart in a box of 3N


# ---- data files ---------------------------------------------------------------------------------

def write_data(path, atoms, bonds=()):
    """atoms: (molecule, type, x, y, z); bonds: (type, a, b)."""
    with open(path, "w") as f:
        f.write(f"fix poly/react test\n\n{len(atoms)} atoms\n4 atom types\n{len(bonds)} bonds\n1 bond types\n")
        f.write(f"0 angles\n1 angle types\n\n0 {3 * N} xlo xhi\n0 {3 * N} ylo yhi\n0 {3 * N} zlo zhi\n\n")
        f.write("Masses\n\n1 1.0\n2 1.0\n3 1.0\n4 1.0\n\nAtoms # full\n\n")
        for i, (mol, t, x, y, z) in enumerate(atoms, 1):
            f.write(f"{i} {mol} {t} 0.0 {x:.4f} {y:.4f} {z:.4f}\n")
        if bonds:
            f.write("\nBonds\n\n" + "".join(f"{k} {t} {a} {b}\n" for k, (t, a, b) in enumerate(bonds, 1)))


def sites():
    for i in range(N):
        for j in range(N):
            for k in range(N):
                yield 3 * i, 3 * j, 3 * k + 0.5


def wrap(x):
    return x % (3 * N)


def make_pairs(path):
    """A at odd IDs, its partner B at the next ID, 0.5 apart; pairs straddle the periodic and the
    processor boundaries."""
    atoms = []
    for x, y, z in sites():
        atoms.append((len(atoms) + 1, 1, wrap(x - 0.25), y + 0.5, z))
        atoms.append((len(atoms) + 1, 2, wrap(x + 0.25), y + 0.5, z))
    write_data(path, atoms)


def make_triplets(path):
    atoms = []
    for x, y, z in sites():
        atoms += [(len(atoms) + 1, 2, x + 0.5, y + 0.5, z), (len(atoms) + 2, 1, x + 1.0, y + 0.5, z),
                  (len(atoms) + 3, 2, x + 1.5, y + 0.5, z)]
    write_data(path, atoms)


def make_squares(path):
    """Dimers 1-2 and 3-4, 0.9 apart: 1-3 and 2-4 may bond, the diagonals are out of reach."""
    atoms, bonds = [], []
    for x, y, z in sites():
        a = len(atoms) + 1
        x, y = wrap(x - 0.45), wrap(y - 0.45)
        atoms += [(a, 1, x, y, z), (a, 1, wrap(x + 0.9), y, z), (a + 2, 1, x, wrap(y + 0.9), z),
                  (a + 2, 1, wrap(x + 0.9), wrap(y + 0.9), z)]
        bonds += [(1, a, a + 1), (1, a + 2, a + 3)]
    write_data(path, atoms, bonds)


def read_data(path):
    """Atoms {id: type label}, bonds [(a, b)] and angles [(a, b, c)] of a file from write_data."""
    types, labels, bonds, angles = {}, {}, [], []
    section = "header"
    with open(path) as f:
        next(f)
        for line in f:
            s = line.split("#")[0].strip()
            if not s:
                continue
            w = s.split()
            if w[0][0].isalpha():
                section = s
            elif section == "Atom Type Labels":
                labels[w[0]] = w[1]
            elif section == "Atoms":
                types[int(w[0])] = w[2]
            elif section == "Bonds":
                bonds.append((int(w[2]), int(w[3])))
            elif section == "Angles":
                angles.append((int(w[2]), int(w[3]), int(w[4])))
    return {i: labels.get(t, t) for i, t in types.items()}, bonds, angles


def thermo(log):
    """Rows of the last thermo table of a LAMMPS log, as dicts."""
    lines = pathlib.Path(log).read_text().splitlines()
    start = max(i for i, l in enumerate(lines) if l.split()[:1] == ["Step"])
    keys = lines[start].split()
    rows = []
    for l in lines[start + 1:]:
        w = l.split()
        if len(w) != len(keys) or not w[0].isdigit():
            break
        rows.append(dict(zip(keys, map(float, w))))
    return rows


# ---- checks -------------------------------------------------------------------------------------

class Failure(Exception):
    pass


def require(condition, message):
    if not condition:
        raise Failure(message)


def check_thermo(log, expected):
    """The last thermo row has these values, e.g. {"f_r[4]": 1000}."""
    row = thermo(log)[-1]
    for key, value in expected.items():
        require(math.isclose(row[key], value, rel_tol=1e-9), f"{key} = {row[key]:g}, expected {value:g}")


def molecules(nbeads, bonds):
    """Number, Mn, Mw, largest and Mw without the largest of the molecules of 2 or more beads."""
    parent = list(range(nbeads + 1))

    def root(t):
        while parent[t] != t:
            parent[t] = parent[parent[t]]
            t = parent[t]
        return t

    for a, b in bonds:
        ra, rb = root(a), root(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
    size = defaultdict(int)
    for t in range(1, nbeads + 1):
        size[root(t)] += 1
    sizes = [s for s in size.values() if s >= 2]
    n, s1, s2, top = len(sizes), sum(sizes), sum(s * s for s in sizes), max(sizes, default=0)
    return {"f_r[8]": n, "f_r[9]": s1 / n if n else 0, "f_r[10]": s2 / s1 if s1 else 0, "f_r[11]": top,
            "f_r[12]": (s2 - top * top) / (s1 - top) if s1 > top else 0}


def check_pairs(path):
    types, bonds, _ = read_data(path)
    bonded = set()
    for a, b in bonds:
        a, b = min(a, b), max(a, b)
        require(a % 2 == 1 and b == a + 1, f"bond {a}-{b} between atoms of different pairs")
        require(types[a] == "C" and types[b] == "D", f"bond {a}-{b} between {types[a]} and {types[b]}")
        bonded |= {a, b}
    for i, t in types.items():
        if i not in bonded:
            require(t in (("A",) if i % 2 else ("B", "D")), f"unbonded atom {i} has type {t}")
    return f"{len(bonds)} bonds"


def check_triplets(path):
    types, bonds, angles = read_data(path)
    neighbors = defaultdict(set)
    for a, b in bonds:
        neighbors[a].add(b)
        neighbors[b].add(a)
    for i, t in types.items():
        if i % 3 == 2:
            require(t == "A" and neighbors[i] == {i - 1, i + 1}, f"center {i}: {t}, bonded to {neighbors[i]}")
        else:
            require(t == "D" and len(neighbors[i]) == 1, f"end {i}: {t}, {len(neighbors[i])} bonds")
    found = {(min(a, c), b, max(a, c)) for a, b, c in angles}
    require(len(angles) == len(found) == N ** 3, f"{len(angles)} angles")
    require(all(b % 3 == 2 and (a, c) == (b - 1, b + 1) for a, b, c in found), "angle not B-A-B")
    return f"{len(bonds)} bonds, {len(angles)} angles"


def check_squares(path, intramolecular):
    types, bonds, _ = read_data(path)
    new = [(min(a, b), max(a, b)) for a, b in bonds if not (abs(a - b) == 1 and min(a, b) % 2 == 1)]
    require(len(new) == (2 if intramolecular else 1) * N ** 3, f"{len(new)} new bonds")
    for a, b in new:
        require(b == a + 2 and (a - 1) % 4 in (0, 1), f"bond {a}-{b} outside a square")
        require(types[a] == types[b] == "C", f"bond {a}-{b} between {types[a]} and {types[b]}")
    return f"{len(new)} new bonds"


# ---- runs ---------------------------------------------------------------------------------------

class Runner:
    def __init__(self, args):
        self.args = args
        self.work = pathlib.Path(args.workdir).resolve()
        self.work.mkdir(parents=True, exist_ok=True)
        self.failures = 0

    def lammps(self, np, script, log, **variables):
        command = [self.args.mpirun, "-np", str(np)] if self.args.mpirun else []
        command += [self.args.lmp, "-in", str(HERE / script), "-log", str(log), "-screen", "none",
                    "-var", "plugin", self.args.plugin]
        for key, value in variables.items():
            command += ["-var", key, str(value)]
        result = subprocess.run(command, cwd=self.work, capture_output=True, text=True)
        if result.returncode != 0:
            text = pathlib.Path(log).read_text() if pathlib.Path(log).exists() else result.stderr
            errors = [l for l in text.splitlines() if "ERROR" in l] or text.splitlines()[-5:]
            raise Failure("LAMMPS failed: " + " | ".join(errors))

    def test(self, name, function):
        try:
            print(f"  {name}: ok, {function()}")
        except Failure as e:
            print(f"  {name}: FAILED, {e}")
            self.failures += 1

    def small(self, np, data, rf, out, seed=7, partner="all", steps=10):
        log = self.work / f"{out}.log"
        self.lammps(np, "small.in", log, data=data, rf=HERE / rf, seed=seed, partner=partner,
                    steps=steps, out=self.work / f"{out}.data")
        return self.work / f"{out}.data", log

    def run(self):
        w = self.work
        make_pairs(w / "pairs.data")
        make_triplets(w / "triplets.data")
        make_squares(w / "squares.data")
        subprocess.run([self.args.preset, "frp", "k_init=0.02", "k_p1=0.2", "k_p2=0.2", "k_term=0.5",
                        "interval=20", "cutoff=1.2", "bond_type=1", "angle_type=1"],
                       stdout=open(w / "frp.rf", "w"), check=True)
        for np in self.args.np:
            print(f"{np} processor(s):")
            self.test("pairs", lambda: check_pairs(self.small(np, "pairs.data", "pairs.rf", f"pairs_{np}")[0]))
            for partner in ("all", "random"):
                self.test(f"triplets, partner {partner}", lambda: self.triplets(np, partner))
            for intra in ("no", "yes"):
                self.test(f"squares, intramolecular {intra}", lambda: self.squares(np, intra))
            for partner in ("all", "random"):
                self.test(f"probabilities, partner {partner}", lambda: self.rates(np, partner))
            for special, newton in (("fene", "on"), ("weights", "on"), ("weights", "off")):
                self.test(f"kg, special {special}, newton_bond {newton}", lambda: self.kg(np, special, newton))
        return self.failures

    def triplets(self, np, partner):
        data, log = self.small(np, "triplets.data", "triplets.rf", f"triplets_{np}_{partner}", partner=partner, steps=4)
        m = N ** 3
        check_thermo(log, {"f_r[3]": 2 * m, "f_r[4]": m, "f_r[5]": m, "f_r[6]": 0, "f_r[7]": 0, "f_r[8]": m,
                           "f_r[9]": 3, "f_r[10]": 3, "f_r[11]": 3, "f_r[12]": 3})
        return check_triplets(data) + ", bonds monomer-monomer then polymer-monomer"

    def squares(self, np, intra):
        data, log = self.small(np, "squares.data", f"squares_{intra}.rf", f"squares_{np}_{intra}", steps=3)
        m = N ** 3
        rings = m if intra == "yes" else 0
        check_thermo(log, {"f_r[4]": 0, "f_r[5]": 0, "f_r[6]": m, "f_r[7]": rings, "f_r[8]": m, "f_r[10]": 4})
        return check_squares(data, intra == "yes") + f", {rings} rings"

    def rates(self, np, partner):
        """Per check, an intact pair reacts with p_ab or B decays with p_dec; when both pass, one of
        the two is taken at random."""
        n = ab = dec = 0.0
        for seed in (11, 22, 33, 44, 55):
            _, log = self.small(np, "pairs.data", "pairs.rf", f"rate_{np}_{partner}_{seed}", seed=seed,
                                partner=partner, steps=20)
            rows = thermo(log)
            for prev, cur in zip(rows, rows[1:]):
                n += prev["c_t[2]"]
                ab += cur["f_r[13]"] - prev["f_r[13]"]
                dec += (cur["f_r[1]"] - cur["f_r[13]"]) - (prev["f_r[1]"] - prev["f_r[13]"])
        report = []
        for name, k, p in (("pair", ab, 0.1 * 0.95 + 0.1 * 0.05 / 2), ("decay", dec, 0.05 * 0.9 + 0.1 * 0.05 / 2)):
            z = (k / n - p) / math.sqrt(p * (1 - p) / n)
            require(abs(z) < 4, f"{name} {k / n:.4f}, expected {p:.4f} ({z:+.1f} sigma)")
            report.append(f"{name} {k / n:.4f} (expected {p:.4f}, {z:+.1f} sigma)")
        return ", ".join(report)

    def kg(self, np, special, newton):
        name = f"kg_{np}_{special}_{newton}"
        data = self.work / f"{name}.data"
        self.lammps(np, "kg.in", self.work / f"{name}.log", rf=self.work / "frp.rf", special=special,
                    newton=newton, out=data)
        self.lammps(1, "reload.in", self.work / f"{name}_reload.log", data=data, special=special, newton=newton)
        run, again = thermo(self.work / f"{name}.log")[-1], thermo(self.work / f"{name}_reload.log")[-1]
        require(run["Bonds"] == again["Bonds"] == run["f_r[3]"], f"bonds {run['Bonds']} {again['Bonds']} {run['f_r[3]']}")
        require(run["Angles"] == again["Angles"], f"angles {run['Angles']} {again['Angles']}")
        for key in ("E_vdwl", "E_bond", "E_angle"):
            require(math.isclose(run[key], again[key], rel_tol=1e-10), f"{key} {run[key]!r} != {again[key]!r}")
        require(sum(run[f"f_r[{k}]"] for k in (4, 5, 6, 7)) == run["f_r[3]"], "the bond tally does not add up")
        types, bonds, _ = read_data(data)
        for key, value in molecules(len(types), bonds).items():
            require(math.isclose(run[key], value, rel_tol=1e-9), f"{key} = {run[key]:g}, data file {value:g}")
        return (f"{run['f_r[1]']:.0f} reactions, {run['Bonds']:.0f} bonds, {run['f_r[7]']:.0f} rings, "
                f"{run['f_r[11]']:.0f} in the largest molecule; energies and molecules equal")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--lmp", required=True, help="LAMMPS executable with plugin support")
    parser.add_argument("--plugin", required=True, help="polyreactplugin.so")
    parser.add_argument("--preset", required=True, help="polyreact_preset")
    parser.add_argument("--mpirun", default="mpirun", help="MPI launcher; empty for a serial LAMMPS")
    parser.add_argument("--np", type=int, nargs="+", default=[1, 4], help="numbers of processors")
    parser.add_argument("--workdir", help="directory for the runs (default: a new temporary directory)")
    args = parser.parse_args()
    # the runs start in the work directory: make paths absolute, find bare names on the PATH
    for key in ("lmp", "plugin", "preset", "mpirun"):
        value = getattr(args, key)
        if value and pathlib.Path(value).exists():
            setattr(args, key, str(pathlib.Path(value).resolve()))
        elif value and key != "plugin" and shutil.which(value):
            setattr(args, key, shutil.which(value))
    if not args.mpirun:
        args.np = [1]
    if not args.workdir:
        args.workdir = tempfile.mkdtemp(prefix="polyreact_plugin_tests_")
        print(f"runs in {args.workdir}")
    failures = Runner(args).run()
    print("all tests passed" if failures == 0 else f"{failures} test(s) failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
