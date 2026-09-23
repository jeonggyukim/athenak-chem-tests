#!/usr/bin/env python
"""Compare an AthenaK binary dump with an Athena++ VTK dump, cell by cell.

Both codes must have started from the same field and integrated the same
interval, which is what the shared-VTK protocol arranges. Differences are then
the ODE solver plus whatever the two network implementations disagree about.

Usage:
    compare_ak_pp.py ak.bin pp.vtk [--out summary.txt]
"""
import argparse
import os
import pathlib
import sys

import numpy as np

ATHENAK_ROOT = pathlib.Path(
    os.environ.get("ATHENAK_ROOT", pathlib.Path.home() / "Projects/athenak-chem")
)
sys.path.insert(0, str(ATHENAK_ROOT / "vis/python"))
import bin_convert  # noqa: E402

SPECIES = [
    "He+", "OHx", "CHx", "CO", "C+", "HCO+",
    "H2", "H+", "H3+", "H2+", "O+", "Si+",
]


def read_athenapp_vtk(fname):
    """Read a legacy Athena++ VTK (RECTILINEAR_GRID, CELL_DATA, big-endian
    float32) into {name: array of shape (nk, nj, ni)}."""
    fields = {}
    with open(fname, "rb") as f:
        dims = None
        ncells = None
        while True:
            line = f.readline()
            if not line:
                break
            parts = line.split()
            if not parts:
                continue
            key = parts[0].decode("ascii", "replace")
            if key == "DIMENSIONS":
                dims = [int(p) - 1 for p in parts[1:4]]
            elif key in ("X_COORDINATES", "Y_COORDINATES", "Z_COORDINATES"):
                n = int(parts[1])
                np.frombuffer(f.read(4 * n), dtype=">f4")
                f.readline()
            elif key == "CELL_DATA":
                ncells = int(parts[1])
            elif key in ("SCALARS", "VECTORS"):
                name = parts[1].decode("ascii", "replace")
                ncomp = 3 if key == "VECTORS" else 1
                if key == "SCALARS":
                    f.readline()  # LOOKUP_TABLE line
                raw = np.frombuffer(f.read(4 * ncomp * ncells), dtype=">f4")
                ni, nj, nk = dims
                if ncomp == 1:
                    fields[name] = raw.reshape(nk, nj, ni).astype(np.float64)
                else:
                    vec = raw.reshape(nk, nj, ni, 3).astype(np.float64)
                    for c, suffix in enumerate(("1", "2", "3")):
                        fields[name + suffix] = vec[..., c]
                f.readline()
    return fields


def stats(ak, pp):
    """Relative difference of AthenaK against Athena++, per cell."""
    denom = np.maximum(np.abs(pp), 1e-30)
    rel = np.abs(ak - pp) / denom
    return np.median(rel), np.percentile(rel, 90), rel.max()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("akfile")
    parser.add_argument("ppfile")
    parser.add_argument("--out")
    args = parser.parse_args()

    ak = bin_convert.read_binary_as_athdf(args.akfile)
    pp = read_athenapp_vtk(args.ppfile)

    gamma = 5.0 / 3.0
    ak_rho = np.asarray(ak["dens"], dtype=np.float64)
    ak_press = np.asarray(ak["eint"], dtype=np.float64) * (gamma - 1.0)
    pp_rho, pp_press = pp["rho"], pp["press"]

    lines = []
    lines.append(f"AthenaK : {args.akfile}  t={float(ak['Time']):.6g}")
    lines.append(f"Athena++: {args.ppfile}")
    lines.append(f"cells   : {ak_rho.size}")
    lines.append("")
    lines.append(f"{'quantity':10s} {'median':>10s} {'90th pct':>10s} {'max':>10s}"
                 f"   {'AthenaK range':>26s}")

    def row(name, a, p):
        med, p90, mx = stats(a, p)
        lines.append(
            f"{name:10s} {med:10.2e} {p90:10.2e} {mx:10.2e}"
            f"   [{a.min():11.4e},{a.max():11.4e}]"
        )

    row("density", ak_rho, pp_rho)
    row("pressure", ak_press, pp_press)
    for i, s in enumerate(SPECIES):
        row(s, np.asarray(ak[f"s_{i:02d}_chem_{s}"], dtype=np.float64), pp["r" + s])

    # Temperature proxy: pressure per particle tracks T when the abundances
    # agree, so a large split between this and the pressure row means the two
    # codes disagree about the composition rather than about the energy.
    lines.append("")
    med, p90, mx = stats(ak_press / ak_rho, pp_press / pp_rho)
    lines.append(f"{'P/rho':10s} {med:10.2e} {p90:10.2e} {mx:10.2e}")

    text = "\n".join(lines)
    print(text)
    if args.out:
        pathlib.Path(args.out).write_text(text + "\n")


if __name__ == "__main__":
    main()
