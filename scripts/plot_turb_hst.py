#!/usr/bin/env python
"""History of the decaying-turbulence runs from the two codes' .hst files.

Both codes write the same history header format (`[n]=name`), with volume
integrals over the box: mass, kinetic energy per direction and total energy.
Thermal energy is total minus kinetic, since neither run has magnetic fields.

Usage:
    plot_turb_hst.py --ak ak.hydro.hst --ap ap.hst --out fig_turb_hst.png
"""
import argparse
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

MYR_PER_CODE = 0.97784  # pc/(km/s) expressed in Myr


def read_hst(fname):
    with open(fname) as f:
        for line in f:
            if "[1]=" in line:
                names = [m.group(1) for m in re.finditer(r"\[\d+\]=(\S+)", line)]
                break
    data = np.loadtxt(fname, ndmin=2)
    cols = {name: data[:, i] for i, name in enumerate(names)}
    ke = cols["1-KE"] + cols["2-KE"] + cols["3-KE"]
    return {
        "t": cols["time"] * MYR_PER_CODE,
        "dt": cols["dt"] * MYR_PER_CODE,
        "mass": cols["mass"],
        "KE": ke,
        "TE": cols["tot-E"] - ke,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ak", required=True, help="AthenaK .hydro.hst")
    parser.add_argument("--ap", required=True, help="Athena++ .hst")
    parser.add_argument("--out", default="fig_turb_hst.png")
    args = parser.parse_args()

    runs = [("Athena++ CVODE", read_hst(args.ap), "C0"),
            ("AthenaK semi-implicit", read_hst(args.ak), "C1")]
    panels = [("KE", "kinetic energy [code]"),
              ("TE", "thermal energy [code]"),
              ("dt", "dt [Myr]"),
              ("mass", "mass [code]")]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for ax, (key, label) in zip(axes.flat, panels):
        for name, h, c in runs:
            ax.plot(h["t"], h[key], color=c, label=name)
        ax.set_ylabel(label)
        if key in ("KE", "TE"):
            ax.set_yscale("log")
    for ax in axes[1]:
        ax.set_xlabel("t [Myr]")
    axes[0, 0].legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=140)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
