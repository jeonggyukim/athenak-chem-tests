#!/usr/bin/env python
"""Volume-integrated history of the decaying-turbulence runs, from the dumps.

Athena++ per-MeshBlock VTK files are summed directly, so a run can be plotted
while it is still going and before its blocks are joined. The timestep comes
from the AthenaK .hst (column 2) and from the Athena++ `cycle=` lines in stdout.

Usage:
    plot_turb_history_dumps.py \
        --ak-si-dir turb64_si_run --ak-si-base turb64_si \
        --ak-bdf-dir turb64_bdf_run --ak-bdf-base turb64_bdf \
        --ap-cvode-dir turb64_cvode_run --ap-cvode-base turb64_cvode --out fig.png

The Kokkos BDF run is optional.
"""
import argparse
import glob
import os
import pathlib
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ATHENAK_ROOT = pathlib.Path(
    os.environ.get("ATHENAK_ROOT", pathlib.Path.home() / "Projects/athenak-chem")
)
sys.path.insert(0, str(ATHENAK_ROOT / "vis/python"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import bin_convert  # noqa: E402
from compare_ak_ap import read_athenapp_vtk  # noqa: E402

GAMMA = 5.0 / 3.0
MYR_PER_CODE = 0.97784
ERG_PER_CODE = 6.882615081124725e31 * 1.0e10  # [mass] (km/s)^2
XC_TOT = 1.6e-4


def integrals(t, dV, rho, v1, v2, v3, press, xH2, xCO):
    m = rho * dV
    M = m.sum()
    return {
        "t": t * MYR_PER_CODE,
        "mass": M,
        "KE": 0.5 * np.sum(m * (v1**2 + v2**2 + v3**2)) * ERG_PER_CODE,
        "Eth": np.sum(press) * dV / (GAMMA - 1.0) * ERG_PER_CODE,
        "vrms": np.sqrt(np.sum(m * (v1**2 + v2**2 + v3**2)) / M),
        "H2": np.sum(m * 2.0 * xH2) / M,
        "CO": np.sum(m * xCO) / M / XC_TOT,
    }


def athenak_history(ak_dir, base):
    rows = []
    for f in sorted(glob.glob(f"{ak_dir}/bin/{base}.hydro_w.*.bin")):
        d = bin_convert.read_binary_as_athdf(f)
        dV = np.prod([np.diff(d[k])[0] for k in ("x1f", "x2f", "x3f")])
        g = {k: np.asarray(d[k], dtype=np.float64) for k in
             ("dens", "velx", "vely", "velz", "eint", "s_06_chem_H2", "s_03_chem_CO")}
        rows.append(integrals(float(d["Time"]), dV, g["dens"], g["velx"], g["vely"],
                              g["velz"], g["eint"] * (GAMMA - 1.0),
                              g["s_06_chem_H2"], g["s_03_chem_CO"]))
    hst = np.loadtxt(glob.glob(f"{ak_dir}/{base}.hydro.hst")[0])
    return rows, hst[:, 0] * MYR_PER_CODE, hst[:, 1]


def vtk_time_and_spacing(fname):
    """Time and uniform cell width from a RECTILINEAR_GRID VTK header."""
    with open(fname, "rb") as fh:
        fh.readline()
        t = float(re.search(r"time=(\S+)", fh.readline().decode("ascii")).group(1))
        line = b""
        while not line.startswith(b"X_COORDINATES"):
            line = fh.readline()
        x = np.frombuffer(fh.read(4 * int(line.split()[1])), dtype=">f4")
    return t, float(x[1] - x[0])


def athenapp_history(ap_dir, base):
    rows = []
    steps = sorted({f.rsplit(".", 2)[-2] for f in
                    glob.glob(f"{ap_dir}/{base}.block*.out1.*.vtk")})
    nblock = len(glob.glob(f"{ap_dir}/{base}.block*.out1.00000.vtk"))
    for n in steps:
        files = sorted(glob.glob(f"{ap_dir}/{base}.block*.out1.{n}.vtk"))
        if len(files) < nblock:
            continue
        t, dx = vtk_time_and_spacing(files[0])
        blocks = [read_athenapp_vtk(f) for f in files]
        cat = {k: np.concatenate([b[k].ravel() for b in blocks]).astype(np.float64)
               for k in ("rho", "vel1", "vel2", "vel3", "press", "rH2", "rCO")}
        rows.append(integrals(t, dx**3, cat["rho"], cat["vel1"], cat["vel2"],
                              cat["vel3"], cat["press"], cat["rH2"], cat["rCO"]))
    t_dt = [(float(m.group(1)), float(m.group(2))) for m in
            re.finditer(r"^cycle=\d+ time=(\S+) dt=(\S+)",
                        pathlib.Path(f"{ap_dir}/run.out").read_text(), re.M)]
    t_dt = np.array(t_dt)
    return rows, t_dt[:, 0] * MYR_PER_CODE, t_dt[:, 1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ak-si-dir", required=True)
    parser.add_argument("--ak-si-base", required=True)
    parser.add_argument("--ap-cvode-dir", required=True)
    parser.add_argument("--ap-cvode-base", required=True)
    parser.add_argument("--ak-si-table-dir", help="optional semi-implicit, tabulated")
    parser.add_argument("--ak-si-table-base")
    parser.add_argument("--ak-bdf-dir")
    parser.add_argument("--ak-bdf-base")
    parser.add_argument("--out", default="fig_turb_history_dumps.png")
    args = parser.parse_args()

    ak, ak_t, ak_dt = athenak_history(args.ak_si_dir, args.ak_si_base)
    ap, ap_t, ap_dt = athenapp_history(args.ap_cvode_dir, args.ap_cvode_base)
    bdf = (athenak_history(args.ak_bdf_dir, args.ak_bdf_base)
           if args.ak_bdf_dir else None)
    tab = (athenak_history(args.ak_si_table_dir, args.ak_si_table_base)
           if args.ak_si_table_dir else None)
    col = lambda rows, k: np.array([r[k] for r in rows])  # noqa: E731
    print(f"AthenaK {len(ak)} dumps to {ak[-1]['t']:.3f} Myr; "
          f"Athena++ {len(ap)} dumps to {ap[-1]['t']:.3f} Myr")

    panels = [
        ("KE", "kinetic energy [erg]", "log"),
        ("Eth", "thermal energy [erg]", "log"),
        ("vrms", r"mass-weighted $v_{\rm rms}$ [km s$^{-1}$]", "linear"),
        ("H2", r"mass-weighted $2x({\rm H_2})$", "log"),
        ("CO", r"mass-weighted $x({\rm CO})/x_{\rm C,tot}$", "log"),
        ("dt", r"$\Delta t$ [Myr]", "log"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for ax, (key, label, scale) in zip(axes.flat, panels):
        if key == "dt":
            ax.plot(ap_t, ap_dt * MYR_PER_CODE, "--", color="#2E8B57", lw=1.6)
            ax.plot(ak_t, ak_dt * MYR_PER_CODE, "-", color="#C1440E", lw=4.0,
                    alpha=0.4)
            if tab:
                ax.plot(tab[1], tab[2] * MYR_PER_CODE, "-.", color="#7A2A06", lw=1.4)
            if bdf:
                ax.plot(bdf[1], bdf[2] * MYR_PER_CODE, ":", color="#3B5BA9", lw=1.6)
        else:
            ax.plot(col(ap, "t"), col(ap, key), "--", color="#2E8B57", lw=1.8,
                    label="Athena++ CVODE")
            ax.plot(col(ak, "t"), col(ak, key), "-", color="#C1440E", lw=4.0,
                    alpha=0.4, label="AthenaK semi-implicit")
            if tab:
                ax.plot(col(tab[0], "t"), col(tab[0], key), "-.", color="#7A2A06",
                        lw=1.4, label="AthenaK semi-implicit, tabulated")
            if bdf:
                ax.plot(col(bdf[0], "t"), col(bdf[0], key), ":", color="#3B5BA9",
                        lw=1.8, label="AthenaK Kokkos BDF")
        ax.set(xlabel="t [Myr]", ylabel=label, yscale=scale,
               xlim=(0, max(ak[-1]["t"], ap[-1]["t"])))
        ax.grid(alpha=0.25, lw=0.5)
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle("Decaying turbulence from a shared snapshot: volume-integrated history")
    fig.savefig(args.out, dpi=140)
    print("wrote", args.out)
    runs = ([("semi-implicit", ak)] + ([("si tabulated", tab[0])] if tab else [])
            + ([("Kokkos BDF", bdf[0])] if bdf else []))
    for name, rows in runs:
        for key in ("mass", "KE", "Eth", "H2", "CO"):
            a, p = col(rows, key), col(ap, key)
            n = min(len(a), len(p))
            rel = (a[n - 1] - p[n - 1]) / p[n - 1]
            print(f"  {name:13s} {key:4s} at t={ap[n-1]['t']:.3f} Myr: "
                  f"AthenaK/Athena++ - 1 = {rel:+.3e}")


if __name__ == "__main__":
    main()
