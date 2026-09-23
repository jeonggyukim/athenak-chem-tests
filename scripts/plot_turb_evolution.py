#!/usr/bin/env python
"""Side-by-side slices, movie and history for the decaying-turbulence runs.

Both codes start from the same VTK snapshot. Panels are laid out with Athena++
(CVODE) on top and AthenaK (semi-implicit) below, sharing one colour scale per
quantity that is fixed across every frame and both codes, so that what changes
in the movie is the gas and not the normalisation.

Usage:
    plot_turb_evolution.py --ak-dir turb_si_run --ak-base turb_si \
        --ap-dir turb_cvode_run --ap-base turb_cvode --outdir figures
"""
import argparse
import glob
import os
import pathlib
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
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
MYR_PER_CODE = 0.97784  # pc/(km/s) expressed in Myr
# Code pressure unit is (1.4 m_H cm^-3)(km/s)^2; this converts it to P/k_B.
PK_PER_CODE = 1.4 * 1.6735575e-24 * 1.0e10 / 1.380649e-16
XC_TOT = 1.6e-4
QUANTITIES = [
    ("nH", r"$n_{\rm H}$ [cm$^{-3}$]", "viridis"),
    ("Pk", r"$P/k_{\rm B}$ [cm$^{-3}$ K]", "inferno"),
    ("H2", r"$2x({\rm H_2})$", "cividis"),
    ("CO", r"$x({\rm CO})/x_{\rm C,tot}$", "magma"),
    ("C+", r"$x({\rm C^+})/x_{\rm C,tot}$", "plasma"),
    ("e", r"$x_{\rm e}$", "viridis"),
    ("H+", r"$x({\rm H^+})$", "cividis"),
    ("H2+", r"$x({\rm H_2^+})$", "magma"),
]
IONS = ["He+", "C+", "HCO+", "H+", "H3+", "H2+", "O+", "Si+"]


def derived(sp, H2, CO):
    """Panel quantities from a {species: abundance} dict."""
    return {
        "H2": 2.0 * H2,
        "CO": CO / XC_TOT,
        "C+": sp["C+"] / XC_TOT,
        "e": sum(sp[s] for s in IONS),
        "H+": sp["H+"],
        "H2+": sp["H2+"],
    }


def athenak_frame(fname):
    d = bin_convert.read_binary_as_athdf(fname)
    nH = np.asarray(d["dens"], dtype=np.float64)
    press = np.asarray(d["eint"], dtype=np.float64) * (GAMMA - 1.0)
    sp = {k.split("_chem_")[1]: np.asarray(v, dtype=np.float64)
          for k, v in d.items() if k.startswith("s_")}
    return {
        "t": float(d["Time"]) * MYR_PER_CODE,
        "nH": nH,
        "Pk": press * PK_PER_CODE,
        **derived(sp, sp["H2"], sp["CO"]),
    }


def athenapp_frame(fname):
    f = read_athenapp_vtk(fname)
    with open(fname, "rb") as fh:
        fh.readline()
        header = fh.readline().decode("ascii", "replace")
    t = np.nan
    for token in header.split():
        if token.startswith("time="):
            t = float(token.split("=")[1]) * MYR_PER_CODE
    return {
        "t": t,
        "nH": f["rho"],
        "Pk": f["press"] * PK_PER_CODE,
        **derived({s: f["r" + s] for s in IONS}, f["rH2"], f["rCO"]),
    }


def fixed_ranges(all_frames):
    """One colour range per quantity, over every frame of both codes."""
    ranges = {}
    for key, _, _ in QUANTITIES:
        # t = 0 holds the uniform r_init abundances, which would set the floor.
        data = np.concatenate([np.asarray(fr[key]).ravel() for fr in all_frames
                               if fr["t"] > 0])
        data = data[np.isfinite(data) & (data > 0)]
        lo, hi = np.percentile(data, [0.5, 99.9])
        ranges[key] = (max(lo, hi * 1e-8), hi)
    return ranges


def history(frames):
    """Mass-weighted means: what a molecular-gas observation would average."""
    out = {k: [] for k in ("t", "H2", "CO", "Pk")}
    for fr in frames:
        w = fr["nH"]
        out["t"].append(fr["t"])
        for key in ("H2", "CO", "Pk"):
            out[key].append(np.sum(w * fr[key]) / np.sum(w))
    return {k: np.array(v) for k, v in out.items()}


def comparison_panel(ap_fr, ak_fr, ranges, outfile, suptitle):
    """Rows alternate Athena++ / AthenaK; each pair of rows holds four quantities."""
    ncol = 4
    nrow = 2 * ((len(QUANTITIES) + ncol - 1) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(17, 4.2 * nrow),
                             constrained_layout=True)
    for q, (key, label, cmap) in enumerate(QUANTITIES):
        lo, hi = ranges[key]
        norm = mcolors.LogNorm(vmin=lo, vmax=hi)
        col = q % ncol
        for code, fr in enumerate((ap_fr, ak_fr)):
            row = 2 * (q // ncol) + code
            ax = axes[row, col]
            k = fr[key].shape[0] // 2
            im = ax.imshow(np.clip(fr[key][k], lo, hi), origin="lower", cmap=cmap,
                           norm=norm, extent=(-16, 16, -16, 16))
            if code == 0:
                ax.set_title(label)
            fig.colorbar(im, ax=ax, fraction=0.046)
    for row in range(nrow):
        name = "Athena++ CVODE" if row % 2 == 0 else "AthenaK semi-implicit"
        axes[row, 0].set_ylabel(f"{name}\ny [pc]")
    for ax in axes[-1]:
        ax.set_xlabel("x [pc]")
    fig.suptitle(suptitle, fontsize=13)
    fig.savefig(outfile, dpi=130)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ak-dir", required=True)
    parser.add_argument("--ak-base", required=True)
    parser.add_argument("--ap-dir", required=True)
    parser.add_argument("--ap-base", required=True)
    parser.add_argument("--outdir", default="figures")
    args = parser.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ak = [athenak_frame(f) for f in
          sorted(glob.glob(f"{args.ak_dir}/bin/{args.ak_base}.hydro_w.*.bin"))]
    ap = [athenapp_frame(f) for f in
          sorted(glob.glob(f"{args.ap_dir}/{args.ap_base}.block0.out1.*.vtk"))]
    print(f"AthenaK {len(ak)} frames to {ak[-1]['t']:.3f} Myr; "
          f"Athena++ {len(ap)} frames to {ap[-1]['t']:.3f} Myr")

    ranges = fixed_ranges(ak + ap)
    for key, _, _ in QUANTITIES:
        print(f"  colour range {key}: {ranges[key][0]:.3e} to {ranges[key][1]:.3e}")

    # Pair each AthenaK frame with the Athena++ frame nearest in time; the two
    # codes choose their own CFL steps, so dump times differ slightly.
    ap_times = np.array([fr["t"] for fr in ap])
    pairs = [(ap[int(np.argmin(np.abs(ap_times - fr["t"])))], fr) for fr in ak]

    framedir = outdir / "compare_frames"
    framedir.mkdir(parents=True, exist_ok=True)
    for i, (ap_fr, ak_fr) in enumerate(pairs):
        comparison_panel(
            ap_fr, ak_fr, ranges, framedir / f"frame_{i:04d}.png",
            f"Decaying turbulence, GOW17 — t = {ak_fr['t']:.3f} Myr",
        )
    movie = outdir / "compare_slices.mp4"
    subprocess.run(
        ["ffmpeg", "-y", "-loglevel", "error", "-framerate", "8",
         "-i", str(framedir / "frame_%04d.png"), "-pix_fmt", "yuv420p",
         "-vf", "scale=1700:-2", str(movie)],
        check=True,
    )
    print("wrote", movie)

    comparison_panel(pairs[-1][0], pairs[-1][1], ranges,
                     outdir / "fig_compare_slices_final.png",
                     f"Decaying turbulence, GOW17 — t = {pairs[-1][1]['t']:.3f} Myr")
    print("wrote", outdir / "fig_compare_slices_final.png")

    ak_h, ap_h = history(ak), history(ap)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    for ax, (key, label) in zip(axes, [
        ("H2", r"mass-weighted $2x({\rm H_2})$"),
        ("CO", r"mass-weighted $x({\rm CO})/x_{\rm C,tot}$"),
        ("Pk", r"mass-weighted $P/k_{\rm B}$ [cm$^{-3}$ K]"),
    ]):
        ax.plot(ap_h["t"], ap_h[key], "--", color="#2E8B57", lw=1.8,
                label="Athena++ CVODE")
        ax.plot(ak_h["t"], ak_h[key], "-", color="#C1440E", lw=2.0,
                label="AthenaK semi-implicit")
        ax.set(xlabel="t [Myr]", ylabel=label, yscale="log")
        ax.grid(alpha=0.25, lw=0.5)
    axes[0].legend(frameon=False, fontsize=9)
    nx = ak[-1]["nH"].shape[-1]
    fig.suptitle(f"Decaying turbulence from a shared snapshot, "
                 f"{ak[-1]['t']:.0f} Myr, {nx}$^3$")
    fig.savefig(outdir / "fig_turb_history.png", dpi=140)
    print("wrote", outdir / "fig_turb_history.png")


if __name__ == "__main__":
    main()
