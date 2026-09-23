# AthenaK chemistry figures

Figures for the AthenaK chemistry module (GOW17 network, semi-implicit solver),
the scripts that make them, and the inputs and run scripts that produce their
data.

This repository is meant to be shareable, so nothing private and nothing bulky
belongs here. Raw simulation output lives in `$CHEMDATA`, by default
`~/Documents/Athenak-chem-tests/turb64`.

| directory | holds |
| --- | --- |
| `reproduce/` | input files, run scripts, the VTK converter, and `log.md` |
| `scripts/` | plotting and comparison scripts |
| `figures/` | the output; tracked only once a figure is near final |

## Decaying 64^3 turbulence: semi-implicit against CVODE

    bash reproduce/run_turb64.sh                  # drive, ic, si, cvode, join, plot
    bash reproduce/run_turb64.sh si cvode join    # rerun the chemistry only
    SMOKE=1 bash reproduce/run_turb64.sh          # every stage at 20 cycles, ~2 min

AthenaK drives hydro-only turbulence to t = 12. Its last dump becomes
`shared_ic64.vtk`, and both codes start from that file and evolve it without
forcing for 2.0454 code units (2 Myr) with GOW17 chemistry. The hydro scheme
matches in both: rk2 time integrator, PLM reconstruction, HLLC Riemann solver.

The plot stage writes three comparisons: `fig_turb_hst.png` (kinetic and
thermal energy, dt and mass from both `.hst` files), `fig_turb_history_dumps.png`
(volume integrals, H2 and CO from the dumps), and the slice comparison
`fig_compare_slices_final.png` with its movie `compare_slices.mp4`.

Machine settings live in `reproduce/machines/{mac,grammar,syntax}.sh`, chosen
by hostname or `MACHINE=`. Every path there can be overridden from the
environment. The Slurm command lines are in the header of `run_turb64.sh`.

| machine | runs | status |
| --- | --- | --- |
| MacBook | every stage, AthenaK on CPU | smoke test passed 2026-09-23 |
| grammar | every stage, AthenaK on CPU | not yet run; needs both codes built there |
| syntax | drive, ic, si with AthenaK on one GPU | not yet run; cvode goes to grammar via the shared `/gpfs` |

Requirements:

| what | where | note |
| --- | --- | --- |
| AthenaK | `$AKBIN` | needs the `GOW17_turb` and `GOW17_read_vtk` pgens |
| Athena++ | `$APBIN` | configured with `--prob=read_vtk --flux=hllc`, GOW17 chemistry, CVODE |
| `join_vtk++` | `$ATHENAPP_ROOT/vis/vtk/` | compiled from `join_vtk++.cpp` |
| python | `$PYATHENA` | numpy, matplotlib, ffmpeg on `PATH`; reads `$ATHENAK_ROOT/vis/python/bin_convert.py` |

Each run writes `versions.txt`, `<code>.diff` and `<code>-untracked.tar` beside
its output, recording the commits and uncommitted changes of both codes.
Add a row to `reproduce/log.md` after every run.
