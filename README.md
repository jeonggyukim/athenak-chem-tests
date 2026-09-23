# AthenaK chemistry tests

Inputs, run scripts and figure scripts that compare the AthenaK chemistry
module (GOW17 network, semi-implicit and Kokkos BDF solvers) with Athena++ and
CVODE. Simulation output is not committed; the scripts rebuild it.

| directory | holds |
| --- | --- |
| `reproduce/inputs/` | the four input files |
| `reproduce/run_turb64.sh` | runs every stage, from the driven box to the figures |
| `reproduce/patches/` | the one-file Athena++ fix the comparison needs |
| `reproduce/log.md` | what was run, on which commits, and what it measured |
| `scripts/` | plotting and comparison scripts |

## The test: 64^3 decaying turbulence

1. AthenaK drives hydro-only turbulence in a periodic 32 pc box at 64^3,
   n_H = 10 cm^-3, from T = 10^4 K to t = 12 code units (`drive64.athinput`).
2. The last dump is converted to a legacy VTK file, `shared_ic64.vtk`.
3. AthenaK and Athena++ both read that file and evolve it for 2.0454 code
   units (2 Myr) with GOW17 chemistry and no forcing (`turb64_si.athinput`,
   `turb64_bdf.athinput`, `turb64_cvode.athinput`).

The two codes have different turbulence drivers, so sharing the file is what
makes a cell-by-cell comparison possible. Both use RK2, PLM and HLLC, split the
chemistry from the hydro once per step, and use the constant H2 grain
formation rate of 3e-17 cm^3 s^-1 (`is_kgrH2_const = true` in AthenaK).

## Getting the codes

**AthenaK**, with the `GOW17_read_vtk` problem generator:

    git clone --recursive -b semi-implicit https://github.com/jeonggyukim/athenak.git ~/athenak
    cd ~/athenak && git checkout 2deb5e40
    cmake -B build -DCMAKE_BUILD_TYPE=Release -DAthena_ENABLE_MPI=ON
    make -C build -j8

**Athena++** at upstream `823614c`, plus a compile fix for `read_vtk.cpp`
under MPI (upstream assigns to an undeclared `ierr`). CVODE comes from
SUNDIALS; 6.7.0 was used.

    git clone https://github.com/PrincetonUniversity/athena.git ~/athena
    cd ~/athena && git checkout 823614c
    git apply <this repo>/reproduce/patches/athena-pp-read_vtk-mpi.diff
    python configure.py --prob=read_vtk --chemistry=gow17 --chem_radiation=const \
        --chem_ode_solver=cvode --cvode_path=<sundials prefix> -mpi
    make -j8
    gcc -Wall -O2 -o vis/vtk/join_vtk++ vis/vtk/join_vtk++.c -lm

**Python** 3 with numpy and matplotlib, and `ffmpeg` on `PATH` for the movie.

## Running

From the repository root:

    bash reproduce/run_turb64.sh                   # every stage
    bash reproduce/run_turb64.sh drive ic          # only the shared initial condition
    bash reproduce/run_turb64.sh si bdf cvode join plot
    SMOKE=1 bash reproduce/run_turb64.sh           # every stage for 20 cycles, ~2 min

| stage | does | time, 8 cores |
| --- | --- | --- |
| `drive` | AthenaK driven turbulence, hydro only, 1 rank | 4.6 min |
| `ic` | last `drive64` dump -> `shared_ic64.vtk` | seconds |
| `si` | AthenaK semi-implicit | 28 s |
| `bdf` | AthenaK Kokkos BDF, rtol 1e-2, atol 1e-15 | ~35-90 min, estimated |
| `cvode` | Athena++ CVODE, rtol 1e-2, atol 1e-15 | 17.5 min |
| `join` | merge Athena++'s per-MeshBlock VTK files | seconds |
| `plot` | history from `.hst` files and dumps, slice comparison and movie | ~1 min |

Times are for an Apple-silicon MacBook Pro.

The script expects the checkouts at `~/athenak` and `~/athena` and writes to
`data/` in this repository (git ignores it). Override any path from the
environment: `ATHENAK_ROOT`, `ATHENAPP_ROOT`, `AKBIN`, `APBIN`, `CHEMDATA`,
`PYATHENA`, `NP`. See `reproduce/machines/default.sh`; the other files there
are the author's machines and are chosen only on those hosts.

Every run writes `versions.txt`, the diff of each code's working tree and its
untracked sources beside its output, so the code behind a result stays
recoverable.

Rerunning `drive` and `ic` with a different compiler and build (g++-15 serial
against g++-16 with MPI, same machine) reproduced `shared_ic64.vtk` byte for
byte. Another CPU architecture may differ in the last bits; that does not
affect the comparison, since both codes read whichever file was made.
