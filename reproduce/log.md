# Run log

Every figure traces to a row here: the code it ran against, the command, and
what the run reported. Raw output is not committed; `run_turb64.sh` rebuilds it
and writes `versions.txt`, the diffs and the untracked sources beside it.

Add a row when a run happens, not when a figure is drawn.

## 2026-09-23

Code: athenak `6f8051cf` on `semi-implicit`, plus the uncommitted
`GOW17_read_vtk` pgen. athena-pp `823614c` on `main`, plus a modified
`src/pgen/read_vtk.cpp`. Machine: MacBook Pro (Apple silicon), OpenMPI.

### Decaying turbulence, first run (superseded)

Output: `~/Documents/Athenak-chem-tests/decay-cvode-compare/turb64_{si,cvode}_run`.
Run by an earlier `run_turb64.sh` with the same stages, before this repository
existed.

The hydro schemes differed: Athena++ ran its default `vl2` integrator and
AthenaK `rk2`. Athena++ wrote no history file. Both are fixed in
`inputs/turb64_cvode.athinput`, so these numbers are not reproducible from the
current inputs.

| stage | ranks | wall |
| --- | --- | --- |
| drive (serial `build-cpu`) | 1 | 4 min 39 s |
| si | 8 | 35.1 s |
| cvode | 8 | 17 min 33 s; 9.37e4 zone-cycles per CPU-second |

`shared_ic64.vtk` is `drive64.hydro_w.00006.bin` (t = 12 code units),
converted by `athenak_bin_to_vtk.py`; reconverting that dump reproduces the
file byte for byte.

### Smoke test, all stages at 20 cycles

    SMOKE=1 CHEMDATA=<scratch>/smoke bash reproduce/run_turb64.sh

Passed, every stage exit 0. The cvode stage took 48.3 s wall on 8 ranks. The
drive ran for only 20 cycles, so the shared field is nearly uniform
(n_H = 9.76 to 10.36 cm^-3), and the numbers below test the pipeline, not the
solver. At t = 0.306 Myr, AthenaK/Athena++ - 1: mass +1.5e-11,
kinetic energy +8.7e-4, thermal energy +8.5 (AthenaK 9.5 times higher),
H2 -0.94, CO +1.4.

### Decaying turbulence, AthenaK reruns

Output: `~/Documents/Athenak-chem-tests/turb64/`. The Athena++ reference is
still the first run above (VL2, no history file). athenak built with
`semi_implicit_adaptive = true` as the default (uncommitted in athenak).

| run | directory | wall (8 ranks) | note |
| --- | --- | --- | --- |
| si, T-dependent H2 grain rate | `turb64_si_kgrT_run` | 29.5 s | AthenaK default `is_kgrH2_const = false` |
| si, constant H2 grain rate | `turb64_si_run` | 28.4 s | `bash reproduce/run_turb64.sh si` |
| si, constant rate, tabulated | `turb64_si_table_run` | 19.6 s | `chemistry/GOW17_thermo_table=true chemistry/semi_implicit_table_deriv=true` on the command line |
| Kokkos BDF, constant rate | `turb64_bdf_run` | stopped at t = 0.196 Myr after ~9 min | rtol 1e-2, atol 1e-15 |

With the T-dependent rate, H2 at t = 0.021 Myr is 0.02 to 0.34 of the
Athena++ value, falling with initial temperature as the ratio of the two rate
coefficients does, and the thermal energy falls more slowly than in Athena++
over the first 0.5 Myr. With the constant rate, AthenaK/Athena++ - 1 at 2 Myr
is -0.9% in thermal energy, -0.9% in kinetic energy, +1.6% in H2 and +26% in
x_CO/x_C,tot. Tabulation changes these by less than 0.5 percentage points.
Kokkos BDF matches Athena++ in CO to within 3% up to 0.196 Myr, where the
semi-implicit runs are 11-22% high.

Zone-cycles per second per core:

| run | chemistry | hydro + MPI + output | whole run |
| --- | --- | --- | --- |
| Athena++ CVODE (VL2) | 1.27e4 | 1.51e5 | 1.17e4 |
| AthenaK semi-implicit | 5.45e5 | 5.18e5 | 2.65e5 |
| AthenaK semi-implicit, tabulated | 1.23e6 | 5.55e5 | 3.82e5 |

Chemistry time is the sum of the per-step `chemistry kernel` lines (AthenaK)
or the per-MeshBlock `chemistry ODE integration` lines (Athena++), divided by
the number of ranks.
