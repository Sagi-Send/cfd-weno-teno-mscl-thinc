# High-Resolution Finite-Volume Solvers for Canonical CFD Tests

This repository contains a collection of Fortran 90 solvers and companion Python
post-processing scripts for validating modern shock-capturing techniques on
classical one- and two-dimensional compressible flow benchmarks. The codebase
implements MUSCL-THINC, WENO-JS, and TENO reconstructions, an HLLC approximate
Riemann solver, and third-order TVD Runge–Kutta time marching to resolve strong
shocks and contact discontinuities that arise in Sod's shock tube, 2D Riemann
problems, and related experiments.

## Repository layout

| Path | Description |
| --- | --- |
| `SOD/` | Baseline 1D Sod shock-tube solver that sweeps over WENO, TENO, and MUSCL-THINC reconstructions and writes solution snapshots for plotting.【F:SOD/sod.f90†L752-L821】 |
| `SOD - Calibrate/` | Parameter-study driver that scans TENO cut-off values (`ct`) and THINC sharpness (`beta`) before writing per-run `.dat` files for calibration plots.【F:SOD - Calibrate/sod-cal.f90†L620-L780】 |
| `SOD- Order/` | Density-wave accuracy test and grid-refinement study with WENO/TENO reconstructions and MUSCL-THINC hybridisation routines shared with the Sod solver.【F:SOD- Order/oa.f90†L200-L512】 |
| `2D Shock Tube/` | Finite-volume solver and plotting utilities for planar shock-tube initial data extended to two dimensions (`st.f90`). |
| `2D Riemann/` | Four-quadrant 2D Riemann-problem setup (`ri-six.f90`) mirroring the 1D architecture for multi-directional discontinuities. |
| `Literature/` | PDFs of the MUSCL-THINC, TENO, and AUSM-D reference papers used to implement the numerical schemes. |
| `*.py` | Matplotlib post-processing utilities for comparing numerical density/temperature traces against exact data and for showing convergence trends.【F:SOD/sod-plot.py†L1-L80】【F:SOD- Order/oa-plot.py†L1-L80】 |

Compiled `.exe` binaries and `.mod` module artifacts are checked in so you can
inspect reference results without rebuilding, but rebuilding from source is
recommended whenever you modify any solver.

## Requirements

* **Fortran toolchain:** Any modern `gfortran` (tested with 11+) or Intel
  compiler capable of Fortran 90 modules.
* **Python 3.9+** with `numpy` and `matplotlib` for the plotting utilities. No
  extra dependencies are required because the scripts rely solely on the Python
  standard library plus those two packages.

Optional: `gnuplot` or another plotting package if you prefer to use the plain
text `.dat` files elsewhere.

## Building and running the solvers

Every solver follows the same pattern: edit the problem parameters near the
bottom of the corresponding `.f90` file, compile with your preferred Fortran
compiler, run the generated executable, and inspect the `.dat` outputs and PNGs
produced by the Python scripts.

### Example: Sod shock tube (1D)

```bash
cd SOD
# Build
gfortran -O3 -std=f2008 -o sod.exe sod.f90
# Run (produces WENO/TENO/MUSCL outputs under names such as TENO_HLLC_grid_806.dat)
./sod.exe
# Plot density traces against the exact solution
python sod-plot.py
```

Key tunables:

* `grid_size` array selects the physical resolution(s) (number of cells in the
  interior domain).
* `scheme` and `solver` arrays determine which reconstruction (TENO, WENO, or
  MUSCL-THINC) and which Riemann solver (currently HLLC) are used per sweep.
* `ct` and `sharpness` configure the TENO cut-off and THINC hyperbolic-tangent
  steepness respectively.【F:SOD/sod.f90†L752-L821】

The solver writes one `.dat` file per scheme/solver pair along with optional
PNG figures via `sod-plot.py`.

### Example: 2D Riemann problem

```bash
cd '2D Riemann'
gfortran -O3 -std=f2008 -o riemann.exe ri-six.f90
./riemann.exe
python ri-plot.py
```

The `init_riemann` routine sets up the four states at the quadrant interfaces.
You can change velocities or pressures in that subroutine to explore different
contact configurations.【F:2D Riemann/ri-six.f90†L60-L118】

### Example: 2D Shock tube

```bash
cd '2D Shock Tube'
gfortran -O3 -std=f2008 -o shock.exe st.f90
./shock.exe
python st-plot.py  # or wall-plot-1.py for boundary layers
```

The `init_shock_tube` subroutine mirrors the 1D setup but fills the entire 2D
mesh, while the shared modules in the same file handle flux construction and
Runge–Kutta stepping.【F:2D Shock Tube/st.f90†L45-L120】

### Calibration and order-study utilities

* `SOD - Calibrate` sweeps through TENO/MUSCL parameters and writes files such
  as `TENO_HLLC_ct_1e-5.dat` or `MUSC_HLLC_beta_1.5.dat` that the included
  `cal-plot.py` script turns into comparison figures for reports.
* `SOD- Order` contains a density-wave driver (`oa.f90`) plus `oa-plot.py` for
  plotting the `L1` error vs. grid spacing alongside a second-order reference
  line. Update `file_paths` in the plotter to toggle between TENO and WENO
  datasets.【F:SOD- Order/oa-plot.py†L1-L60】

## Output files

Each solver writes ASCII `.dat` files with columns of `x`, `rho`, `u`, and `p`
(or the 2D analogue). These files double as checkpoints for post-processing,
automated regression tests, and archival data for publications. Python plotting
scripts expect these filenames by default, so keep the naming convention if you
want to reuse the provided figures.

PNG figures such as `Sod Tube.png`, `Temperature.png`, and `Calibration-*.png`
are stored alongside the generating scripts for convenience.

## Reference material

The `Literature/` folder contains PDFs of the papers the implementation follows:

* **TENO:** template for adaptive stencil selection and cut-off parameter usage.
* **MUSCL-THINC:** blending of monotone upstream-centered schemes with THINC
  limiters for better interface capturing.
* **AUSMD:** source material for pressure-based flux splitting.

Consult those papers when you need deeper theoretical background or want to
cross-check coefficients.

## Extending the code

* Add additional Riemann solvers by implementing a new subroutine next to
  `riemann_HLLC` and extending `riemann_solver`'s selector block.
* Modify `activateScheme` to introduce other reconstructions (e.g., WENO-Z or
  MP5) while keeping the same finite-volume infrastructure.
* Adjust `ApplyBoundaryConditions` for inflow/outflow variations—the 1D solver
  currently supports periodic (order test) and outflow (calibration) options in
  separate directories, so you can copy whichever variant suits your study.【F:SOD- Order/oa.f90†L512-L620】【F:SOD - Calibrate/sod-cal.f90†L620-L666】

Contributions and experiment reports are welcome—open an issue describing the
scenario and attach the `.dat` files or figures that demonstrate the behaviour
under investigation.
