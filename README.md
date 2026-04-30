# CFD WENO, TENO, MUSCL-THINC

Compact Fortran/Python project for comparing high-resolution finite-volume
schemes on compressible Euler and Navier-Stokes shock problems.

The code compares:

- TENO
- WENO-JS
- MUSCL-THINC
- HLLC approximate Riemann fluxes

## Results

### 2D Riemann Problem at t = 0.8

![2D Riemann density comparison at t = 0.8](2D%20Riemann/Euler.png)

### Sod Shock Tube

![Sod shock tube density comparison](SOD/Sod%20Tube.png)

### 2D Viscous Shock Tube

![Bottom wall density comparison](2D%20Shock%20Tube/bottom_wall.png)

## Cases

| Folder | Purpose |
| --- | --- |
| `SOD/` | 1D Sod shock-tube validation |
| `SOD - Calibrate/` | Parameter sweeps for TENO and MUSCL-THINC |
| `SOD- Order/` | Grid-refinement and order study |
| `2D Riemann/` | 2D inviscid Riemann problem and animations |
| `2D Shock Tube/` | 2D viscous shock-tube post-processing |

## Run

Compile and run a Fortran case from its folder, for example:

```powershell
gfortran -O3 -fopenmp -o solver.exe ri-six.f90
.\solver.exe
```

For the 2D Riemann animation workflow:

```powershell
.\run_riemann_parallel.ps1 -Mode preview -BoundaryMode REFL
```

Python plotting scripts read the generated `.dat` or snapshot files and write
the figures shown above.

## Notes

Generated `.dat`, `.png`, `.mp4`, `frames/`, and `build/` outputs are result
artifacts. The main source files are the `.f90`, `.py`, and `.ps1` scripts in
each case folder.
