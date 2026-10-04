# Quarter-Sector 2D Finite Element Analysis

This is a MATLAB App Designer for static finite element analysis project completed collaboratively with my teammates from Group 9 of the Finite Element Methods course. The project performs a static finite element analysis on a pressurized annular sector. The solver supports 3-node triangular (`T3`) and 4-node quadrilateral (`Q4`) elements and outputs results for displacement and stress.

## Features

- Generate a mesh for a quarter annulus in the first quadrant.
- Choose `T3` or `Q4` elements.
- Assemble the global stiffness matrix and equivalent nodal pressure loads.
- Solve the static plane-strain problem with roller constraints on the two straight symmetry boundaries.
- Display mesh and force plots, numeric result summaries, and selectable displacement/stress contours.
- Provide an experimental modal-analysis interface.

## Requirements

- MATLAB with App Designer UI components, including `uifigure`, `uiaxes`, and `uidropdown`.
- MATLAB's `eigs` function for modal analysis.
- No third-party MATLAB toolboxes are intentionally required by the static solver.
- Keep the project `.m` files in the same folder or add their folder to the MATLAB path.

## Run the App

From MATLAB, set the Current Folder to this project directory and run:

```matlab
app = FEM_Interface;
```

The app opens with default inputs. The main defaults are:

| Input | Default | Meaning |
| --- | ---: | --- |
| Young's modulus, `E` | `2e11` | Elastic modulus |
| Poisson's ratio, `nu` | `0.3` | Poisson ratio |
| Thickness | `1` | Out-of-plane thickness |
| Inner radius, `r` | `0.5` | Inner radius |
| Outer radius, `R` | `1` | Outer radius |
| Inner pressure | `2e7` | Pressure applied at the inner boundary |
| Outer pressure | `-1e7` | Pressure applied at the outer boundary; sign is used as entered |
| Radial nodes | `11` | Number of nodes in the radial direction |
| Angular nodes | `21` | Number of nodes along the quarter arc |
| Element type | `T3` | Element formulation (`T3` or `Q4`) |

Use a consistent unit system. The defaults are consistent with SI units: lengths in metres, elastic modulus and pressure in pascals, and thickness in metres.

## Typical Workflow

1. Enter material properties, geometry, pressures, and mesh density.
2. Select `T3` or `Q4` in the element-type switch.
3. Select **MESHING** to display the generated mesh, or **APPLIED FORCE** to display the equivalent nodal force vectors.
4. Select **RUN ANALYSIS** to solve the static system and show displacement and stress summaries.
5. Choose a result in the first `TYPE` dropdown to display its contour on the axes.
6. Select **RESET** to clear the plot and result text.

The contour-selection callback currently rebuilds and solves the model when a result is selected; it does not reuse the static solution from **RUN ANALYSIS**.

## Model and Numerical Formulation

- The geometry spans radii `r` to `R` and angles from 0 to 90 degrees.
- Each node has two translational degrees of freedom, `u_x` and `u_y`.
- The straight edges at `y = 0` and `x = 0` are constrained as roller/symmetry boundaries in the static solve.
- The constitutive matrix in the app uses isotropic linear elasticity under plane strain.
- `Q4` elements use Gaussian quadrature; `T3` elements use a constant-strain triangle formulation.
- The stress postprocessor computes Cartesian, radial, circumferential, out-of-plane, shear, and von Mises stress values.

## Project Files

| File | Purpose |
| --- | --- |
| `FEM_Interface.m` | App Designer interface and analysis, plotting, and modal callbacks |
| `mesh2D.m` | Generates node coordinates and `T3`/`Q4` connectivity for the quarter annulus |
| `stiffnessmat2D.m` | Assembles the static global stiffness matrix |
| `stiffnessmat2D_modal.m` | Assembles global stiffness and mass matrices for modal analysis |
| `force2D.m` | Converts inner and outer boundary pressures into an equivalent nodal load vector |
| `stress_cal.m` | Calculates element stress components and von Mises stress |
| `quadrature2D.m` | Supplies integration points and weights |
| `Larange2D.m` | Evaluates element shape functions and derivatives |
| `plot_meshing.m` | Standalone mesh-generation and plotting example |

## Current Limitations

- The modal callback in `FEM_Interface.m` requests both stiffness and mass matrices from `stiffnessmat2D.m`, although that function returns only the stiffness matrix. `stiffnessmat2D_modal.m` is the routine that returns both matrices, but it is not currently called by the callback.
- Several stress-contour branches in the result dropdown refer to variable names that do not match the outputs returned by `stress_cal.m`. Some stress selections may therefore fail even though the static analysis summary can be produced.
- The modal interface and some stress-contour selections should be verified or corrected before relying on them for results.
