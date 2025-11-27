# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
# Initial build
mkdir build && cd build
cmake ..
make -j4

# Rebuild after changes
cd build && make -j4
```

## Running Simulations

```bash
cd build

# Cavity flow (lid-driven cavity)
./cavity_flow [nx] [U_lid] [end_time]
# Example: ./cavity_flow 128 0.01 20.0

# Channel flow (Poiseuille flow)
./channel_flow [nx] [ny] [U_in] [end_time]
# Example: ./channel_flow 256 32 0.01 5.0
```

## Visualization (from build directory)

```bash
# Final state plot (PDF output)
python3 ../scripts/visualize.py output_channel --plot-final

# Animation
python3 ../scripts/visualize.py output_channel --animation

# Convergence analysis
python3 ../scripts/convergence.py output_channel

# Validation against Ghia benchmark (cavity flow)
python3 ../scripts/validation.py output_cavity --Re 100
```

## Architecture

This is a 2D incompressible Navier-Stokes solver using **MAC method** (staggered grid) with **Projection method** for time integration.

### Core Classes (in `fluid` namespace)

- **Grid**: Staggered grid with pressure at cell centers, velocities at cell faces. Arrays include ghost cells for boundary conditions.
- **Solver**: Time integration using Projection method (predictor-corrector). Handles convection (1st-order upwind) and diffusion (2nd-order central).
- **PressureSolver**: SOR iterative solver for pressure Poisson equation.
- **BoundaryCondition**: Supports NoSlip, Inflow, Outflow. Factory methods: `cavityFlow()`, `channelFlow()`.
- **CSVWriter**: Outputs field data to CSV files in `output_*/data/` directories.

### Simulation Loop (Projection Method)

1. Compute intermediate velocity u* (without pressure)
2. Solve pressure Poisson equation: ∇²p = (ρ/Δt)∇·u*
3. Correct velocity: u^{n+1} = u* - (Δt/ρ)∇p

### Key Implementation Details

- Grid indexing: `p[i][j]` at cell center, `u[i][j]` at left face, `v[i][j]` at bottom face
- Ghost cells: p is (nx+2)×(ny+2), u is (nx+1)×(ny+2), v is (nx+2)×(ny+1)
- OpenMP parallelization available via `USE_OPENMP` compile flag (auto-detected by CMake)
- All physical quantities use SI units (meters, seconds, kg/m³, Pa)

## Python Dependencies

Scripts require: numpy, matplotlib, pandas, scipy
