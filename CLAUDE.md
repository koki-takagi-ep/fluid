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

# Cavity flow (lid-driven cavity) - Projection method
./cavity_flow [nx] [U_lid] [end_time]
# Example: ./cavity_flow 128 0.01 20.0

# Channel flow (Poiseuille flow) - Projection method
./channel_flow [nx] [ny] [U_in] [end_time]
# Example: ./channel_flow 256 32 0.01 5.0

# Channel flow - SIMPLE method
./channel_flow_simple [nx] [ny] [U_in] [end_time]
# Example: ./channel_flow_simple 64 16 0.01 2.0
```

## Visualization (from build directory)

```bash
python3 ../scripts/visualize.py output/channel_projection --plot-final   # Projection法
python3 ../scripts/visualize.py output/channel_simple --plot-final       # SIMPLE法
python3 ../scripts/visualize.py output/cavity_projection --plot-final    # Cavity flow
python3 ../scripts/validation.py output/cavity_projection --Re 100       # Ghia benchmark
```

## Architecture

2D incompressible Navier-Stokes solver using **MAC method** (staggered grid) with two time integration schemes:
- **Projection method** (Chorin, 1968) - `Solver` class
- **SIMPLE method** (Patankar & Spalding, 1972) - `SimpleSolver` class

### Core Classes (in `fluid` namespace)

| Class | File | Description |
|-------|------|-------------|
| `Grid` | Grid.hpp/cpp | Staggered grid: p at cell centers, u/v at faces. Ghost cells for BCs. |
| `Solver` | Solver.hpp/cpp | Projection method time integration |
| `SimpleSolver` | SimpleSolver.hpp/cpp | SIMPLE method time integration |
| `PressureSolver` | PressureSolver.hpp/cpp | SOR solver for pressure Poisson equation |
| `BoundaryCondition` | BoundaryCondition.hpp/cpp | NoSlip, Inflow, Outflow. Factory: `cavityFlow()`, `channelFlow()` |
| `CSVWriter` | CSVWriter.hpp/cpp | CSV output to `output_*/data/` |

### Constants

Numerical constants are defined in `include/Constants.hpp` to avoid magic numbers:
- `ATMOSPHERIC_PRESSURE`, `WATER_DENSITY`, `WATER_KINEMATIC_VISCOSITY`
- `DEFAULT_SOR_OMEGA`, `DEFAULT_PRESSURE_TOLERANCE`
- `LAPLACIAN_CENTER_COEFF`, `FOUR_POINT_AVERAGE_COEFF`

### Simulation Loop

**Projection Method:**
1. Compute intermediate velocity u* (without pressure gradient)
2. Solve pressure Poisson: ∇²p = (ρ/Δt)∇·u*
3. Correct velocity: u^{n+1} = u* - (Δt/ρ)∇p

**SIMPLE Method:**
Same structure but with under-relaxation for pressure (α_p ≈ 0.3) and velocity (α_u ≈ 0.7).

### Grid Indexing

- `p[i][j]`: cell center, array size (nx+2)×(ny+2)
- `u[i][j]`: left face of cell (i,j), array size (nx+1)×(ny+2)
- `v[i][j]`: bottom face of cell (i,j), array size (nx+2)×(ny+1)
- Internal cells: `1 <= i <= nx`, `1 <= j <= ny`
- Ghost cells: index 0 and nx+1 (or ny+1)

### Spatial Discretization

- Convection: 1st-order upwind
- Diffusion: 2nd-order central difference
- Pressure: Red-Black SOR

### Physical Units

All quantities in SI units: meters, seconds, kg/m³, Pa.

## Output Structure

```
output/
├── cavity_projection/      # Cavity flow (Projection method)
├── channel_projection/     # Channel flow (Projection method)
└── channel_simple/         # Channel flow (SIMPLE method)
    ├── data/
    │   ├── field_000000.csv    # x, y, u, v, p, magnitude
    │   ├── field_000001.csv
    │   └── metadata.csv        # Grid parameters
    ├── figures/
    │   └── result.svg          # Visualization output (SVG format)
    └── simulation.log          # Computation time and settings log
```

## Python Dependencies

numpy, matplotlib, pandas, scipy

## Visualization Scripts

| Script | Purpose |
|--------|---------|
| `visualize.py` | Plot velocity/pressure fields, streamlines (SVG output) |
| `validation.py` | Compare with Ghia et al. (1982) benchmark data |
| `convergence.py` | Analyze solver convergence history |
