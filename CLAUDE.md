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

# Channel flow - PISO method
./channel_flow_piso [nx] [ny] [U_in] [end_time] [nCorrectors]
# Example: ./channel_flow_piso 64 16 0.01 2.0 2

# Cavity flow - PISO method
./cavity_flow_piso [nx] [U_lid] [end_time] [nCorrectors]
# Example: ./cavity_flow_piso 64 0.01 10.0 2

# Grid convergence study (multiple solvers/limiters)
./channel_flow_convergence <solver> <limiter> <nx> <ny> <U_max> <end_time> <output_dir>
# solver: projection, simple, piso
# limiter: none, minmod, superbee, vanleer, mc
# Example: ./channel_flow_convergence projection vanleer 128 32 0.015 5.0 output/conv_test

# TVD scheme benchmark (compares all limiters)
./tvd_benchmark [nx] [U_lid] [end_time]
# Example: ./tvd_benchmark 64 0.01 10.0
```

## Visualization (from build directory)

```bash
python3 ../scripts/visualize.py output/channel_projection --plot-final   # Projection法
python3 ../scripts/visualize.py output/channel_simple --plot-final       # SIMPLE法
python3 ../scripts/visualize.py output/channel_piso --plot-final         # PISO法
python3 ../scripts/visualize.py output/cavity_projection --plot-final    # Cavity flow
python3 ../scripts/visualize.py output/cavity_piso --plot-final          # Cavity PISO
python3 ../scripts/validation.py output/cavity_projection --Re 100       # Ghia benchmark
```

## Architecture

2D incompressible Navier-Stokes solver using **MAC method** (staggered grid) with three time integration schemes:
- **Projection method** (Chorin, 1968) - `Solver` class
- **SIMPLE method** (Patankar & Spalding, 1972) - `SimpleSolver` class
- **PISO method** (Issa, 1986) - `PisoSolver` class

### Class Hierarchy

```
SolverBase (abstract)
├── Solver        (Projection method)
├── SimpleSolver  (SIMPLE method)
└── PisoSolver    (PISO method)
```

### Core Classes (in `fluid` namespace)

| Class | File | Description |
|-------|------|-------------|
| `Grid` | Grid.hpp/cpp | Staggered grid: p at cell centers, u/v at faces. Ghost cells for BCs. |
| `SolverBase` | SolverBase.hpp/cpp | Abstract base: time step, convection/diffusion, velocity correction |
| `Solver` | Solver.hpp/cpp | Projection method (inherits SolverBase) |
| `SimpleSolver` | SimpleSolver.hpp/cpp | SIMPLE method (inherits SolverBase) |
| `PisoSolver` | PisoSolver.hpp/cpp | PISO method (inherits SolverBase) |
| `PressureSolver` | PressureSolver.hpp/cpp | SOR solver for pressure Poisson equation |
| `BoundaryCondition` | BoundaryCondition.hpp/cpp | NoSlip, Inflow, Outflow. Factory: `cavityFlow()`, `channelFlow()` |
| `CSVWriter` | CSVWriter.hpp/cpp | CSV output to `output_*/data/` |

### Constants

Numerical constants are defined in `include/Constants.hpp` to avoid magic numbers:
- Physical: `WATER_DENSITY`, `WATER_DYNAMIC_VISCOSITY`, `WATER_KINEMATIC_VISCOSITY`
- Solver: `DEFAULT_CFL`, `DEFAULT_SOR_OMEGA`, `DEFAULT_PRESSURE_TOLERANCE`
- SIMPLE: `DEFAULT_VELOCITY_RELAXATION`, `DEFAULT_PRESSURE_RELAXATION`
- Simulation: `DEFAULT_CAVITY_END_TIME`, `DEFAULT_CHANNEL_END_TIME`
- Discretization: `LAPLACIAN_CENTER_COEFF`, `POISEUILLE_MEAN_MAX_RATIO`

### Simulation Loop

**Projection Method:**
1. Compute intermediate velocity u* (without pressure gradient)
2. Solve pressure Poisson: ∇²p = (ρ/Δt)∇·u*
3. Correct velocity: u^{n+1} = u* - (Δt/ρ)∇p

**SIMPLE Method:**
Same structure but with under-relaxation for pressure (α_p ≈ 0.3) and velocity (α_u ≈ 0.7).

**PISO Method:**
Non-iterative predictor-corrector scheme with multiple correction steps:
1. Predictor: Compute u* (same as projection)
2. First Corrector: Solve ∇²p' = (ρ/Δt)∇·u*, correct u** = u* - (Δt/ρ)∇p'
3. Second Corrector: Solve ∇²p'' = (ρ/Δt)∇·u**, correct u^{n+1} = u** - (Δt/ρ)∇p''

PISO is particularly suited for transient flows (no outer iteration required).

**Important:** All solvers retain the previous pressure field (`grid.p`) as the initial guess for SOR iteration. This enables 2nd-order grid convergence for all three methods.

### Grid Convergence (Verified)

Channel flow (Poiseuille) L2 error convergence orders:
- Projection: p = 1.95
- SIMPLE: p = 1.95
- PISO: p = 1.98

### Grid Indexing

- `p[i][j]`: cell center, array size (nx+2)×(ny+2)
- `u[i][j]`: left face of cell (i,j), array size (nx+1)×(ny+2)
- `v[i][j]`: bottom face of cell (i,j), array size (nx+2)×(ny+1)
- Internal cells: `1 <= i <= nx`, `1 <= j <= ny`
- Ghost cells: index 0 and nx+1 (or ny+1)

### Spatial Discretization

- Convection: 1st-order upwind (default) or TVD schemes
- Diffusion: 2nd-order central difference
- Pressure: Red-Black SOR

### TVD Schemes (FluxLimiter.hpp)

Optional TVD flux limiters for convection terms:
- `LimiterType::None` - 1st-order upwind (default)
- `LimiterType::Minmod` - Most diffusive, very stable
- `LimiterType::Superbee` - Least diffusive, can be oscillatory
- `LimiterType::VanLeer` - Good balance (recommended)
- `LimiterType::MC` - Monotonized Central

Usage: `solver->setLimiter(fluid::LimiterType::VanLeer);`

### Physical Units

All quantities in SI units: meters, seconds, kg/m³, Pa.

## Output Structure

```
output/
├── cavity_projection/      # Cavity flow (Projection method)
├── cavity_piso/            # Cavity flow (PISO method)
├── channel_projection/     # Channel flow (Projection method)
├── channel_simple/         # Channel flow (SIMPLE method)
└── channel_piso/           # Channel flow (PISO method)
    ├── data/
    │   ├── field_000000.csv    # x, y, u, v, p, magnitude
    │   ├── field_000001.csv
    │   └── metadata.csv        # Grid parameters
    ├── figures/
    │   ├── result.svg          # Visualization output (SVG format)
    │   ├── result.pdf          # PDF format
    │   └── result.png          # PNG format (300 dpi)
    └── simulation.log          # Computation time and settings log
```

## Python Dependencies

numpy, matplotlib, pandas, scipy

**Note:** Use Python 3.11 with Homebrew packages:
```bash
/opt/homebrew/opt/python@3.11/bin/python3.11 scripts/visualize.py ...
```

## Visualization Scripts

| Script | Purpose |
|--------|---------|
| `visualize.py` | Plot velocity/pressure fields, streamlines (SVG/PDF/PNG output) |
| `validation.py` | Compare with Ghia et al. (1982) benchmark data |
| `validation_comparison.py` | Compare all 3 solvers against Ghia benchmark |
| `convergence.py` | Analyze solver convergence history |
| `convergence_all_schemes.py` | Grid convergence study for multiple solver/limiter combinations |
| `tvd_validation.py` | Compare TVD limiters against Ghia benchmark |
| `colors.py` | Wong's colorblind-friendly color palette for plots |

## Git/PR ルール

- **ブランチ名**: プレフィックス付きのシンプルな名前を使用
  - `feature/` - 新機能（例: `feature/tvd-scheme`）
  - `fix/` - バグ修正（例: `fix/convergence-plot`）
  - `docs/` - ドキュメント（例: `docs/readme-update`）
  - `refactor/` - リファクタリング（例: `refactor/solver-cleanup`）
  - `test/` - テスト追加（例: `test/boundary-condition`）
- **PRタイトル・本文**: 日本語で記述
- **コミットメッセージ**: 英語（Conventional Commits形式）
