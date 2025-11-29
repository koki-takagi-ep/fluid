#!/usr/bin/env python3
"""
Grid convergence study for all solver/limiter combinations

This script runs channel flow simulations with multiple solver types and
TVD limiters, then plots the L2 error convergence for comparison.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed


def poiseuille_velocity(y: np.ndarray, H: float, U_max: float) -> np.ndarray:
    """Analytical Hagen-Poiseuille velocity profile"""
    eta = 2.0 * y / H - 1.0
    return U_max * (1.0 - eta**2)


def load_simulation_data(output_dir: str):
    """Load simulation results"""
    data_dir = os.path.join(output_dir, "data")
    if not os.path.exists(data_dir):
        data_dir = output_dir

    metadata_file = os.path.join(data_dir, "metadata.csv")
    df_meta = pd.read_csv(metadata_file)
    meta = dict(zip(df_meta['parameter'], df_meta['value']))
    nx, ny = int(meta['nx']), int(meta['ny'])
    lx, ly = float(meta['lx']), float(meta['ly'])
    dx, dy = lx / nx, ly / ny

    files = sorted(glob.glob(os.path.join(data_dir, "field_*.csv")))
    if not files:
        raise FileNotFoundError("No field files found")

    last_file = files[-1]
    df = pd.read_csv(last_file, comment='#')

    x = df['x'].values.reshape(nx, ny)
    y = df['y'].values.reshape(nx, ny)
    u = df['u'].values.reshape(nx, ny)

    return {
        'nx': nx, 'ny': ny, 'lx': lx, 'ly': ly,
        'dx': dx, 'dy': dy,
        'x': x, 'y': y, 'u': u
    }


def compute_l2_error(data: dict, U_max: float, x_position: float = 0.9) -> dict:
    """Compute L2 error norm against analytical solution"""
    nx, ny = data['nx'], data['ny']
    lx, ly = data['lx'], data['ly']
    h = data['dy']

    i = int(round(x_position * (nx - 1)))
    i = min(max(i, 0), nx - 1)

    y_sim = data['y'][i, :]
    u_sim = data['u'][i, :]

    u_exact = poiseuille_velocity(y_sim, ly, U_max)
    error = u_sim - u_exact

    l2_error = np.sqrt(np.mean(error**2))
    linf_error = np.max(np.abs(error))

    return {
        'h': h,
        'ny': ny,
        'l2_error': l2_error,
        'linf_error': linf_error
    }


def run_single_simulation(args):
    """Run a single simulation (for parallel execution)"""
    solver, limiter, ny, U_max, end_time, build_dir = args
    nx = ny * 4
    output_dir = f"output/conv_{solver}_{limiter}_ny{ny}"
    full_output_dir = os.path.join(build_dir, output_dir)

    # Remove existing output
    subprocess.run(f"rm -rf {full_output_dir}", shell=True)

    # Run simulation
    cmd = f"./channel_flow_convergence {solver} {limiter} {nx} {ny} {U_max} {end_time} {output_dir}"
    result = subprocess.run(
        cmd, shell=True, cwd=build_dir,
        capture_output=True, text=True
    )

    if result.returncode != 0:
        return None, f"Error: {result.stderr}"

    # Compute error
    try:
        data = load_simulation_data(full_output_dir)
        errors = compute_l2_error(data, U_max)
        errors['solver'] = solver
        errors['limiter'] = limiter
        errors['nx'] = nx
        return errors, None
    except Exception as e:
        return None, str(e)


def run_convergence_study(schemes: list, ny_list: list, U_max: float = 0.015,
                          end_time: float = 5.0, build_dir: str = '.', parallel: bool = True):
    """Run convergence study for multiple schemes"""
    results = []

    # Prepare all simulation parameters
    all_args = []
    for solver, limiter in schemes:
        for ny in ny_list:
            all_args.append((solver, limiter, ny, U_max, end_time, build_dir))

    total = len(all_args)
    print(f"\nRunning {total} simulations...")

    if parallel:
        # Parallel execution
        with ProcessPoolExecutor(max_workers=4) as executor:
            futures = {executor.submit(run_single_simulation, args): args for args in all_args}
            for i, future in enumerate(as_completed(futures)):
                args = futures[future]
                solver, limiter, ny, _, _, _ = args
                errors, err_msg = future.result()
                if errors:
                    results.append(errors)
                    print(f"[{i+1}/{total}] {solver}+{limiter} ny={ny}: L2={errors['l2_error']:.2e}")
                else:
                    print(f"[{i+1}/{total}] {solver}+{limiter} ny={ny}: FAILED - {err_msg}")
    else:
        # Sequential execution
        for i, args in enumerate(all_args):
            solver, limiter, ny, _, _, _ = args
            errors, err_msg = run_single_simulation(args)
            if errors:
                results.append(errors)
                print(f"[{i+1}/{total}] {solver}+{limiter} ny={ny}: L2={errors['l2_error']:.2e}")
            else:
                print(f"[{i+1}/{total}] {solver}+{limiter} ny={ny}: FAILED - {err_msg}")

    return pd.DataFrame(results)


def compute_convergence_order(h_values: np.ndarray, errors: np.ndarray) -> tuple:
    """Compute convergence order via linear regression"""
    log_h = np.log(h_values)
    log_e = np.log(errors)
    coeffs = np.polyfit(log_h, log_e, 1)
    p = coeffs[0]
    C = np.exp(coeffs[1])
    return p, C


def setup_axis_style(ax):
    """Set axis style (inward tick marks)"""
    ax.tick_params(axis='both', which='major', direction='in', length=6, width=1,
                   labelsize=10, top=True, right=True)
    ax.tick_params(axis='both', which='minor', direction='in', length=3, width=0.8,
                   top=True, right=True)

    for spine in ax.spines.values():
        spine.set_linewidth(1)


def plot_all_schemes(df: pd.DataFrame, save_file: str = None):
    """Plot grid convergence for all schemes with distinct line styles"""

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))

    # Define styles for each scheme
    # Colors: projection=blue, simple=red, piso=green
    # Line styles: none=solid, minmod=dashed, superbee=dashdot, vanleer=dotted, mc=loosely dashed
    solver_colors = {
        'projection': '#1f77b4',  # blue
        'simple': '#d62728',       # red
        'piso': '#2ca02c'          # green
    }

    limiter_styles = {
        'none': '-',
        'minmod': '--',
        'superbee': '-.',
        'vanleer': ':',
        'mc': (0, (5, 1))  # loosely dashed
    }

    limiter_markers = {
        'none': 'o',
        'minmod': 's',
        'superbee': '^',
        'vanleer': 'D',
        'mc': 'v'
    }

    # Get unique schemes
    schemes = df.groupby(['solver', 'limiter']).groups.keys()

    convergence_orders = {}

    for solver, limiter in sorted(schemes):
        subset = df[(df['solver'] == solver) & (df['limiter'] == limiter)]
        subset = subset.sort_values('h')

        h_values = subset['h'].values
        l2_errors = subset['l2_error'].values

        if len(h_values) < 2:
            continue

        p, C = compute_convergence_order(h_values, l2_errors)
        convergence_orders[(solver, limiter)] = p

        label = f'{solver.capitalize()} + {limiter} (p={p:.2f})'
        color = solver_colors.get(solver, 'gray')
        linestyle = limiter_styles.get(limiter, '-')
        marker = limiter_markers.get(limiter, 'o')

        ax.loglog(h_values, l2_errors, color=color, linestyle=linestyle,
                  marker=marker, markersize=5, linewidth=1, label=label)

    # Reference lines
    h_all = df['h'].values
    h_ref = np.array([h_all.min() * 0.8, h_all.max() * 1.2])
    l2_median = df['l2_error'].median()

    ax.loglog(h_ref, l2_median * (h_ref / h_ref.mean())**1, 'k--',
              linewidth=1, alpha=0.4, label='1st order')
    ax.loglog(h_ref, l2_median * (h_ref / h_ref.mean())**2, 'k:',
              linewidth=1, alpha=0.4, label='2nd order')

    ax.set_xlabel(r'Grid spacing $h$ [m]', fontsize=11)
    ax.set_ylabel(r'$L_2$ error [m/s]', fontsize=11)
    ax.set_title('Grid Convergence: All Solver/Limiter Combinations', fontsize=12)
    ax.legend(loc='lower right', fontsize=7, framealpha=0.9, ncol=2)

    setup_axis_style(ax)
    ax.set_box_aspect(1)

    plt.tight_layout()

    # Save
    if save_file is None:
        save_file = 'channel_convergence_all.svg'

    base_path = os.path.splitext(save_file)[0]
    for ext in ['.svg', '.pdf', '.png']:
        out_file = base_path + ext
        dpi = 300 if ext == '.png' else None
        plt.savefig(out_file, bbox_inches='tight', dpi=dpi)
        print(f"Saved: {out_file}")

    return convergence_orders


def print_convergence_table(df: pd.DataFrame):
    """Print convergence results as table"""

    print("\n" + "=" * 80)
    print("Grid Convergence Study Results - All Schemes")
    print("=" * 80)

    schemes = df.groupby(['solver', 'limiter']).groups.keys()

    for solver, limiter in sorted(schemes):
        subset = df[(df['solver'] == solver) & (df['limiter'] == limiter)]
        subset = subset.sort_values('ny')

        if len(subset) < 2:
            continue

        h_values = subset['h'].values
        l2_errors = subset['l2_error'].values
        p, _ = compute_convergence_order(h_values, l2_errors)

        print(f"\n{solver.upper()} + {limiter.upper()} (order = {p:.2f})")
        print("-" * 50)
        print(f"{'ny':>6} | {'h [m]':>12} | {'L2 error':>12}")
        print("-" * 50)

        for _, row in subset.iterrows():
            print(f"{int(row['ny']):>6} | {row['h']:>12.6f} | {row['l2_error']:>12.6e}")

    print("\n" + "=" * 80)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Grid convergence study for all solver/limiter combinations')
    parser.add_argument('--ny', nargs='+', type=int,
                        default=[8, 16, 32, 64],
                        help='List of ny values')
    parser.add_argument('--U_max', type=float, default=0.015,
                        help='Maximum velocity [m/s]')
    parser.add_argument('--end_time', type=float, default=5.0,
                        help='Simulation end time [s]')
    parser.add_argument('--build_dir', type=str, default='build',
                        help='Build directory path')
    parser.add_argument('--save', type=str, default=None,
                        help='Save figure to file')
    parser.add_argument('--load', type=str, default=None,
                        help='Load existing results from CSV')
    parser.add_argument('--sequential', action='store_true',
                        help='Run simulations sequentially (not in parallel)')

    args = parser.parse_args()

    # Define schemes to test
    # Solver: projection, simple, piso
    # Limiter: none, minmod, superbee, vanleer, mc
    schemes = [
        ('projection', 'none'),
        ('projection', 'minmod'),
        ('projection', 'vanleer'),
        ('simple', 'none'),
        ('simple', 'minmod'),
        ('simple', 'vanleer'),
        ('piso', 'none'),
        ('piso', 'minmod'),
        ('piso', 'vanleer'),
    ]

    if args.load:
        df = pd.read_csv(args.load)
    else:
        df = run_convergence_study(
            schemes, args.ny, args.U_max, args.end_time,
            args.build_dir, parallel=not args.sequential
        )

        # Save results
        csv_file = os.path.join(args.build_dir, 'output', 'convergence_all_results.csv')
        os.makedirs(os.path.dirname(csv_file), exist_ok=True)
        df.to_csv(csv_file, index=False)
        print(f"\nResults saved to {csv_file}")

    # Print table
    print_convergence_table(df)

    # Plot
    save_file = args.save or 'docs/images/channel_convergence.svg'
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    convergence_orders = plot_all_schemes(df, save_file)

    print("\n=== Convergence Orders ===")
    for (solver, limiter), order in sorted(convergence_orders.items()):
        print(f"  {solver}+{limiter}: {order:.2f}")
