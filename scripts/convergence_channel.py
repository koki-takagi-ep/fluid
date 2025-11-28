#!/usr/bin/env python3
"""
Channel Flow (Hagen-Poiseuille) 格子収束性解析

L2誤差ノルムを計算し、格子解像度に対する収束次数を評価する。

L2誤差の定義:
    ||e||_L2 = sqrt( (1/|Omega|) * integral |u_num - u_exact|^2 dOmega )
            ≈ sqrt( (1/N) * sum_i (u_num,i - u_exact,i)^2 )

格子収束性:
    ||e||_L2 = C * h^p

    ここで h は格子幅、p は収束次数（スキームの精度次数）

    両対数グラフで:
    log(||e||_L2) = log(C) + p * log(h)

    傾き p が収束次数を示す。
    - 1次精度スキーム: p ≈ 1
    - 2次精度スキーム: p ≈ 2
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import subprocess
import sys


def poiseuille_velocity(y: np.ndarray, H: float, U_max: float) -> np.ndarray:
    """Hagen-Poiseuille流れの速度分布（理論解）

    u(y) = U_max * (1 - (2y/H - 1)^2)
    """
    eta = 2.0 * y / H - 1.0
    return U_max * (1.0 - eta**2)


def load_simulation_data(output_dir: str):
    """シミュレーション結果を読み込む"""
    data_dir = os.path.join(output_dir, "data")
    if not os.path.exists(data_dir):
        data_dir = output_dir

    # メタデータ
    metadata_file = os.path.join(data_dir, "metadata.csv")
    df_meta = pd.read_csv(metadata_file)
    meta = dict(zip(df_meta['parameter'], df_meta['value']))
    nx, ny = int(meta['nx']), int(meta['ny'])
    lx, ly = float(meta['lx']), float(meta['ly'])
    dx, dy = lx / nx, ly / ny

    # 最終時刻のフィールドデータ
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
    """L2誤差ノルムを計算

    Args:
        data: シミュレーションデータ
        U_max: 最大速度
        x_position: 誤差を計算するx位置（無次元、0-1）

    Returns:
        dict: L2誤差、Linf誤差、格子幅など
    """
    nx, ny = data['nx'], data['ny']
    lx, ly = data['lx'], data['ly']
    h = data['dy']  # y方向の格子幅（収束解析用）

    # 指定x位置でのプロファイル
    i = int(round(x_position * (nx - 1)))
    i = min(max(i, 0), nx - 1)

    y_sim = data['y'][i, :]  # 次元付きy座標
    u_sim = data['u'][i, :]

    # 理論解
    u_exact = poiseuille_velocity(y_sim, ly, U_max)

    # 誤差
    error = u_sim - u_exact

    # L2ノルム: sqrt( (1/N) * sum(e^2) )
    l2_error = np.sqrt(np.mean(error**2))

    # L∞ノルム（最大誤差）
    linf_error = np.max(np.abs(error))

    # 相対L2誤差
    l2_relative = l2_error / U_max

    return {
        'h': h,
        'ny': ny,
        'l2_error': l2_error,
        'l2_relative': l2_relative,
        'linf_error': linf_error,
        'linf_relative': linf_error / U_max
    }


def run_convergence_study(ny_list: list, U_max: float = 0.015,
                          end_time: float = 5.0, build_dir: str = '.'):
    """複数の格子解像度でシミュレーションを実行し収束性を解析

    Args:
        ny_list: y方向格子数のリスト（例: [8, 16, 32, 64, 128]）
        U_max: 最大速度 [m/s]
        end_time: シミュレーション終了時刻 [s]
        build_dir: ビルドディレクトリのパス

    Returns:
        DataFrame: 各格子での誤差データ
    """
    results = []

    for ny in ny_list:
        # アスペクト比を維持（Lx/Ly = 10）
        nx = ny * 4  # チャネル長さ/高さ = 10 なので概ね4倍

        print(f"\n=== Running simulation: nx={nx}, ny={ny} ===")

        # 出力ディレクトリを設定
        output_dir = f"output/convergence_ny{ny}"

        # 既存の出力を削除
        subprocess.run(f"rm -rf {output_dir}", shell=True, cwd=build_dir)

        # シミュレーション実行
        cmd = f"./channel_flow {nx} {ny} {U_max} {end_time}"
        print(f"Command: {cmd}")

        # channel_flow の出力先を一時的に変更するため、環境変数やコードを修正する必要がある
        # 代わりに、標準出力ディレクトリを使用してからリネーム
        result = subprocess.run(
            cmd, shell=True, cwd=build_dir,
            capture_output=True, text=True
        )

        if result.returncode != 0:
            print(f"Error running simulation: {result.stderr}")
            continue

        # 出力ディレクトリをリネーム
        subprocess.run(
            f"mv output/channel_projection {output_dir}",
            shell=True, cwd=build_dir
        )

        # 誤差計算
        try:
            data = load_simulation_data(os.path.join(build_dir, output_dir))
            errors = compute_l2_error(data, U_max)
            errors['nx'] = nx
            results.append(errors)
            print(f"  h = {errors['h']:.6f}, L2 error = {errors['l2_error']:.6e}")
        except Exception as e:
            print(f"Error computing L2 error: {e}")
            continue

    return pd.DataFrame(results)


def compute_convergence_order(h_values: np.ndarray, errors: np.ndarray) -> tuple:
    """収束次数を計算（最小二乗法による線形回帰）

    log(error) = log(C) + p * log(h)

    Returns:
        (p, C): 収束次数 p と定数 C
    """
    log_h = np.log(h_values)
    log_e = np.log(errors)

    # 線形回帰
    coeffs = np.polyfit(log_h, log_e, 1)
    p = coeffs[0]  # 傾き = 収束次数
    C = np.exp(coeffs[1])  # 切片

    return p, C


def plot_convergence(df: pd.DataFrame, save_file: str = None):
    """格子収束性プロットを作成（両対数グラフ）"""

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    h_values = df['h'].values
    l2_errors = df['l2_error'].values
    linf_errors = df['linf_error'].values

    # 収束次数を計算
    p_l2, C_l2 = compute_convergence_order(h_values, l2_errors)
    p_linf, C_linf = compute_convergence_order(h_values, linf_errors)

    # 参照線用のデータ
    h_ref = np.array([h_values.min() * 0.8, h_values.max() * 1.2])

    # === 左: L2誤差 ===
    ax = axes[0]
    ax.loglog(h_values, l2_errors, 'bo-', markersize=8, linewidth=2,
              label=f'$L_2$ error (slope = {p_l2:.2f})')

    # 1次・2次の参照線
    ax.loglog(h_ref, C_l2 * h_ref**1, 'k--', linewidth=1, alpha=0.5,
              label='1st order ($h^1$)')
    ax.loglog(h_ref, C_l2 * h_ref**2, 'k:', linewidth=1, alpha=0.5,
              label='2nd order ($h^2$)')

    ax.set_xlabel(r'Grid spacing $h$ [m]', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'$L_2$ error [m/s]', fontsize=11, fontweight='bold')
    ax.set_title(r'$L_2$ Error Convergence', fontsize=12, fontweight='bold')
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, which='both', alpha=0.3, linestyle='--')

    # 収束次数をテキストで表示
    ax.text(0.05, 0.95, f'Convergence order: {p_l2:.2f}',
            transform=ax.transAxes, fontsize=11, fontweight='bold',
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # === 右: L∞誤差 ===
    ax = axes[1]
    ax.loglog(h_values, linf_errors, 'rs-', markersize=8, linewidth=2,
              label=f'$L_\\infty$ error (slope = {p_linf:.2f})')

    ax.loglog(h_ref, C_linf * h_ref**1, 'k--', linewidth=1, alpha=0.5,
              label='1st order ($h^1$)')
    ax.loglog(h_ref, C_linf * h_ref**2, 'k:', linewidth=1, alpha=0.5,
              label='2nd order ($h^2$)')

    ax.set_xlabel(r'Grid spacing $h$ [m]', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'$L_\infty$ error [m/s]', fontsize=11, fontweight='bold')
    ax.set_title(r'$L_\infty$ Error Convergence', fontsize=12, fontweight='bold')
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, which='both', alpha=0.3, linestyle='--')

    ax.text(0.05, 0.95, f'Convergence order: {p_linf:.2f}',
            transform=ax.transAxes, fontsize=11, fontweight='bold',
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # 全体タイトル
    fig.suptitle('Grid Convergence Study: Hagen-Poiseuille Flow',
                 fontsize=14, fontweight='bold', y=1.02)

    plt.tight_layout()

    # 保存
    if save_file is None:
        save_file = 'channel_convergence.svg'

    base_path = os.path.splitext(save_file)[0]
    for ext in ['.svg', '.pdf', '.png']:
        out_file = base_path + ext
        dpi = 300 if ext == '.png' else None
        plt.savefig(out_file, bbox_inches='tight', dpi=dpi)
        print(f"Figure saved to {out_file}")

    return p_l2, p_linf


def print_convergence_table(df: pd.DataFrame):
    """収束性の結果をテーブル形式で出力"""

    print("\n" + "=" * 70)
    print("Grid Convergence Study Results")
    print("=" * 70)
    print(f"{'ny':>6} | {'h [m]':>12} | {'L2 error':>12} | {'L∞ error':>12} | {'L2 ratio':>10}")
    print("-" * 70)

    prev_l2 = None
    for _, row in df.iterrows():
        ratio_str = ""
        if prev_l2 is not None:
            ratio = prev_l2 / row['l2_error']
            ratio_str = f"{ratio:.2f}"
        print(f"{int(row['ny']):>6} | {row['h']:>12.6f} | {row['l2_error']:>12.6e} | "
              f"{row['linf_error']:>12.6e} | {ratio_str:>10}")
        prev_l2 = row['l2_error']

    print("=" * 70)

    # 収束次数
    h_values = df['h'].values
    l2_errors = df['l2_error'].values
    p_l2, _ = compute_convergence_order(h_values, l2_errors)

    print(f"\nL2 convergence order: {p_l2:.2f}")
    print("(Expected: 1.0 for 1st-order upwind, 2.0 for 2nd-order central)")
    print("=" * 70)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Grid convergence study for channel flow')
    parser.add_argument('--ny', nargs='+', type=int,
                        default=[8, 16, 32, 64],
                        help='List of ny values (default: 8 16 32 64)')
    parser.add_argument('--U_max', type=float, default=0.015,
                        help='Maximum velocity [m/s]')
    parser.add_argument('--end_time', type=float, default=5.0,
                        help='Simulation end time [s]')
    parser.add_argument('--build_dir', type=str, default='.',
                        help='Build directory path')
    parser.add_argument('--save', type=str, default=None,
                        help='Save figure to file')
    parser.add_argument('--load', type=str, default=None,
                        help='Load existing results from CSV')

    args = parser.parse_args()

    if args.load:
        # 既存の結果を読み込み
        df = pd.read_csv(args.load)
    else:
        # シミュレーション実行
        df = run_convergence_study(
            args.ny, args.U_max, args.end_time, args.build_dir
        )

        # 結果を保存
        csv_file = os.path.join(args.build_dir, 'convergence_results.csv')
        df.to_csv(csv_file, index=False)
        print(f"\nResults saved to {csv_file}")

    # テーブル出力
    print_convergence_table(df)

    # プロット
    save_file = args.save or os.path.join(args.build_dir, 'output', 'channel_convergence.svg')
    os.makedirs(os.path.dirname(save_file), exist_ok=True)
    p_l2, p_linf = plot_convergence(df, save_file)

    print(f"\nFinal convergence orders:")
    print(f"  L2:   {p_l2:.2f}")
    print(f"  L∞:   {p_linf:.2f}")
