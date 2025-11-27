#!/usr/bin/env python3
"""
収束確認用のグラフを作成するスクリプト

時間発展に伴う速度場の変化を確認
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


def get_data_dir(output_dir: str) -> str:
    """dataディレクトリのパスを取得"""
    data_dir = os.path.join(output_dir, "data")
    if os.path.exists(data_dir):
        return data_dir
    return output_dir


def get_figures_dir(output_dir: str) -> str:
    """figuresディレクトリのパスを取得"""
    figures_dir = os.path.join(output_dir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir, exist_ok=True)
    return figures_dir


def load_metadata(output_dir: str) -> dict:
    """メタデータを読み込む"""
    data_dir = get_data_dir(output_dir)
    metadata_file = os.path.join(data_dir, "metadata.csv")
    df = pd.read_csv(metadata_file)
    return dict(zip(df['parameter'], df['value']))


def load_field(filepath: str, nx: int, ny: int) -> dict:
    """フィールドデータを読み込む"""
    with open(filepath, 'r') as f:
        first_line = f.readline()
        time = float(first_line.split('=')[1])

    df = pd.read_csv(filepath, comment='#')

    u = df['u'].values.reshape(nx, ny)
    v = df['v'].values.reshape(nx, ny)
    mag = df['magnitude'].values.reshape(nx, ny)

    return {'time': time, 'u': u, 'v': v, 'magnitude': mag}


def setup_axis_style(ax, xlabel='', ylabel='', title=''):
    """軸のスタイルを設定"""
    ax.set_xlabel(xlabel, fontsize=11, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=11, fontweight='bold')
    if title:
        ax.set_title(title, fontsize=12, fontweight='bold')

    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(which='major', direction='out', length=6, width=1, labelsize=10)
    ax.tick_params(which='minor', direction='out', length=3, width=0.5)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))


def plot_convergence(output_dir: str, save_file: str = None, vel_unit='mm/s'):
    """収束確認グラフを作成"""
    metadata = load_metadata(output_dir)
    nx, ny = int(metadata['nx']), int(metadata['ny'])

    data_dir = get_data_dir(output_dir)
    files = sorted(glob.glob(os.path.join(data_dir, "field_*.csv")))

    if len(files) < 2:
        print("Need at least 2 field files for convergence check")
        return

    # 単位変換係数
    V_scale = 1000.0 if vel_unit == 'mm/s' else 1.0

    # データ収集
    times = []
    max_velocities = []
    centerline_u = []  # 中心線でのu速度
    changes = []       # 前ステップからの変化量

    prev_mag = None

    for f in files:
        field = load_field(f, nx, ny)
        times.append(field['time'])
        max_velocities.append(np.max(field['magnitude']) * V_scale)

        # 中心線での速度（出口付近）
        i_exit = int(nx * 0.9)  # 出口近く
        j_mid = ny // 2
        centerline_u.append(field['u'][i_exit, j_mid] * V_scale)

        # 変化量
        if prev_mag is not None:
            change = np.max(np.abs(field['magnitude'] - prev_mag)) * V_scale
            changes.append(change)
        prev_mag = field['magnitude'].copy()

    # プロット作成
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    # 左上: 最大速度の時間発展
    ax = axes[0, 0]
    ax.plot(times, max_velocities, 'b-', linewidth=1.5, marker='o', markersize=4)
    setup_axis_style(ax, xlabel='Time (s)', ylabel=f'Max velocity ({vel_unit})',
                     title='Maximum Velocity vs Time')
    ax.grid(True, alpha=0.3, linestyle='--')

    # 右上: 中心線速度の時間発展
    ax = axes[0, 1]
    ax.plot(times, centerline_u, 'r-', linewidth=1.5, marker='s', markersize=4)
    setup_axis_style(ax, xlabel='Time (s)', ylabel=f'Centerline $u$ ({vel_unit})',
                     title='Centerline Velocity (near exit)')
    ax.grid(True, alpha=0.3, linestyle='--')

    # 左下: 変化量の時間発展（収束確認）
    ax = axes[1, 0]
    if len(changes) > 0:
        ax.semilogy(times[1:], changes, 'g-', linewidth=1.5, marker='^', markersize=4)
        ax.set_xlabel('Time (s)', fontsize=11, fontweight='bold')
        ax.set_ylabel(f'Max change ({vel_unit})', fontsize=11, fontweight='bold')
        ax.set_title('Velocity Change (Convergence)', fontsize=12, fontweight='bold')
        ax.tick_params(which='major', direction='out', length=6, width=1, labelsize=10)
        ax.grid(True, alpha=0.3, linestyle='--', which='both')
    else:
        ax.text(0.5, 0.5, 'Not enough data', ha='center', va='center', transform=ax.transAxes)

    # 右下: 定常状態への収束率
    ax = axes[1, 1]
    if len(max_velocities) > 1:
        final_vel = max_velocities[-1]
        convergence_ratio = [abs(v - final_vel) / final_vel * 100 for v in max_velocities]
        ax.semilogy(times, convergence_ratio, 'm-', linewidth=1.5, marker='d', markersize=4)
        ax.set_xlabel('Time (s)', fontsize=11, fontweight='bold')
        ax.set_ylabel('Deviation from final (%)', fontsize=11, fontweight='bold')
        ax.set_title('Convergence to Steady State', fontsize=12, fontweight='bold')
        ax.tick_params(which='major', direction='out', length=6, width=1, labelsize=10)
        ax.grid(True, alpha=0.3, linestyle='--', which='both')

    plt.tight_layout()

    if save_file is None:
        save_file = os.path.join(get_figures_dir(output_dir), 'convergence.pdf')

    plt.savefig(save_file, bbox_inches='tight')
    print(f"Convergence plot saved to {save_file}")

    # 収束情報を表示
    print(f"\n=== Convergence Summary ===")
    print(f"Initial max velocity: {max_velocities[0]:.4f} {vel_unit}")
    print(f"Final max velocity:   {max_velocities[-1]:.4f} {vel_unit}")
    if len(changes) > 0:
        print(f"Final change rate:    {changes[-1]:.6f} {vel_unit}")
    print(f"Simulation time:      {times[-1]:.2f} s")
    print(f"===========================")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot convergence of fluid simulation')
    parser.add_argument('output_dir', help='Output directory containing CSV files')
    parser.add_argument('--save', type=str, default=None, help='Save figure to file')

    args = parser.parse_args()

    plot_convergence(args.output_dir, args.save)
