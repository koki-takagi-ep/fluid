#!/usr/bin/env python3
"""
TVDスキームのGhiaベンチマーク比較
各リミッター（Upwind, Minmod, Superbee, Van Leer, MC）の精度を比較
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# Ghia et al. (1982) benchmark data for Re=100
GHIA_RE100 = {
    'u_y': np.array([0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813,
                     0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609,
                     0.9688, 0.9766, 1.0000]),
    'u': np.array([0.0000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150,
                   -0.15662, -0.21090, -0.20581, -0.13641, 0.00332, 0.23151,
                   0.68717, 0.73722, 0.78871, 0.84123, 1.00000]),
    'v_x': np.array([0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266,
                     0.2344, 0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531,
                     0.9609, 0.9688, 1.0000]),
    'v': np.array([0.00000, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077,
                   0.17507, 0.17527, 0.05454, -0.24533, -0.22445, -0.16914,
                   -0.10313, -0.08864, -0.07391, -0.05906, 0.00000])
}


def load_simulation_data(output_dir):
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

    # 最終時刻のフィールドデータ
    files = sorted(glob.glob(os.path.join(data_dir, "field_*.csv")))
    if not files:
        raise FileNotFoundError("No field files found")

    last_file = files[-1]
    df = pd.read_csv(last_file, comment='#')

    x = df['x'].values.reshape(nx, ny)
    y = df['y'].values.reshape(nx, ny)
    u = df['u'].values.reshape(nx, ny)
    v = df['v'].values.reshape(nx, ny)

    return {
        'nx': nx, 'ny': ny, 'lx': lx, 'ly': ly,
        'x': x, 'y': y, 'u': u, 'v': v
    }


def extract_centerline_data(data, U_lid):
    """中心線データを抽出"""
    nx, ny = data['nx'], data['ny']
    lx, ly = data['lx'], data['ly']

    # 垂直中心線のu速度
    i_center = nx // 2
    y_line = data['y'][i_center, :] / ly
    u_line = data['u'][i_center, :] / U_lid

    # 水平中心線のv速度
    j_center = ny // 2
    x_line = data['x'][:, j_center] / lx
    v_line = data['v'][:, j_center] / U_lid

    return {
        'y': y_line, 'u': u_line,
        'x': x_line, 'v': v_line
    }


def compute_rms_error(sim_data, ghia_data, U_lid):
    """Ghiaデータとの二乗平均平方根誤差を計算"""
    centerline = extract_centerline_data(sim_data, U_lid)

    # u速度の誤差（補間）
    u_interp = np.interp(GHIA_RE100['u_y'], centerline['y'], centerline['u'])
    u_rms = np.sqrt(np.mean((u_interp - GHIA_RE100['u'])**2))

    # v速度の誤差（補間）
    v_interp = np.interp(GHIA_RE100['v_x'], centerline['x'], centerline['v'])
    v_rms = np.sqrt(np.mean((v_interp - GHIA_RE100['v'])**2))

    return u_rms, v_rms


def setup_axes_style(ax):
    """軸のスタイルを設定（ユーザー指定のデザイン）"""
    # 目盛りを内向きに
    ax.tick_params(axis='both', which='major', direction='in', length=6, width=1,
                   labelsize=10, top=True, right=True)
    ax.tick_params(axis='both', which='minor', direction='in', length=3, width=0.8,
                   top=True, right=True)

    # 軸の線幅
    for spine in ax.spines.values():
        spine.set_linewidth(1)


def main():
    build_dir = os.path.dirname(os.path.abspath(__file__))
    build_dir = os.path.join(os.path.dirname(build_dir), 'build')

    # リミッター設定（描画順序を考慮）
    limiters = {
        'Upwind (1st-order)': 'tvd_upwind',
        'Minmod': 'tvd_minmod',
        'Superbee': 'tvd_superbee',
        'Van Leer': 'tvd_vanleer',
        'MC': 'tvd_mc'
    }

    # 色: 実線は濃い色、破線は区別しやすい色
    colors = {
        'Upwind (1st-order)': '#1f77b4',  # blue
        'Minmod': '#d62728',               # red
        'Superbee': '#2ca02c',             # green
        'Van Leer': '#7f7f7f',             # gray (dotted)
        'MC': '#9467bd'                    # purple
    }

    linestyles = {
        'Upwind (1st-order)': '-',
        'Minmod': '--',
        'Superbee': '-.',
        'Van Leer': ':',
        'MC': '-'
    }

    # zorder: 破線を前面に（大きい値が前面）
    zorders = {
        'Upwind (1st-order)': 2,
        'Minmod': 5,
        'Superbee': 4,
        'Van Leer': 6,
        'MC': 3
    }

    U_lid = 0.01  # m/s

    # フィギュアサイズ（ユーザー指定に近いサイズ）
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    rms_errors = {}

    # 左: u-velocity along vertical centerline
    ax = axes[0]

    # まずシミュレーション結果をプロット（背面）
    for name, dirname in limiters.items():
        path = os.path.join(build_dir, f'output/{dirname}')
        if os.path.exists(path):
            try:
                data = load_simulation_data(path)
                centerline = extract_centerline_data(data, U_lid)
                ax.plot(centerline['u'], centerline['y'],
                       color=colors[name], linestyle=linestyles[name],
                       linewidth=1, label=name, zorder=zorders[name])

                u_rms, v_rms = compute_rms_error(data, GHIA_RE100, U_lid)
                rms_errors[name] = (u_rms, v_rms)
            except Exception as e:
                print(f"Error loading {name}: {e}")

    # Ghiaデータを最前面に
    ax.plot(GHIA_RE100['u'], GHIA_RE100['u_y'], 'ko', markersize=6,
            label='Ghia et al. (1982)', zorder=10)

    ax.set_xlabel(r'$u/U_{lid}$', fontsize=11)
    ax.set_ylabel(r'$y/L$', fontsize=11)
    ax.set_title(r'$u$-velocity along vertical centerline', fontsize=12)
    ax.legend(loc='lower right', fontsize=8, framealpha=0.9)
    ax.set_xlim(-0.4, 1.0)
    ax.set_ylim(0, 1)

    # マイナー目盛り
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    setup_axes_style(ax)

    # 右: v-velocity along horizontal centerline
    ax = axes[1]

    for name, dirname in limiters.items():
        path = os.path.join(build_dir, f'output/{dirname}')
        if os.path.exists(path):
            try:
                data = load_simulation_data(path)
                centerline = extract_centerline_data(data, U_lid)
                ax.plot(centerline['x'], centerline['v'],
                       color=colors[name], linestyle=linestyles[name],
                       linewidth=1, label=name, zorder=zorders[name])
            except Exception as e:
                print(f"Error loading {name}: {e}")

    # Ghiaデータを最前面に
    ax.plot(GHIA_RE100['v_x'], GHIA_RE100['v'], 'ko', markersize=6,
            label='Ghia et al. (1982)', zorder=10)

    ax.set_xlabel(r'$x/L$', fontsize=11)
    ax.set_ylabel(r'$v/U_{lid}$', fontsize=11)
    ax.set_title(r'$v$-velocity along horizontal centerline', fontsize=12)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, 0.4)

    # マイナー目盛り
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    setup_axes_style(ax)

    fig.suptitle('TVD Limiter Comparison: Cavity Flow Validation (Re = 100)',
                 fontsize=13, y=0.98)
    plt.tight_layout()

    # RMS誤差を表示
    if rms_errors:
        print("\n=== RMS Error vs Ghia et al. (1982) ===")
        print(f"{'Limiter':<25} {'u-RMS':<12} {'v-RMS':<12}")
        print("-" * 49)
        for name, (u_rms, v_rms) in rms_errors.items():
            print(f"{name:<25} {u_rms:<12.6f} {v_rms:<12.6f}")

    # 保存
    output_dir = os.path.join(os.path.dirname(build_dir), 'docs', 'images')
    os.makedirs(output_dir, exist_ok=True)

    for ext in ['.svg', '.pdf', '.png']:
        out_file = os.path.join(output_dir, f'tvd_validation_comparison{ext}')
        dpi = 300 if ext == '.png' else None
        plt.savefig(out_file, bbox_inches='tight', dpi=dpi)
        print(f"Saved: {out_file}")


if __name__ == '__main__':
    main()
