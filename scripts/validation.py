#!/usr/bin/env python3
"""
Lid-Driven Cavity Flow の検証スクリプト

Ghia, Ghia & Shin (1982) のベンチマークデータとの比較
"High-Re solutions for incompressible flow using the Navier-Stokes equations
 and a multigrid method", J. Comput. Phys., 48, 387-411

参照データ: 垂直・水平中心線での速度プロファイル
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


# Ghia et al. (1982) のベンチマークデータ
# 垂直中心線 (x=0.5) での u 速度
GHIA_DATA_U = {
    # Re = 100
    100: {
        'y': [0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813,
              0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609,
              0.9688, 0.9766, 1.0000],
        'u': [0.0000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662,
              -0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73722,
              0.78871, 0.84123, 1.0000]
    },
    # Re = 400
    400: {
        'y': [0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813,
              0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609,
              0.9688, 0.9766, 1.0000],
        'u': [0.0000, -0.08186, -0.09266, -0.10338, -0.14612, -0.24299, -0.32726,
              -0.17119, -0.11477, 0.02135, 0.16256, 0.29093, 0.55892, 0.61756,
              0.68439, 0.75837, 1.0000]
    },
    # Re = 1000
    1000: {
        'y': [0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813,
              0.4531, 0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609,
              0.9688, 0.9766, 1.0000],
        'u': [0.0000, -0.18109, -0.20196, -0.22220, -0.29730, -0.38289, -0.27805,
              -0.10648, -0.06080, 0.05702, 0.18719, 0.33304, 0.46604, 0.51117,
              0.57492, 0.65928, 1.0000]
    }
}

# 水平中心線 (y=0.5) での v 速度
GHIA_DATA_V = {
    # Re = 100
    100: {
        'x': [0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266,
              0.2344, 0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531,
              0.9609, 0.9688, 1.0000],
        'v': [0.0000, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507,
              0.17527, 0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864,
              -0.07391, -0.05906, 0.0000]
    },
    # Re = 400
    400: {
        'x': [0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266,
              0.2344, 0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531,
              0.9609, 0.9688, 1.0000],
        'v': [0.0000, 0.18360, 0.19713, 0.20920, 0.22965, 0.28124, 0.30203,
              0.30174, 0.05186, -0.38598, -0.44993, -0.23827, -0.22847, -0.19254,
              -0.15663, -0.12146, 0.0000]
    },
    # Re = 1000
    1000: {
        'x': [0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266,
              0.2344, 0.5000, 0.8047, 0.8594, 0.9063, 0.9453, 0.9531,
              0.9609, 0.9688, 1.0000],
        'v': [0.0000, 0.27485, 0.29012, 0.30353, 0.32627, 0.37095, 0.33075,
              0.32235, 0.02526, -0.31966, -0.42665, -0.51550, -0.39188, -0.33714,
              -0.27669, -0.21388, 0.0000]
    }
}


def get_data_dir(output_dir: str) -> str:
    """dataディレクトリのパスを取得（新旧両方の構造に対応）"""
    data_dir = os.path.join(output_dir, "data")
    if os.path.exists(data_dir):
        return data_dir
    return output_dir  # 旧構造の場合


def get_figures_dir(output_dir: str) -> str:
    """figuresディレクトリのパスを取得（なければ作成）"""
    figures_dir = os.path.join(output_dir, "figures")
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir, exist_ok=True)
    return figures_dir


def load_simulation_data(output_dir: str):
    """シミュレーション結果を読み込む"""
    data_dir = get_data_dir(output_dir)

    # メタデータ
    metadata_file = os.path.join(data_dir, "metadata.csv")
    df_meta = pd.read_csv(metadata_file)
    meta = dict(zip(df_meta['parameter'], df_meta['value']))
    nx, ny = int(meta['nx']), int(meta['ny'])
    lx, ly = float(meta['lx']), float(meta['ly'])

    # 最終時刻のフィールドデータを取得
    import glob
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


def extract_centerline_profiles(data: dict, U_lid: float = 1.0):
    """中心線での速度プロファイルを抽出（無次元化）"""
    nx, ny = data['nx'], data['ny']
    lx, ly = data['lx'], data['ly']

    # 垂直中心線 (x = 0.5*lx) での u 速度
    i_mid = nx // 2
    y_vertical = data['y'][i_mid, :] / ly  # 無次元化
    u_vertical = data['u'][i_mid, :] / U_lid  # 無次元化

    # 水平中心線 (y = 0.5*ly) での v 速度
    j_mid = ny // 2
    x_horizontal = data['x'][:, j_mid] / lx  # 無次元化
    v_horizontal = data['v'][:, j_mid] / U_lid  # 無次元化

    return {
        'y_vertical': y_vertical,
        'u_vertical': u_vertical,
        'x_horizontal': x_horizontal,
        'v_horizontal': v_horizontal
    }


def setup_axis_style(ax):
    """軸のスタイルを設定（内向き目盛り）"""
    ax.tick_params(axis='both', which='major', direction='in', length=6, width=1,
                   labelsize=10, top=True, right=True)
    ax.tick_params(axis='both', which='minor', direction='in', length=3, width=0.8,
                   top=True, right=True)

    for spine in ax.spines.values():
        spine.set_linewidth(1)


def plot_validation(output_dir: str, Re: int, U_lid: float, save_file: str = None):
    """検証プロットを作成（単一ケース、正方形グラフ）"""

    # デフォルトの出力先はfiguresディレクトリ
    if save_file is None:
        save_file = os.path.join(get_figures_dir(output_dir), f'validation_Re{Re}.pdf')

    if Re not in GHIA_DATA_U:
        print(f"Warning: No Ghia data available for Re={Re}")
        available = list(GHIA_DATA_U.keys())
        print(f"Available Re: {available}")
        return

    # シミュレーションデータを読み込み
    data = load_simulation_data(output_dir)
    profiles = extract_centerline_profiles(data, U_lid)

    # Ghiaのデータ
    ghia_u = GHIA_DATA_U[Re]
    ghia_v = GHIA_DATA_V[Re]

    # プロット作成（正方形の2パネル）
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # 左: 垂直中心線での u 速度
    ax = axes[0]
    ax.plot(profiles['u_vertical'], profiles['y_vertical'],
            'b-', linewidth=1, label='Present (CFD)', zorder=2)
    ax.plot(ghia_u['u'], ghia_u['y'],
            'ko', markersize=6, label='Ghia et al. (1982)', zorder=10)

    ax.set_xlabel(r'$u/U_{lid}$', fontsize=11)
    ax.set_ylabel(r'$y/L$', fontsize=11)
    ax.set_title(r'$u$-velocity along vertical centerline', fontsize=12)
    ax.set_xlim(-0.4, 1.0)
    ax.set_ylim(0, 1)
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.legend(loc='lower right', fontsize=8, framealpha=0.9)
    setup_axis_style(ax)
    ax.set_box_aspect(1)

    # 右: 水平中心線での v 速度
    ax = axes[1]
    ax.plot(profiles['x_horizontal'], profiles['v_horizontal'],
            'b-', linewidth=1, label='Present (CFD)', zorder=2)
    ax.plot(ghia_v['x'], ghia_v['v'],
            'ko', markersize=6, label='Ghia et al. (1982)', zorder=10)

    ax.set_xlabel(r'$x/L$', fontsize=11)
    ax.set_ylabel(r'$v/U_{lid}$', fontsize=11)
    ax.set_title(r'$v$-velocity along horizontal centerline', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, 0.4)
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
    setup_axis_style(ax)
    ax.set_box_aspect(1)

    # 全体タイトル
    fig.suptitle(f'Lid-Driven Cavity Flow Validation (Re = {Re})',
                 fontsize=13, y=0.98)

    plt.tight_layout()

    # 保存
    base_path = os.path.splitext(save_file)[0]
    for ext in ['.svg', '.pdf', '.png']:
        out_file = base_path + ext
        dpi = 300 if ext == '.png' else None
        plt.savefig(out_file, bbox_inches='tight', dpi=dpi)
        print(f"Saved: {out_file}")

    # 誤差を計算
    compute_error(profiles, ghia_u, ghia_v, Re)


def plot_multi_case_validation(output_dirs: list, labels: list, Re: int, U_lid: float,
                                save_file: str = None):
    """複数ケース（異なるセル数・スキーム）の検証プロットを作成

    Ghiaデータ: 黒ドット（1回だけ表示）
    数値解: 色付きライン（ケースごと）
    """

    if Re not in GHIA_DATA_U:
        print(f"Warning: No Ghia data available for Re={Re}")
        available = list(GHIA_DATA_U.keys())
        print(f"Available Re: {available}")
        return

    # Ghiaのデータ
    ghia_u = GHIA_DATA_U[Re]
    ghia_v = GHIA_DATA_V[Re]

    # カラーパレット
    colors = ['#1f77b4', '#d62728', '#2ca02c', '#ff7f0e', '#9467bd']
    linestyles = ['-', '--', '-.', ':']

    # プロット作成（正方形の2パネル）
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # 各ケースのデータを読み込みプロット
    for idx, (output_dir, label) in enumerate(zip(output_dirs, labels)):
        try:
            data = load_simulation_data(output_dir)
            profiles = extract_centerline_profiles(data, U_lid)

            color = colors[idx % len(colors)]
            ls = linestyles[idx % len(linestyles)]

            # 左: 垂直中心線での u 速度
            axes[0].plot(profiles['u_vertical'], profiles['y_vertical'],
                        color=color, linestyle=ls, linewidth=1, label=label, zorder=3+idx)

            # 右: 水平中心線での v 速度
            axes[1].plot(profiles['x_horizontal'], profiles['v_horizontal'],
                        color=color, linestyle=ls, linewidth=1, label=label, zorder=3+idx)

            # 誤差計算
            from scipy import interpolate
            f_u = interpolate.interp1d(profiles['y_vertical'], profiles['u_vertical'],
                                        kind='linear', fill_value='extrapolate')
            u_interp = f_u(ghia_u['y'])
            u_rms = np.sqrt(np.mean((u_interp - np.array(ghia_u['u']))**2))

            f_v = interpolate.interp1d(profiles['x_horizontal'], profiles['v_horizontal'],
                                        kind='linear', fill_value='extrapolate')
            v_interp = f_v(ghia_v['x'])
            v_rms = np.sqrt(np.mean((v_interp - np.array(ghia_v['v']))**2))

            print(f"{label}: u_RMS = {u_rms:.5f}, v_RMS = {v_rms:.5f}")

        except Exception as e:
            print(f"Error loading {output_dir}: {e}")
            continue

    # Ghiaデータを黒ドットでプロット（最後に描画して前面に）
    axes[0].plot(ghia_u['u'], ghia_u['y'],
                'ko', markersize=6, label='Ghia et al. (1982)', zorder=10)
    axes[1].plot(ghia_v['x'], ghia_v['v'],
                'ko', markersize=6, label='Ghia et al. (1982)', zorder=10)

    # 左パネルのスタイル
    ax = axes[0]
    ax.set_xlabel(r'$u/U_{lid}$', fontsize=11)
    ax.set_ylabel(r'$y/L$', fontsize=11)
    ax.set_title(r'$u$-velocity along vertical centerline', fontsize=12)
    ax.set_xlim(-0.4, 1.0)
    ax.set_ylim(0, 1)
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.legend(loc='lower right', fontsize=8, framealpha=0.9)
    setup_axis_style(ax)
    ax.set_box_aspect(1)

    # 右パネルのスタイル
    ax = axes[1]
    ax.set_xlabel(r'$x/L$', fontsize=11)
    ax.set_ylabel(r'$v/U_{lid}$', fontsize=11)
    ax.set_title(r'$v$-velocity along horizontal centerline', fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.5, 0.4)
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)
    setup_axis_style(ax)
    ax.set_box_aspect(1)

    # 全体タイトル
    fig.suptitle(f'Lid-Driven Cavity Flow Validation (Re = {Re})',
                 fontsize=13, y=0.98)

    plt.tight_layout()

    # 保存（SVG, PDF, PNG）
    if save_file is None:
        save_file = f'validation_Re{Re}_comparison.svg'

    base_path = os.path.splitext(save_file)[0]
    for ext in ['.svg', '.pdf', '.png']:
        out_file = base_path + ext
        dpi = 300 if ext == '.png' else None
        plt.savefig(out_file, bbox_inches='tight', dpi=dpi)
        print(f"Saved: {out_file}")


def compute_error(profiles: dict, ghia_u: dict, ghia_v: dict, Re: int):
    """Ghiaデータとの誤差を計算"""

    # 補間して比較
    from scipy import interpolate

    # u速度の誤差
    f_u = interpolate.interp1d(profiles['y_vertical'], profiles['u_vertical'],
                                kind='linear', fill_value='extrapolate')
    u_interp = f_u(ghia_u['y'])
    u_error = np.sqrt(np.mean((u_interp - np.array(ghia_u['u']))**2))
    u_max_error = np.max(np.abs(u_interp - np.array(ghia_u['u'])))

    # v速度の誤差
    f_v = interpolate.interp1d(profiles['x_horizontal'], profiles['v_horizontal'],
                                kind='linear', fill_value='extrapolate')
    v_interp = f_v(ghia_v['x'])
    v_error = np.sqrt(np.mean((v_interp - np.array(ghia_v['v']))**2))
    v_max_error = np.max(np.abs(v_interp - np.array(ghia_v['v'])))

    print(f"\n=== Validation Results (Re = {Re}) ===")
    print(f"u-velocity (vertical centerline):")
    print(f"  RMS error:  {u_error:.6f}")
    print(f"  Max error:  {u_max_error:.6f}")
    print(f"v-velocity (horizontal centerline):")
    print(f"  RMS error:  {v_error:.6f}")
    print(f"  Max error:  {v_max_error:.6f}")
    print("=" * 40)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Validate cavity flow simulation')
    parser.add_argument('output_dirs', nargs='+', help='Output directories containing CSV files')
    parser.add_argument('--labels', nargs='+', default=None,
                        help='Labels for each case (e.g., "32x32 Projection" "64x64 SIMPLE")')
    parser.add_argument('--Re', type=int, default=100, help='Reynolds number')
    parser.add_argument('--U_lid', type=float, default=0.01, help='Lid velocity [m/s]')
    parser.add_argument('--save', type=str, default=None, help='Save figure to file')

    args = parser.parse_args()

    if len(args.output_dirs) == 1:
        # 単一ケース
        plot_validation(args.output_dirs[0], args.Re, args.U_lid, args.save)
    else:
        # 複数ケース比較
        labels = args.labels if args.labels else [os.path.basename(d) for d in args.output_dirs]
        if len(labels) != len(args.output_dirs):
            print("Error: Number of labels must match number of output directories")
            exit(1)
        plot_multi_case_validation(args.output_dirs, labels, args.Re, args.U_lid, args.save)
