#!/usr/bin/env python3
"""
Channel Flow (Hagen-Poiseuille) の検証スクリプト

平行平板間の完全発達層流（Hagen-Poiseuille流れ）の理論解との比較

理論解:
  u(y) = U_max * (1 - (2y/H - 1)^2)

  ここで:
  - U_max: 中心（y = H/2）での最大速度
  - H: チャネル高さ
  - y: 壁面からの距離 (0 <= y <= H)

圧力勾配との関係:
  dp/dx = -8 * mu * U_max / H^2

平均流速:
  U_mean = (2/3) * U_max
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


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
    files = sorted(glob.glob(os.path.join(data_dir, "field_*.csv")))
    if not files:
        raise FileNotFoundError("No field files found")

    last_file = files[-1]
    df = pd.read_csv(last_file, comment='#')

    x = df['x'].values.reshape(nx, ny)
    y = df['y'].values.reshape(nx, ny)
    u = df['u'].values.reshape(nx, ny)
    v = df['v'].values.reshape(nx, ny)
    p = df['p'].values.reshape(nx, ny)

    return {
        'nx': nx, 'ny': ny, 'lx': lx, 'ly': ly,
        'x': x, 'y': y, 'u': u, 'v': v, 'p': p
    }


def poiseuille_velocity(y: np.ndarray, H: float, U_max: float) -> np.ndarray:
    """Hagen-Poiseuille流れの速度分布（理論解）

    u(y) = U_max * (1 - (2y/H - 1)^2)

    ここで y ∈ [0, H], 最大速度は y = H/2 で発生
    """
    eta = 2.0 * y / H - 1.0  # eta ∈ [-1, 1]
    return U_max * (1.0 - eta**2)


def poiseuille_pressure_gradient(mu: float, U_max: float, H: float) -> float:
    """Hagen-Poiseuille流れの圧力勾配（理論解）

    dp/dx = -8 * mu * U_max / H^2
    """
    return -8.0 * mu * U_max / (H**2)


def setup_axis_style(ax):
    """軸のスタイルを設定（内向き目盛り）"""
    ax.tick_params(axis='both', which='major', direction='in', length=6, width=1,
                   labelsize=10, top=True, right=True)
    ax.tick_params(axis='both', which='minor', direction='in', length=3, width=0.8,
                   top=True, right=True)

    for spine in ax.spines.values():
        spine.set_linewidth(1)


def extract_velocity_profiles(data: dict, x_positions: list = None):
    """指定したx位置での速度プロファイルを抽出

    Args:
        data: シミュレーションデータ
        x_positions: 抽出するx位置のリスト（無次元）[0, 1]
                     Noneの場合は [0.25, 0.5, 0.75, 1.0] を使用

    Returns:
        dict: 各x位置での速度プロファイル
    """
    if x_positions is None:
        x_positions = [0.25, 0.5, 0.75, 1.0]

    nx, ny = data['nx'], data['ny']
    lx, ly = data['lx'], data['ly']

    profiles = {}
    for x_norm in x_positions:
        # x位置に最も近いインデックスを見つける
        x_target = x_norm * lx
        i = int(round(x_norm * (nx - 1)))
        i = min(max(i, 0), nx - 1)

        y_profile = data['y'][i, :] / ly  # 無次元化
        u_profile = data['u'][i, :]

        profiles[x_norm] = {
            'y': y_profile,
            'u': u_profile,
            'x_actual': data['x'][i, 0] / lx
        }

    return profiles


def compute_errors(y_sim: np.ndarray, u_sim: np.ndarray,
                   H: float, U_max: float) -> dict:
    """理論解との誤差を計算

    Args:
        y_sim: シミュレーションのy座標（次元付き）
        u_sim: シミュレーションのu速度
        H: チャネル高さ
        U_max: 最大速度

    Returns:
        dict: 各種誤差指標
    """
    # 理論解
    u_theory = poiseuille_velocity(y_sim * H, H, U_max)

    # 誤差計算
    error = u_sim - u_theory

    # 内部点のみで評価（壁境界を除く）
    interior = (y_sim > 0.05) & (y_sim < 0.95)
    error_interior = error[interior]

    return {
        'rms': np.sqrt(np.mean(error**2)),
        'rms_interior': np.sqrt(np.mean(error_interior**2)),
        'max': np.max(np.abs(error)),
        'max_interior': np.max(np.abs(error_interior)),
        'l2_norm': np.sqrt(np.sum(error**2)) / len(error),
        'relative_rms': np.sqrt(np.mean(error**2)) / U_max
    }


def plot_validation(output_dir: str, U_max: float, H: float = None,
                    save_file: str = None):
    """検証プロットを作成（単一ケース、正方形グラフ）

    Args:
        output_dir: 出力ディレクトリ
        U_max: 最大速度（チャネル中心での速度）[m/s]
        H: チャネル高さ [m]（Noneの場合はメタデータから取得）
        save_file: 保存先ファイルパス
    """
    # データ読み込み
    data = load_simulation_data(output_dir)

    if H is None:
        H = data['ly']

    # 速度プロファイル抽出
    x_positions = [0.1, 0.25, 0.5, 0.75, 0.9]
    profiles = extract_velocity_profiles(data, x_positions)

    # 理論解
    y_theory = np.linspace(0, 1, 100)
    u_theory = poiseuille_velocity(y_theory * H, H, U_max) / U_max  # 無次元化

    # プロット作成（正方形パネル）
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # 左: 速度プロファイルの比較
    ax = axes[0]
    colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(x_positions)))

    # 理論解（黒線）
    ax.plot(u_theory, y_theory, 'k-', linewidth=1, label='Theory (Poiseuille)', zorder=10)

    # 数値解（各x位置）
    for idx, (x_norm, prof) in enumerate(profiles.items()):
        u_norm = prof['u'] / U_max  # 無次元化
        ax.plot(u_norm, prof['y'], 'o-', color=colors[idx],
                markersize=3, linewidth=1,
                label=f'CFD ($x/L_x$ = {prof["x_actual"]:.2f})', zorder=2)

    ax.set_xlabel(r'$u/U_{max}$', fontsize=11)
    ax.set_ylabel(r'$y/H$', fontsize=11)
    ax.set_title('Velocity Profile Comparison', fontsize=12)
    ax.set_xlim(-0.05, 1.1)
    ax.set_ylim(0, 1)
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.legend(loc='center left', fontsize=7, framealpha=0.9)
    setup_axis_style(ax)
    ax.set_box_aspect(1)

    # 右: 出口での詳細比較
    ax = axes[1]

    # 出口に最も近いプロファイル
    x_outlet = max(profiles.keys())
    prof_outlet = profiles[x_outlet]
    y_sim = prof_outlet['y']
    u_sim = prof_outlet['u'] / U_max

    # 理論解との比較
    u_theory_at_sim = poiseuille_velocity(y_sim * H, H, U_max) / U_max

    ax.plot(u_theory_at_sim, y_sim, 'k-', linewidth=1, label='Theory', zorder=10)
    ax.plot(u_sim, y_sim, 'bo', markersize=5,
            label=f'CFD ($x/L_x$ = {prof_outlet["x_actual"]:.2f})', zorder=2)

    # 誤差を計算
    errors = compute_errors(y_sim, prof_outlet['u'], H, U_max)

    ax.set_xlabel(r'$u/U_{max}$', fontsize=11)
    ax.set_ylabel(r'$y/H$', fontsize=11)
    ax.set_title('Outlet Velocity Profile', fontsize=12)
    ax.set_xlim(-0.05, 1.1)
    ax.set_ylim(0, 1)
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.legend(loc='center left', fontsize=9, framealpha=0.9)
    setup_axis_style(ax)
    ax.set_box_aspect(1)

    # 誤差情報をテキストで追加
    error_text = (f"RMS error: {errors['relative_rms']*100:.2f}%\n"
                  f"Max error: {errors['max']/U_max*100:.2f}%")
    ax.text(0.95, 0.05, error_text, transform=ax.transAxes, fontsize=9,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    # 全体タイトル
    fig.suptitle('Hagen-Poiseuille Flow Validation', fontsize=13, y=0.98)

    plt.tight_layout()

    # 保存
    if save_file is None:
        figures_dir = get_figures_dir(output_dir)
        save_file = os.path.join(figures_dir, 'channel_validation.svg')

    base_path = os.path.splitext(save_file)[0]
    for ext in ['.svg', '.pdf', '.png']:
        out_file = base_path + ext
        dpi = 300 if ext == '.png' else None
        plt.savefig(out_file, bbox_inches='tight', dpi=dpi)
        print(f"Saved: {out_file}")

    # 誤差の詳細出力
    print("\n=== Validation Results ===")
    print(f"Outlet profile (x/Lx = {prof_outlet['x_actual']:.2f}):")
    print(f"  RMS error (relative): {errors['relative_rms']*100:.4f}%")
    print(f"  Max error (relative): {errors['max']/U_max*100:.4f}%")
    print(f"  RMS error (interior): {errors['rms_interior']/U_max*100:.4f}%")
    print("=" * 30)

    return errors


def plot_multi_case_validation(output_dirs: list, labels: list,
                               U_max: float, H: float = None,
                               save_file: str = None):
    """複数ケース（異なるセル数・スキーム）の検証プロット（正方形グラフ）

    Args:
        output_dirs: 出力ディレクトリのリスト
        labels: 各ケースのラベル
        U_max: 最大速度 [m/s]
        H: チャネル高さ [m]
        save_file: 保存先ファイルパス
    """
    # 理論解
    y_theory = np.linspace(0, 1, 100)

    # プロット作成（正方形パネル）
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    colors = ['#1f77b4', '#d62728', '#2ca02c', '#ff7f0e', '#9467bd']
    linestyles = ['-', '--', '-.', ':']

    errors_all = {}
    H_case = None

    for idx, (output_dir, label) in enumerate(zip(output_dirs, labels)):
        try:
            data = load_simulation_data(output_dir)
            if H is None:
                H_case = data['ly']
            else:
                H_case = H

            # 出口付近のプロファイル
            profiles = extract_velocity_profiles(data, [0.9])
            prof = profiles[0.9]

            color = colors[idx % len(colors)]
            ls = linestyles[idx % len(linestyles)]

            # 左パネル: 速度プロファイル
            u_norm = prof['u'] / U_max
            axes[0].plot(u_norm, prof['y'], color=color, linestyle=ls,
                        linewidth=1, label=label, zorder=3+idx)

            # 誤差計算
            errors = compute_errors(prof['y'], prof['u'], H_case, U_max)
            errors_all[label] = errors

            # 右パネル: 理論解との差
            u_theory_at_sim = poiseuille_velocity(prof['y'] * H_case, H_case, U_max)
            error_percent = (prof['u'] - u_theory_at_sim) / U_max * 100
            axes[1].plot(error_percent, prof['y'], color=color, linestyle=ls,
                        linewidth=1, label=label, zorder=3+idx)

            print(f"{label}: RMS error = {errors['relative_rms']*100:.4f}%")

        except Exception as e:
            print(f"Error loading {output_dir}: {e}")
            continue

    # 理論解を追加（左パネル）
    if H_case is not None:
        u_theory_norm = poiseuille_velocity(y_theory * H_case, H_case, U_max) / U_max
        axes[0].plot(u_theory_norm, y_theory, 'k-', linewidth=1,
                     label='Theory', zorder=10)

    # 左パネルのスタイル
    ax = axes[0]
    ax.set_xlabel(r'$u/U_{max}$', fontsize=11)
    ax.set_ylabel(r'$y/H$', fontsize=11)
    ax.set_title('Velocity Profile at Outlet', fontsize=12)
    ax.set_xlim(-0.05, 1.1)
    ax.set_ylim(0, 1)
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.legend(loc='center left', fontsize=8, framealpha=0.9)
    setup_axis_style(ax)
    ax.set_box_aspect(1)

    # 右パネルのスタイル
    ax = axes[1]
    ax.axvline(x=0, color='k', linestyle='-', linewidth=0.5)
    ax.set_xlabel('Error [%]', fontsize=11)
    ax.set_ylabel(r'$y/H$', fontsize=11)
    ax.set_title('Deviation from Theory', fontsize=12)
    ax.set_ylim(0, 1)
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.legend(loc='best', fontsize=8, framealpha=0.9)
    setup_axis_style(ax)
    ax.set_box_aspect(1)

    # 全体タイトル
    fig.suptitle('Hagen-Poiseuille Flow Validation', fontsize=13, y=0.98)

    plt.tight_layout()

    # 保存
    if save_file is None:
        save_file = 'channel_validation_comparison.svg'

    base_path = os.path.splitext(save_file)[0]
    for ext in ['.svg', '.pdf', '.png']:
        out_file = base_path + ext
        dpi = 300 if ext == '.png' else None
        plt.savefig(out_file, bbox_inches='tight', dpi=dpi)
        print(f"Saved: {out_file}")

    return errors_all


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Validate channel flow (Hagen-Poiseuille) simulation')
    parser.add_argument('output_dirs', nargs='+',
                        help='Output directories containing CSV files')
    parser.add_argument('--labels', nargs='+', default=None,
                        help='Labels for each case')
    parser.add_argument('--U_max', type=float, default=0.015,
                        help='Maximum velocity (center) [m/s]')
    parser.add_argument('--H', type=float, default=None,
                        help='Channel height [m] (default: from metadata)')
    parser.add_argument('--save', type=str, default=None,
                        help='Save figure to file')

    args = parser.parse_args()

    if len(args.output_dirs) == 1:
        # 単一ケース
        plot_validation(args.output_dirs[0], args.U_max, args.H, args.save)
    else:
        # 複数ケース比較
        labels = args.labels if args.labels else \
                 [os.path.basename(d) for d in args.output_dirs]
        if len(labels) != len(args.output_dirs):
            print("Error: Number of labels must match number of output directories")
            exit(1)
        plot_multi_case_validation(args.output_dirs, labels,
                                   args.U_max, args.H, args.save)
