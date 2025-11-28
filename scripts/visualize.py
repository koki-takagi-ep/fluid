#!/usr/bin/env python3
"""
非圧縮流体シミュレーション結果の可視化スクリプト

使用方法:
    python visualize.py <output_dir> [options]

例:
    python visualize.py output_cavity --plot-final
    python visualize.py output_cavity --animation
    python visualize.py output_cavity --streamlines
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import matplotlib.cm as cm


def setup_axis_style(ax, xlabel='', ylabel='', title=''):
    """軸のスタイルを設定（論文調）"""
    # 軸ラベル
    ax.set_xlabel(xlabel, fontsize=10, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=10, fontweight='bold')
    if title:
        ax.set_title(title, fontsize=11, fontweight='bold', pad=12)  # タイトルとグラフの間隔

    # 目盛りを外向きに設定
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.tick_params(which='major', direction='out', length=6, width=1, labelsize=9)
    ax.tick_params(which='minor', direction='out', length=3, width=0.5)

    # マイナー目盛りを追加
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))


def setup_channel_flow_axes(ax, x_max, y_max):
    """チャネルフロー用の軸設定

    x軸: 0から(x_maxを含む5の倍数)まで、5mm間隔のメジャーティック、1mm間隔のマイナーティック
    y軸: 0, 1, 2, ... のメジャーティックのみ（マイナーティックなし）
    """
    # x軸の設定
    x_tick_max = int(np.ceil(x_max / 5) * 5)  # 5の倍数に切り上げ
    ax.set_xlim(0, x_tick_max)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(MultipleLocator(1))

    # y軸の設定
    y_tick_max = int(np.ceil(y_max))
    ax.set_ylim(0, y_tick_max)
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(plt.NullLocator())  # マイナーティックなし


def setup_colorbar_style(cbar, label=''):
    """カラーバーのスタイルを設定"""
    cbar.set_label(label, fontsize=10, fontweight='bold')
    cbar.ax.tick_params(direction='out', length=4, width=1, labelsize=9)


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


def load_metadata(output_dir: str) -> dict:
    """メタデータを読み込む"""
    data_dir = get_data_dir(output_dir)
    metadata_file = os.path.join(data_dir, "metadata.csv")
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

    df = pd.read_csv(metadata_file)
    return dict(zip(df['parameter'], df['value']))


def load_field(filepath: str, nx: int, ny: int) -> dict:
    """フィールドデータを読み込む"""
    with open(filepath, 'r') as f:
        first_line = f.readline()
        time = float(first_line.split('=')[1])

    df = pd.read_csv(filepath, comment='#')

    x = df['x'].values.reshape(nx, ny)
    y = df['y'].values.reshape(nx, ny)
    u = df['u'].values.reshape(nx, ny)
    v = df['v'].values.reshape(nx, ny)
    p = df['p'].values.reshape(nx, ny)
    mag = df['magnitude'].values.reshape(nx, ny)

    return {
        'time': time,
        'x': x, 'y': y,
        'u': u, 'v': v,
        'p': p, 'magnitude': mag
    }


def get_field_files(output_dir: str) -> list:
    """フィールドファイルのリストを取得（時間順にソート）"""
    data_dir = get_data_dir(output_dir)
    pattern = os.path.join(data_dir, "field_*.csv")
    files = sorted(glob.glob(pattern))
    return files


def plot_velocity_field(field: dict, ax=None, title=None, show_streamlines=True,
                        length_unit='mm', vel_unit='mm/s'):
    """速度場をプロット"""
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4.5))

    # 単位変換係数
    if length_unit == 'mm':
        L_scale = 1000.0
    elif length_unit == 'm':
        L_scale = 1.0
    else:
        L_scale = 1.0

    if vel_unit == 'mm/s':
        V_scale = 1000.0
    elif vel_unit == 'm/s':
        V_scale = 1.0
    else:
        V_scale = 1.0

    x, y = field['x'] * L_scale, field['y'] * L_scale
    u, v = field['u'] * V_scale, field['v'] * V_scale
    mag = field['magnitude'] * V_scale

    # 速度の大きさをカラーマップで表示
    cf = ax.contourf(x, y, mag, levels=50, cmap='viridis')
    cbar = plt.colorbar(cf, ax=ax)
    setup_colorbar_style(cbar, label=f'Velocity ({vel_unit})')

    # ベクトル場（細かく配置、小さい黒矢印）
    skip_x = max(1, x.shape[0] // 32)  # x方向: 32本程度
    skip_y = max(1, x.shape[1] // 16)  # y方向: 16本程度
    ax.quiver(x[::skip_x, ::skip_y], y[::skip_x, ::skip_y],
              u[::skip_x, ::skip_y], v[::skip_x, ::skip_y],
              color='black', alpha=0.6, width=0.002, scale=None,
              headwidth=2.5, headlength=3, headaxislength=2.5)

    # スタイル設定
    title_text = title if title else f't = {field["time"]:.4f} s'
    setup_axis_style(ax, xlabel=f'$x$ ({length_unit})', ylabel=f'$y$ ({length_unit})', title=title_text)
    ax.set_aspect('equal')

    return ax


def plot_pressure_field(field: dict, ax=None, title=None, length_unit='mm'):
    """圧力場をプロット"""
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4.5))

    # 単位変換係数
    if length_unit == 'mm':
        L_scale = 1000.0
    else:
        L_scale = 1.0

    x, y = field['x'] * L_scale, field['y'] * L_scale
    p = field['p']

    cf = ax.contourf(x, y, p, levels=50, cmap='RdBu_r')
    cbar = plt.colorbar(cf, ax=ax)
    setup_colorbar_style(cbar, label='Pressure (Pa)')

    # スタイル設定
    title_text = title if title else f'Pressure at t = {field["time"]:.4f} s'
    setup_axis_style(ax, xlabel=f'$x$ ({length_unit})', ylabel=f'$y$ ({length_unit})', title=title_text)
    ax.set_aspect('equal')

    return ax


def plot_streamlines(field: dict, ax=None, title=None, length_unit='mm', vel_unit='mm/s'):
    """流線をプロット（手動積分による正確な流線）"""
    from scipy.interpolate import RegularGridInterpolator
    from scipy.integrate import solve_ivp

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4.5))

    # 単位変換係数
    if length_unit == 'mm':
        L_scale = 1000.0
    else:
        L_scale = 1.0

    if vel_unit == 'mm/s':
        V_scale = 1000.0
    else:
        V_scale = 1.0

    # グリッド情報
    nx, ny = field['x'].shape
    x_1d = np.linspace(field['x'].min(), field['x'].max(), nx) * L_scale
    y_1d = np.linspace(field['y'].min(), field['y'].max(), ny) * L_scale

    u = field['u'] * V_scale
    v = field['v'] * V_scale
    mag = field['magnitude'] * V_scale

    # 背景に速度の大きさをプロット
    X, Y = np.meshgrid(x_1d, y_1d)
    cf = ax.contourf(X, Y, mag.T, levels=50, cmap='viridis')
    cbar = plt.colorbar(cf, ax=ax)
    setup_colorbar_style(cbar, label=f'Velocity ({vel_unit})')

    # 速度場の補間関数を作成（境界外挿を許可）
    u_interp = RegularGridInterpolator((x_1d, y_1d), u, bounds_error=False, fill_value=None)
    v_interp = RegularGridInterpolator((x_1d, y_1d), v, bounds_error=False, fill_value=None)

    def velocity(t, pos):
        x, y = pos
        ux = u_interp([x, y])[0]
        vy = v_interp([x, y])[0]
        # NaNの場合は最近傍の値を使用
        if np.isnan(ux) or np.isnan(vy):
            return [0.0, 0.0]
        return [ux, vy]

    # 流線の開始点（y方向に均等配置）
    n_streams = 15
    y_margin = (y_1d[-1] - y_1d[0]) * 0.03  # マージンを小さくして壁に近い流線も描画
    y_starts = np.linspace(y_1d[0] + y_margin, y_1d[-1] - y_margin, n_streams)
    x_start = x_1d[0]  # 入口の左端から開始（途切れを防ぐ）

    # 流線を積分で計算
    x_length = x_1d[-1] - x_1d[0]
    # 平均速度を基準にして十分な積分時間を確保
    mean_u = max(np.mean(np.abs(u)), 1e-6)
    t_max = (x_length / mean_u) * 5.0  # 十分な余裕を持つ
    t_span = [0, t_max]
    # 細かいサンプリングで滑らかな線を描画
    n_points = 5000
    t_eval = np.linspace(0, t_max, n_points)

    # 領域の境界（少し内側にマージンを取る）
    x_min_bound = x_1d[0]
    x_max_bound = x_1d[-1]
    y_min_bound = y_1d[0]
    y_max_bound = y_1d[-1]

    for y0 in y_starts:
        try:
            # 領域外に出たら停止するイベント関数
            def out_of_bounds(t, pos):
                x, y = pos
                if x < x_min_bound or x > x_max_bound:
                    return -1
                if y < y_min_bound or y > y_max_bound:
                    return -1
                return 1
            out_of_bounds.terminal = True
            out_of_bounds.direction = -1

            sol = solve_ivp(velocity, t_span, [x_start, y0], t_eval=t_eval,
                           method='RK45', dense_output=True, events=out_of_bounds,
                           max_step=x_length/100)  # 細かいステップで精度向上

            # 有効な点をプロット
            if len(sol.y[0]) > 1:
                # 領域内の点のみをマスク
                mask = ((sol.y[0] >= x_min_bound) & (sol.y[0] <= x_max_bound) &
                        (sol.y[1] >= y_min_bound) & (sol.y[1] <= y_max_bound))
                if np.sum(mask) > 1:
                    ax.plot(sol.y[0][mask], sol.y[1][mask], 'k-', linewidth=0.6, alpha=0.8)
        except Exception:
            pass

    # スタイル設定
    title_text = title if title else f'Streamlines at t = {field["time"]:.4f} s'
    setup_axis_style(ax, xlabel=f'$x$ ({length_unit})', ylabel=f'$y$ ({length_unit})', title=title_text)
    ax.set_aspect('equal')

    return ax


def plot_centerline_velocity(field: dict, ax=None, length_unit='mm', vel_unit='mm/s'):
    """中心線での速度プロファイル"""
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 4.5))

    # 単位変換係数
    if length_unit == 'mm':
        L_scale = 1000.0
    else:
        L_scale = 1.0

    if vel_unit == 'mm/s':
        V_scale = 1000.0
    else:
        V_scale = 1.0

    nx, ny = field['x'].shape

    # 垂直中心線（x = 0.5 付近）
    i_mid = nx // 2
    y_line = field['y'][i_mid, :] * L_scale
    u_line = field['u'][i_mid, :] * V_scale

    ax.plot(u_line, y_line, 'b-', linewidth=1.5, label=r'$u$ velocity')

    # スタイル設定
    setup_axis_style(ax, xlabel=f'$u$ ({vel_unit})', ylabel=f'$y$ ({length_unit})',
                     title=f'Centerline velocity (t = {field["time"]:.2f} s)')

    ax.legend(loc='best', fontsize=9, frameon=True, edgecolor='black', fancybox=False)

    return ax


def create_animation(output_dir: str, output_file: str = None, fps: int = 10,
                     length_unit='mm', vel_unit='mm/s'):
    """アニメーションを作成"""
    metadata = load_metadata(output_dir)
    nx, ny = int(metadata['nx']), int(metadata['ny'])

    # デフォルトの出力先はfiguresディレクトリ
    if output_file is None:
        output_file = os.path.join(get_figures_dir(output_dir), 'animation.mp4')

    # 単位変換係数
    if length_unit == 'mm':
        L_scale = 1000.0
    else:
        L_scale = 1.0

    if vel_unit == 'mm/s':
        V_scale = 1000.0
    else:
        V_scale = 1.0

    files = get_field_files(output_dir)
    if len(files) == 0:
        print("No field files found!")
        return

    print(f"Creating animation from {len(files)} frames...")

    # 全フレームの最大速度を計算（カラーバー範囲を固定するため）
    max_mag = 0
    for f in files:
        field = load_field(f, nx, ny)
        max_mag = max(max_mag, np.max(field['magnitude']) * V_scale)

    # 固定のレベルを作成
    levels = np.linspace(0, max_mag, 51)

    # 初期フレーム
    field = load_field(files[0], nx, ny)
    x, y = field['x'] * L_scale, field['y'] * L_scale
    skip = max(1, nx // 16)

    # Figure と Axes を作成
    fig, ax = plt.subplots(figsize=(6, 5.5))

    # 初期プロット（カラーバーを作成するため）
    cf = ax.contourf(x, y, field['magnitude'] * V_scale, levels=levels,
                     cmap='viridis', vmin=0, vmax=max_mag)
    cbar = fig.colorbar(cf, ax=ax)
    setup_colorbar_style(cbar, label=f'Velocity ({vel_unit})')

    ax.quiver(x[::skip, ::skip], y[::skip, ::skip],
              field['u'][::skip, ::skip] * V_scale, field['v'][::skip, ::skip] * V_scale,
              color='white', alpha=0.8)

    setup_axis_style(ax, xlabel=f'$x$ ({length_unit})', ylabel=f'$y$ ({length_unit})',
                     title=f't = {field["time"]:.4f} s')
    ax.set_aspect('equal')

    def update(frame):
        ax.clear()
        field = load_field(files[frame], nx, ny)

        # 固定レベルでcontourf（カラーバーの範囲が一定に保たれる）
        ax.contourf(x, y, field['magnitude'] * V_scale, levels=levels,
                    cmap='viridis', vmin=0, vmax=max_mag)

        ax.quiver(x[::skip, ::skip], y[::skip, ::skip],
                  field['u'][::skip, ::skip] * V_scale, field['v'][::skip, ::skip] * V_scale,
                  color='white', alpha=0.8)

        setup_axis_style(ax, xlabel=f'$x$ ({length_unit})', ylabel=f'$y$ ({length_unit})',
                         title=f't = {field["time"]:.4f} s')
        ax.set_aspect('equal')

        return []

    anim = FuncAnimation(fig, update, frames=len(files), interval=1000//fps, blit=False)

    # 保存
    try:
        writer = FFMpegWriter(fps=fps, metadata={'title': 'Fluid Simulation'})
        anim.save(output_file, writer=writer)
        print(f"Animation saved to {output_file}")
    except Exception as e:
        # FFmpegがない場合はGIF形式で保存
        gif_file = output_file.replace('.mp4', '.gif')
        anim.save(gif_file, writer='pillow', fps=fps)
        print(f"Animation saved to {gif_file} (FFmpeg not available)")

    plt.close()


def plot_final_state(output_dir: str, save_file: str = None):
    """最終状態を4パネルで表示"""
    metadata = load_metadata(output_dir)
    nx, ny = int(metadata['nx']), int(metadata['ny'])

    files = get_field_files(output_dir)
    if len(files) == 0:
        print("No field files found!")
        return

    # デフォルトの出力先はfiguresディレクトリ（PDF形式）
    if save_file is None:
        save_file = os.path.join(get_figures_dir(output_dir), 'result.pdf')

    field = load_field(files[-1], nx, ny)

    # チャネルフローかどうかを判定（アスペクト比で判断）
    x_range = field['x'].max() - field['x'].min()
    y_range = field['y'].max() - field['y'].min()
    is_channel_flow = (x_range / y_range) > 3  # アスペクト比が3以上ならチャネルフロー

    # mm単位での最大値
    L_scale = 1000.0  # mm
    x_max_mm = field['x'].max() * L_scale
    y_max_mm = field['y'].max() * L_scale

    fig, axes = plt.subplots(2, 2, figsize=(10, 9))

    # 速度場（ベクトル＋カラー）
    ax = axes[0, 0]
    plot_velocity_field(field, ax=ax, show_streamlines=False,
                        title='Velocity Field')
    if is_channel_flow:
        setup_channel_flow_axes(ax, x_max_mm, y_max_mm)

    # 流線
    ax = axes[0, 1]
    plot_streamlines(field, ax=ax, title='Streamlines')
    if is_channel_flow:
        setup_channel_flow_axes(ax, x_max_mm, y_max_mm)

    # 圧力場
    ax = axes[1, 0]
    plot_pressure_field(field, ax=ax, title='Pressure Field')
    if is_channel_flow:
        setup_channel_flow_axes(ax, x_max_mm, y_max_mm)

    # 中心線速度プロファイル
    ax = axes[1, 1]
    plot_centerline_velocity(field, ax=ax)

    plt.tight_layout()

    plt.savefig(save_file, bbox_inches='tight')
    print(f"Figure saved to {save_file}")


def main():
    parser = argparse.ArgumentParser(description='Visualize fluid simulation results')
    parser.add_argument('output_dir', help='Output directory containing CSV files')
    parser.add_argument('--plot-final', action='store_true', help='Plot final state')
    parser.add_argument('--animation', action='store_true', help='Create animation')
    parser.add_argument('--streamlines', action='store_true', help='Plot streamlines only')
    parser.add_argument('--save', type=str, default=None, help='Save figure to file')
    parser.add_argument('--fps', type=int, default=10, help='Animation FPS')

    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        print(f"Error: Directory {args.output_dir} not found!")
        return

    if args.animation:
        create_animation(args.output_dir, args.save, args.fps)
    elif args.streamlines:
        metadata = load_metadata(args.output_dir)
        nx, ny = int(metadata['nx']), int(metadata['ny'])
        files = get_field_files(args.output_dir)
        if files:
            field = load_field(files[-1], nx, ny)
            fig, ax = plt.subplots(figsize=(10, 8))
            plot_streamlines(field, ax=ax)
            save_file = args.save if args.save else os.path.join(get_figures_dir(args.output_dir), 'streamlines.pdf')
            plt.savefig(save_file, bbox_inches='tight')
            print(f"Figure saved to {save_file}")
    else:
        # デフォルトは最終状態の4パネル表示
        plot_final_state(args.output_dir, args.save)


if __name__ == '__main__':
    main()
