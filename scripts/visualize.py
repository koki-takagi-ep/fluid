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


def format_colorbar_ticks(cbar, n_ticks=5):
    """カラーバーのティックをキリの良い数字にフォーマット（有効数字2桁）"""
    import matplotlib.ticker as ticker

    vmin, vmax = cbar.vmin, cbar.vmax

    # 適切な間隔を計算（キリの良い数字に）
    range_val = vmax - vmin
    if range_val == 0:
        return

    # 桁数を計算
    order = np.floor(np.log10(range_val))
    step_base = 10 ** order

    # キリの良い間隔を選択（1, 2, 2.5, 5の倍数）
    nice_steps = [1, 2, 2.5, 5, 10]
    for step_mult in nice_steps:
        step = step_base * step_mult / 10
        n_steps = range_val / step
        if n_steps <= n_ticks + 1:
            break

    # ティック位置を計算
    tick_min = np.ceil(vmin / step) * step
    tick_max = np.floor(vmax / step) * step
    ticks = np.arange(tick_min, tick_max + step / 2, step)

    # 有効数字2桁でフォーマット
    cbar.set_ticks(ticks)

    # フォーマット関数
    def format_tick(x, pos):
        if x == 0:
            return '0'
        abs_x = abs(x)
        if abs_x >= 1000:
            return f'{x:.2g}'
        elif abs_x >= 1:
            return f'{x:.3g}'
        else:
            return f'{x:.2g}'

    cbar.ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_tick))
    cbar.ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_tick))


def setup_colorbar_style(cbar, label='', horizontal=False):
    """カラーバーのスタイルを設定"""
    cbar.set_label(label, fontsize=10, fontweight='bold')
    cbar.ax.tick_params(direction='out', length=4, width=1, labelsize=9)
    format_colorbar_ticks(cbar)


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
                        length_unit='mm', vel_unit='mm/s', horizontal_cbar=False):
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

    # カラーバー（横向き対応）
    if horizontal_cbar:
        cbar = plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.15, aspect=30)
    else:
        cbar = plt.colorbar(cf, ax=ax)
    setup_colorbar_style(cbar, label=f'Velocity $|\\vec{{u}}|$ ({vel_unit})')

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


def plot_pressure_field(field: dict, ax=None, title=None, length_unit='mm', horizontal_cbar=False):
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

    # カラーバー（横向き対応）
    if horizontal_cbar:
        cbar = plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.15, aspect=30)
    else:
        cbar = plt.colorbar(cf, ax=ax)
    setup_colorbar_style(cbar, label='Pressure $p$ (Pa)')

    # スタイル設定
    title_text = title if title else f'Pressure at t = {field["time"]:.4f} s'
    setup_axis_style(ax, xlabel=f'$x$ ({length_unit})', ylabel=f'$y$ ({length_unit})', title=title_text)
    ax.set_aspect('equal')

    return ax


def plot_streamlines(field: dict, ax=None, title=None, length_unit='mm', vel_unit='mm/s', horizontal_cbar=False):
    """流線をプロット（流れの種類に応じて最適な方法を選択）"""
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

    # カラーバー（横向き対応）
    if horizontal_cbar:
        cbar = plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.15, aspect=30)
    else:
        cbar = plt.colorbar(cf, ax=ax)
    setup_colorbar_style(cbar, label=f'Velocity $|\\vec{{u}}|$ ({vel_unit})')

    # アスペクト比で流れの種類を判定
    x_range = x_1d[-1] - x_1d[0]
    y_range = y_1d[-1] - y_1d[0]
    is_channel_flow = (x_range / y_range) > 2  # アスペクト比が2以上ならチャネルフロー

    if is_channel_flow:
        # チャネルフロー: 左端から右へ流れる流線を積分で計算
        u_interp = RegularGridInterpolator((x_1d, y_1d), u, bounds_error=False, fill_value=None)
        v_interp = RegularGridInterpolator((x_1d, y_1d), v, bounds_error=False, fill_value=None)

        def velocity(t, pos):
            x, y = pos
            ux = u_interp([x, y])[0]
            vy = v_interp([x, y])[0]
            if np.isnan(ux) or np.isnan(vy):
                return [0.0, 0.0]
            return [ux, vy]

        n_streams = 15
        y_margin = y_range * 0.03
        y_starts = np.linspace(y_1d[0] + y_margin, y_1d[-1] - y_margin, n_streams)
        x_start = x_1d[0]

        mean_u = max(np.mean(np.abs(u)), 1e-6)
        t_max = (x_range / mean_u) * 5.0
        t_span = [0, t_max]
        n_points = 5000
        t_eval = np.linspace(0, t_max, n_points)

        x_min_bound, x_max_bound = x_1d[0], x_1d[-1]
        y_min_bound, y_max_bound = y_1d[0], y_1d[-1]

        for y0 in y_starts:
            try:
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
                               max_step=x_range/100)

                if len(sol.y[0]) > 1:
                    mask = ((sol.y[0] >= x_min_bound) & (sol.y[0] <= x_max_bound) &
                            (sol.y[1] >= y_min_bound) & (sol.y[1] <= y_max_bound))
                    if np.sum(mask) > 1:
                        x_line = sol.y[0][mask]
                        y_line = sol.y[1][mask]
                        if x_line[0] > x_min_bound:
                            x_line = np.insert(x_line, 0, x_min_bound)
                            y_line = np.insert(y_line, 0, y_line[0])
                        ax.plot(x_line, y_line, 'k-', linewidth=0.3, alpha=1.0)
            except Exception:
                pass
    else:
        # キャビティフロー等: matplotlibのstreamplotを使用
        # 速度場をより細かい格子に補間して滑らかな流線を描画
        n_fine = 100
        x_fine = np.linspace(x_1d[0], x_1d[-1], n_fine)
        y_fine = np.linspace(y_1d[0], y_1d[-1], n_fine)
        X_fine, Y_fine = np.meshgrid(x_fine, y_fine)

        u_interp = RegularGridInterpolator((x_1d, y_1d), u, bounds_error=False, fill_value=0)
        v_interp = RegularGridInterpolator((x_1d, y_1d), v, bounds_error=False, fill_value=0)

        U_fine = np.zeros((n_fine, n_fine))
        V_fine = np.zeros((n_fine, n_fine))
        for i in range(n_fine):
            for j in range(n_fine):
                U_fine[j, i] = u_interp([x_fine[i], y_fine[j]])[0]
                V_fine[j, i] = v_interp([x_fine[i], y_fine[j]])[0]

        # 流線を描画（密度を高めに、極細の黒線）
        ax.streamplot(X_fine, Y_fine, U_fine, V_fine,
                      color='black', linewidth=0.3, density=2.0, arrowsize=0.5)

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


def plot_final_state(output_dir: str, save_file: str = None, dpi: int = 600):
    """最終状態を縦一列で表示（4パネル）"""
    metadata = load_metadata(output_dir)
    nx, ny = int(metadata['nx']), int(metadata['ny'])

    files = get_field_files(output_dir)
    if len(files) == 0:
        print("No field files found!")
        return

    # デフォルトの出力先はfiguresディレクトリ（SVG形式）
    if save_file is None:
        save_file = os.path.join(get_figures_dir(output_dir), 'result.svg')

    field = load_field(files[-1], nx, ny)

    # チャネルフローかどうかを判定（アスペクト比で判断）
    x_range = field['x'].max() - field['x'].min()
    y_range = field['y'].max() - field['y'].min()
    aspect_ratio = x_range / y_range
    is_channel_flow = aspect_ratio > 3  # アスペクト比が3以上ならチャネルフロー

    # mm単位での最大値
    L_scale = 1000.0  # mm
    x_max_mm = field['x'].max() * L_scale
    y_max_mm = field['y'].max() * L_scale

    # 縦一列レイアウト（4行1列）
    # カラーマップ図の高さはアスペクト比に応じて調整
    # centerline velocityは正方形
    if is_channel_flow:
        # チャネルフロー: 横長の図
        panel_width = 10
        panel_height = panel_width / aspect_ratio + 1.5  # カラーバー分の余白
        square_size = 5  # 正方形プロットのサイズ
        fig_height = panel_height * 3 + square_size + 2  # 余白
        fig = plt.figure(figsize=(panel_width, fig_height))

        # GridSpecで柔軟なレイアウト
        from matplotlib.gridspec import GridSpec
        gs = GridSpec(4, 1, figure=fig, height_ratios=[panel_height, panel_height, panel_height, square_size],
                      hspace=0.4)
    else:
        # キャビティフロー: 正方形に近い図
        panel_size = 6
        fig_height = panel_size * 4 + 3  # 余白
        fig = plt.figure(figsize=(panel_size, fig_height))

        from matplotlib.gridspec import GridSpec
        gs = GridSpec(4, 1, figure=fig, height_ratios=[1, 1, 1, 1], hspace=0.35)

    # 速度場（ベクトル＋カラー）
    ax1 = fig.add_subplot(gs[0])
    plot_velocity_field(field, ax=ax1, show_streamlines=False,
                        title=r'Velocity Field $|\vec{u}|$',
                        horizontal_cbar=True)
    if is_channel_flow:
        setup_channel_flow_axes(ax1, x_max_mm, y_max_mm)

    # 流線
    ax2 = fig.add_subplot(gs[1])
    plot_streamlines(field, ax=ax2, title=r'Streamlines $|\vec{u}|$',
                     horizontal_cbar=True)
    if is_channel_flow:
        setup_channel_flow_axes(ax2, x_max_mm, y_max_mm)

    # 圧力場
    ax3 = fig.add_subplot(gs[2])
    plot_pressure_field(field, ax=ax3, title=r'Pressure Field $p$',
                        horizontal_cbar=True)
    if is_channel_flow:
        setup_channel_flow_axes(ax3, x_max_mm, y_max_mm)

    # 中心線速度プロファイル（正方形）
    ax4 = fig.add_subplot(gs[3])
    plot_centerline_velocity(field, ax=ax4)
    ax4.set_box_aspect(1)  # 正方形のaxes box

    # レイアウト調整（GridSpecを使用しているためsubplots_adjustを使用）
    fig.subplots_adjust(left=0.12, right=0.95, top=0.97, bottom=0.05)

    # SVG, PDF, PNG の3形式で保存
    base_path = os.path.splitext(save_file)[0]

    # SVG
    svg_file = base_path + '.svg'
    plt.savefig(svg_file, bbox_inches='tight')
    print(f"Figure saved to {svg_file}")

    # PDF
    pdf_file = base_path + '.pdf'
    plt.savefig(pdf_file, bbox_inches='tight')
    print(f"Figure saved to {pdf_file}")

    # PNG
    png_file = base_path + '.png'
    plt.savefig(png_file, bbox_inches='tight', dpi=dpi)
    print(f"Figure saved to {png_file} (dpi={dpi})")


def main():
    parser = argparse.ArgumentParser(description='Visualize fluid simulation results')
    parser.add_argument('output_dir', help='Output directory containing CSV files')
    parser.add_argument('--plot-final', action='store_true', help='Plot final state')
    parser.add_argument('--animation', action='store_true', help='Create animation')
    parser.add_argument('--streamlines', action='store_true', help='Plot streamlines only')
    parser.add_argument('--save', type=str, default=None, help='Save figure to file')
    parser.add_argument('--fps', type=int, default=10, help='Animation FPS')
    parser.add_argument('--dpi', type=int, default=600, help='Output DPI for PNG images (default: 600)')

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
            save_file = args.save if args.save else os.path.join(get_figures_dir(args.output_dir), 'streamlines.svg')
            plt.savefig(save_file, bbox_inches='tight')
            print(f"Figure saved to {save_file}")
    else:
        # デフォルトは最終状態の4パネル表示
        plot_final_state(args.output_dir, args.save, dpi=args.dpi)


if __name__ == '__main__':
    main()
