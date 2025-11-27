/**
 * @file Grid.hpp
 * @brief スタガード格子（MAC格子）の定義
 *
 * MAC法（Marker-and-Cell）で使用するスタガード格子を実装。
 * 非圧縮性流れの数値解法において、圧力-速度の結合問題を
 * 安定的に解くための格子配置。
 *
 * @see Harlow, F. H., & Welch, J. E. (1965). Physics of Fluids, 8(12), 2182-2189.
 */

#pragma once

#include <vector>
#include <cmath>

namespace fluid {

/**
 * @class Grid
 * @brief 2次元スタガード格子クラス
 *
 * 変数の配置:
 * - 圧力 p[i][j]: セル(i,j)の中心
 * - 速度 u[i][j]: セル(i-1,j)と(i,j)の境界面（x方向）
 * - 速度 v[i][j]: セル(i,j-1)と(i,j)の境界面（y方向）
 *
 * インデックスの範囲（ゴーストセルを含む）:
 * - p: [0, nx+1] x [0, ny+1]
 * - u: [0, nx] x [0, ny+1]
 * - v: [0, nx+1] x [0, ny]
 *
 * @code
 *     v[i][j+1]
 *        |
 *  u[i][j] -- p[i][j] -- u[i+1][j]
 *        |
 *     v[i][j]
 * @endcode
 */
class Grid {
public:
    //==========================================================================
    // 格子パラメータ
    //==========================================================================

    int nx;       ///< x方向のセル数
    int ny;       ///< y方向のセル数
    double lx;    ///< x方向の領域サイズ [m]
    double ly;    ///< y方向の領域サイズ [m]
    double dx;    ///< x方向の格子幅 [m]
    double dy;    ///< y方向の格子幅 [m]

    //==========================================================================
    // 変数場（ゴーストセルを含む）
    //==========================================================================

    std::vector<std::vector<double>> p;       ///< 圧力場 [Pa]: (nx+2) x (ny+2)
    std::vector<std::vector<double>> u;       ///< x方向速度 [m/s]: (nx+1) x (ny+2)
    std::vector<std::vector<double>> v;       ///< y方向速度 [m/s]: (nx+2) x (ny+1)

    std::vector<std::vector<double>> u_star;  ///< 中間速度場 u* [m/s]
    std::vector<std::vector<double>> v_star;  ///< 中間速度場 v* [m/s]

    //==========================================================================
    // コンストラクタ
    //==========================================================================

    /**
     * @brief 格子を初期化
     * @param nx x方向のセル数
     * @param ny y方向のセル数
     * @param lx x方向の領域サイズ [m]
     * @param ly y方向の領域サイズ [m]
     */
    Grid(int nx, int ny, double lx, double ly);

    //==========================================================================
    // 発散計算
    //==========================================================================

    /**
     * @brief セル(i,j)での速度場の発散を計算
     *
     * ∇・u = ∂u/∂x + ∂v/∂y ≈ (u[i+1][j] - u[i][j])/dx + (v[i][j+1] - v[i][j])/dy
     *
     * @param i x方向のセルインデックス [1, nx]
     * @param j y方向のセルインデックス [1, ny]
     * @return 発散値 [1/s]
     */
    double computeDivergence(int i, int j) const;

    /**
     * @brief 領域全体での最大発散を計算
     * @return 最大発散値 [1/s]
     */
    double maxDivergence() const;

    //==========================================================================
    // 座標取得
    //==========================================================================

    /**
     * @brief セル(i,j)の中心のx座標
     * @param i セルインデックス
     * @return x座標 [m]
     */
    double cellCenterX(int i) const { return (i + 0.5) * dx; }

    /**
     * @brief セル(i,j)の中心のy座標
     * @param j セルインデックス
     * @return y座標 [m]
     */
    double cellCenterY(int j) const { return (j + 0.5) * dy; }

    /// u速度点のx座標
    double uPositionX(int i) const { return i * dx; }

    /// u速度点のy座標
    double uPositionY(int j) const { return (j + 0.5) * dy; }

    /// v速度点のx座標
    double vPositionX(int i) const { return (i + 0.5) * dx; }

    /// v速度点のy座標
    double vPositionY(int j) const { return j * dy; }

    //==========================================================================
    // 補間
    //==========================================================================

    /**
     * @brief セル中心でのu速度を補間
     *
     * u_center = 0.5 * (u[i][j+1] + u[i+1][j+1])
     *
     * @param i セルインデックス
     * @param j セルインデックス
     * @return セル中心でのu速度 [m/s]
     */
    double uAtCellCenter(int i, int j) const;

    /**
     * @brief セル中心でのv速度を補間
     *
     * v_center = 0.5 * (v[i+1][j] + v[i+1][j+1])
     *
     * @param i セルインデックス
     * @param j セルインデックス
     * @return セル中心でのv速度 [m/s]
     */
    double vAtCellCenter(int i, int j) const;

    /**
     * @brief セル中心での速度の大きさを計算
     * @param i セルインデックス
     * @param j セルインデックス
     * @return 速度の大きさ |u| [m/s]
     */
    double velocityMagnitude(int i, int j) const;

private:
    /// 2次元配列を初期化
    static std::vector<std::vector<double>> make2DArray(int rows, int cols);
};

} // namespace fluid
