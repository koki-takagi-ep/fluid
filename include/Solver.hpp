/**
 * @file Solver.hpp
 * @brief 非圧縮性Navier-Stokes方程式の時間発展ソルバー
 *
 * Projection法（Chorin's Method）による時間積分を実装。
 *
 * @see Chorin, A. J. (1968). Math. Comp., 22(104), 745-762.
 */

#pragma once

#include "Grid.hpp"
#include "BoundaryCondition.hpp"
#include "PressureSolver.hpp"
#include <string>
#include <functional>

namespace fluid {

/**
 * @class Solver
 * @brief 非圧縮性Navier-Stokes方程式ソルバー
 *
 * Projection法（分離速度場法）による時間積分:
 *
 * 1. **予測ステップ**: 圧力項を除いた運動量方程式
 *    @f[
 *    \mathbf{u}^* = \mathbf{u}^n + \Delta t \left[ -(\mathbf{u}^n \cdot \nabla)\mathbf{u}^n + \nu \nabla^2 \mathbf{u}^n \right]
 *    @f]
 *
 * 2. **圧力Poisson方程式**: 連続の式を満たす圧力場
 *    @f[
 *    \nabla^2 p^{n+1} = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*
 *    @f]
 *
 * 3. **修正ステップ**: 速度場の更新
 *    @f[
 *    \mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla p^{n+1}
 *    @f]
 *
 * 空間離散化:
 * - 移流項: 1次精度風上差分
 * - 拡散項: 2次精度中心差分
 */
class Solver {
public:
    //==========================================================================
    // 物理パラメータ（SI単位系）
    //==========================================================================

    double rho;     ///< 密度 [kg/m³]
    double nu;      ///< 動粘性係数 [m²/s]

    //==========================================================================
    // 数値パラメータ
    //==========================================================================

    double dt;          ///< 時間刻み [s]
    double cfl;         ///< CFL数（自動時間刻み用、典型値: 0.5）
    bool autoTimeStep;  ///< 適応的時間刻みを使用するか

    //==========================================================================
    // ソルバーコンポーネント
    //==========================================================================

    PressureSolver pressureSolver;  ///< 圧力Poisson方程式ソルバー

    //==========================================================================
    // コンストラクタ
    //==========================================================================

    /**
     * @brief ソルバーを初期化
     * @param rho 密度 [kg/m³]（デフォルト: 1000 = 水）
     * @param nu 動粘性係数 [m²/s]（デフォルト: 1e-6 = 水）
     * @param dt 時間刻み [s]
     */
    Solver(double rho = 1000.0, double nu = 1.0e-6, double dt = 0.001);

    //==========================================================================
    // 時間発展
    //==========================================================================

    /**
     * @brief 1タイムステップ進める
     * @param grid 格子
     * @param bc 境界条件
     * @return 圧力ソルバーの反復回数（-1: 非収束）
     */
    int step(Grid& grid, const BoundaryCondition& bc);

    /**
     * @brief シミュレーションを実行
     * @param grid 格子
     * @param bc 境界条件
     * @param endTime 終了時刻 [s]
     * @param outputInterval 出力間隔（ステップ数）
     * @param outputCallback 出力時に呼ばれるコールバック関数
     */
    void run(Grid& grid, const BoundaryCondition& bc, double endTime,
             int outputInterval = 100,
             std::function<void(int, double, const Grid&)> outputCallback = nullptr);

    //==========================================================================
    // アクセサ
    //==========================================================================

    /// 現在時刻 [s]
    double getTime() const { return time; }

    /// 現在のステップ数
    int getStepCount() const { return stepCount; }

private:
    double time = 0.0;    ///< 現在時刻 [s]
    int stepCount = 0;    ///< ステップカウンタ

    //==========================================================================
    // 内部メソッド
    //==========================================================================

    /**
     * @brief 適応時間刻みを計算
     *
     * CFL条件と粘性条件を満たす時間刻みを返す:
     * @f[
     * \Delta t \leq \min\left( \frac{\text{CFL} \cdot \Delta x}{|u|_{\max}}, \frac{\Delta x^2}{4\nu} \right)
     * @f]
     */
    double computeTimeStep(const Grid& grid) const;

    /// 中間速度場 u*, v* を計算（予測ステップ）
    void computeIntermediateVelocity(Grid& grid, const BoundaryCondition& bc);

    /// 速度場を補正（修正ステップ）
    void correctVelocity(Grid& grid);

    /**
     * @brief u速度の移流項を計算（1次風上差分）
     * @f[
     * (\mathbf{u} \cdot \nabla) u \approx u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y}
     * @f]
     */
    double convectionU(const Grid& grid, int i, int j) const;

    /// v速度の移流項を計算
    double convectionV(const Grid& grid, int i, int j) const;

    /**
     * @brief u速度の拡散項を計算（2次中心差分）
     * @f[
     * \nu \nabla^2 u = \nu \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right)
     * @f]
     */
    double diffusionU(const Grid& grid, int i, int j) const;

    /// v速度の拡散項を計算
    double diffusionV(const Grid& grid, int i, int j) const;
};

} // namespace fluid
