/**
 * @file Solver.hpp
 * @brief Projection法による非圧縮性Navier-Stokes方程式ソルバー
 *
 * Projection法（Chorin's Method）による時間積分を実装。
 *
 * @see Chorin, A. J. (1968). Math. Comp., 22(104), 745-762.
 */

#pragma once

#include "SolverBase.hpp"
#include "PressureSolver.hpp"

namespace fluid {

/**
 * @class Solver
 * @brief Projection法による非圧縮性Navier-Stokes方程式ソルバー
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
 */
class Solver : public SolverBase {
public:
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
    int step(Grid& grid, const BoundaryCondition& bc) override;
};

} // namespace fluid
