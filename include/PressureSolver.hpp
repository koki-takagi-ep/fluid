/**
 * @file PressureSolver.hpp
 * @brief 圧力Poisson方程式ソルバー
 *
 * SOR法（Successive Over-Relaxation）による反復解法を実装。
 *
 * @see Young, D. M. (1954). Trans. Amer. Math. Soc., 76(1), 92-111.
 */

#pragma once

#include "Grid.hpp"
#include "BoundaryCondition.hpp"

namespace fluid {

/**
 * @class PressureSolver
 * @brief 圧力Poisson方程式のSORソルバー
 *
 * Projection法の第2ステップで解く圧力Poisson方程式:
 * @f[
 * \nabla^2 p = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*
 * @f]
 *
 * 離散化（5点ステンシル）:
 * @f[
 * \frac{p_{i+1,j} - 2p_{i,j} + p_{i-1,j}}{\Delta x^2} + \frac{p_{i,j+1} - 2p_{i,j} + p_{i,j-1}}{\Delta y^2} = \text{RHS}
 * @f]
 *
 * SOR法による反復:
 * @f[
 * p_{i,j}^{(k+1)} = (1-\omega) p_{i,j}^{(k)} + \frac{\omega}{a_{i,j}} \left( b_{i,j} - \sum_{m \neq i,j} a_{m} p_m \right)
 * @f]
 *
 * 最適緩和係数（正方格子）:
 * @f[
 * \omega_{\text{opt}} = \frac{2}{1 + \sin(\pi / N)}
 * @f]
 */
class PressureSolver {
public:
    //==========================================================================
    // パラメータ
    //==========================================================================

    double omega;       ///< SOR緩和係数（1.0 < ω < 2.0、典型値: 1.8）
    int maxIterations;  ///< 最大反復回数
    double tolerance;   ///< 収束判定の許容残差

    //==========================================================================
    // コンストラクタ
    //==========================================================================

    /**
     * @brief ソルバーを初期化
     * @param omega SOR緩和係数（デフォルト: 1.8）
     * @param maxIter 最大反復回数（デフォルト: 10000）
     * @param tol 収束許容値（発散の閾値、デフォルト: 1e-6）
     */
    PressureSolver(double omega = 1.8, int maxIter = 10000, double tol = 1e-6);

    //==========================================================================
    // ソルバー
    //==========================================================================

    /**
     * @brief 圧力Poisson方程式を解く
     * @param grid 格子（圧力場が更新される）
     * @param bc 境界条件
     * @param dt 時間刻み [s]
     * @param rho 密度 [kg/m³]
     * @return 収束に要した反復回数（-1: 非収束）
     */
    int solve(Grid& grid, const BoundaryCondition& bc, double dt, double rho);

    //==========================================================================
    // アクセサ
    //==========================================================================

    /// 最後の反復での残差を取得
    double getLastResidual() const { return lastResidual; }

private:
    double lastResidual = 0.0;  ///< 最終残差

    /**
     * @brief Poisson方程式の右辺を計算
     *
     * RHS = (ρ/Δt) * ∇・u*
     */
    double computeRHS(const Grid& grid, int i, int j, double dt, double rho) const;
};

} // namespace fluid
