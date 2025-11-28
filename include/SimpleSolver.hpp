/**
 * @file SimpleSolver.hpp
 * @brief SIMPLE法による非圧縮性Navier-Stokes方程式ソルバー
 *
 * Semi-Implicit Method for Pressure-Linked Equations (SIMPLE)
 *
 * @see Patankar, S. V. & Spalding, D. B. (1972). Int. J. Heat Mass Transfer, 15, 1787-1806.
 */

#pragma once

#include "Grid.hpp"
#include "BoundaryCondition.hpp"
#include <functional>

namespace fluid {

/**
 * @class SimpleSolver
 * @brief SIMPLE法による非圧縮性Navier-Stokes方程式ソルバー
 *
 * SIMPLE法のアルゴリズム:
 *
 * 1. **運動量方程式の離散化**:
 *    @f[
 *    a_P u_P = \sum_{nb} a_{nb} u_{nb} + (p_w - p_e) A + b
 *    @f]
 *
 * 2. **推定速度場の計算** (圧力 p* を使用):
 *    @f[
 *    u^* = \frac{\sum a_{nb} u_{nb}^* + b}{a_P} + \frac{p_w^* - p_e^*}{a_P} A
 *    @f]
 *
 * 3. **圧力補正方程式**:
 *    @f[
 *    \nabla^2 p' = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*
 *    @f]
 *
 * 4. **速度と圧力の補正**:
 *    @f[
 *    p = p^* + \alpha_p p', \quad u = u^* - \frac{\Delta t}{\rho} \frac{\partial p'}{\partial x}
 *    @f]
 *
 * 5. 収束まで繰り返し
 */
class SimpleSolver {
public:
    //==========================================================================
    // 物理パラメータ（SI単位系）
    //==========================================================================

    double rho;     ///< 密度 [kg/m³]
    double nu;      ///< 動粘性係数 [m²/s]

    //==========================================================================
    // 数値パラメータ
    //==========================================================================

    double dt;              ///< 時間刻み [s]
    double cfl;             ///< CFL数（自動時間刻み用）
    bool autoTimeStep;      ///< 適応的時間刻みを使用するか

    double alpha_u;         ///< 速度の緩和係数（典型値: 0.7）
    double alpha_p;         ///< 圧力の緩和係数（典型値: 0.3）

    double omega;           ///< SOR緩和係数
    int maxPressureIter;    ///< 圧力補正の最大反復回数
    double pressureTol;     ///< 圧力収束許容値

    int maxOuterIter;       ///< SIMPLE外部反復の最大回数
    double convergenceTol;  ///< 外部反復の収束判定閾値

    //==========================================================================
    // コンストラクタ
    //==========================================================================

    /**
     * @brief ソルバーを初期化
     * @param rho 密度 [kg/m³]（デフォルト: 1000 = 水）
     * @param nu 動粘性係数 [m²/s]（デフォルト: 1e-6 = 水）
     * @param dt 時間刻み [s]
     */
    SimpleSolver(double rho = 1000.0, double nu = 1.0e-6, double dt = 0.001);

    //==========================================================================
    // 時間発展
    //==========================================================================

    /**
     * @brief 1タイムステップ進める（SIMPLE法）
     * @param grid 格子
     * @param bc 境界条件
     * @return 収束に要した外部反復回数（-1: 非収束）
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

    /// 最後のステップの外部反復回数
    int getLastOuterIterations() const { return lastOuterIter; }

private:
    double time = 0.0;
    int stepCount = 0;
    int lastOuterIter = 0;

    // 圧力補正場
    std::vector<std::vector<double>> p_prime;

    //==========================================================================
    // 内部メソッド
    //==========================================================================

    /// 適応時間刻みを計算
    double computeTimeStep(const Grid& grid) const;

    /// 運動量方程式の係数を計算
    void computeMomentumCoefficients(const Grid& grid, int i, int j,
                                     double& a_P, double& a_E, double& a_W,
                                     double& a_N, double& a_S) const;

    /// 推定速度場 u*, v* を計算
    void computeStarVelocity(Grid& grid, const BoundaryCondition& bc);

    /// 圧力補正方程式を解く
    int solvePressureCorrection(Grid& grid, const BoundaryCondition& bc);

    /// 速度と圧力を補正
    void correctVelocityAndPressure(Grid& grid);

    /// 収束判定（速度残差）
    double computeResidual(const Grid& grid) const;

    /// u速度の移流項
    double convectionU(const Grid& grid, int i, int j) const;

    /// v速度の移流項
    double convectionV(const Grid& grid, int i, int j) const;

    /// u速度の拡散項
    double diffusionU(const Grid& grid, int i, int j) const;

    /// v速度の拡散項
    double diffusionV(const Grid& grid, int i, int j) const;
};

} // namespace fluid
