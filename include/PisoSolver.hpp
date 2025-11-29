/**
 * @file PisoSolver.hpp
 * @brief PISO法による非圧縮性Navier-Stokes方程式ソルバー
 *
 * Pressure-Implicit with Splitting of Operators (PISO)
 *
 * @see Issa, R. I. (1986). J. Comput. Phys., 62, 40-65.
 */

#pragma once

#include "SolverBase.hpp"
#include <vector>

namespace fluid {

/**
 * @class PisoSolver
 * @brief PISO法による非圧縮性Navier-Stokes方程式ソルバー
 *
 * PISO法はSIMPLE法を拡張した非反復式の手法で、特に非定常問題に適している。
 * 1タイムステップ内で複数回の圧力補正を行うことで、外部反復なしに
 * 高精度な解を得る。
 *
 * PISO法のアルゴリズム:
 *
 * 1. **Predictor（予測ステップ）**:
 *    圧力項を除いた運動量方程式を解き、中間速度 u* を計算:
 *    @f[
 *    \frac{u^* - u^n}{\Delta t} = -(\mathbf{u}^n \cdot \nabla)\mathbf{u}^n + \nu \nabla^2 \mathbf{u}^n
 *    @f]
 *
 * 2. **First Corrector（第1補正ステップ）**:
 *    第1圧力補正方程式を解く:
 *    @f[
 *    \nabla^2 p' = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*
 *    @f]
 *    速度と圧力を補正:
 *    @f[
 *    \mathbf{u}^{**} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla p', \quad p^* = p^n + p'
 *    @f]
 *
 * 3. **Second Corrector（第2補正ステップ）**:
 *    第2圧力補正方程式を解く:
 *    @f[
 *    \nabla^2 p'' = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^{**}
 *    @f]
 *    速度と圧力を補正:
 *    @f[
 *    \mathbf{u}^{n+1} = \mathbf{u}^{**} - \frac{\Delta t}{\rho} \nabla p'', \quad p^{n+1} = p^* + p''
 *    @f]
 *
 * PISO法 vs SIMPLE法:
 * - SIMPLE: 外部反復が必要、緩和係数を使用
 * - PISO: 非反復式（1タイムステップで1回のみ実行）、緩和係数不要
 * - PISO: 非定常問題に適している
 * - PISO: 追加の補正ステップにより精度向上
 */
class PisoSolver : public SolverBase {
public:
    //==========================================================================
    // PISO法固有のパラメータ
    //==========================================================================

    double omega;           ///< SOR緩和係数
    int maxPressureIter;    ///< 圧力補正の最大反復回数
    double pressureTol;     ///< 圧力収束許容値

    int nCorrectors;        ///< 補正ステップ数（デフォルト: 2）

    //==========================================================================
    // コンストラクタ
    //==========================================================================

    /**
     * @brief ソルバーを初期化
     * @param rho 密度 [kg/m³]（デフォルト: 1000 = 水）
     * @param nu 動粘性係数 [m²/s]（デフォルト: 1e-6 = 水）
     * @param dt 時間刻み [s]
     * @param nCorrectors 補正ステップ数（デフォルト: 2）
     */
    PisoSolver(double rho = 1000.0, double nu = 1.0e-6, double dt = 0.001, int nCorrectors = 2);

    //==========================================================================
    // 時間発展
    //==========================================================================

    /**
     * @brief 1タイムステップ進める（PISO法）
     * @param grid 格子
     * @param bc 境界条件
     * @return 圧力補正ステップで使用した総反復回数（-1: 非収束）
     */
    int step(Grid& grid, const BoundaryCondition& bc) override;

    //==========================================================================
    // アクセサ
    //==========================================================================

    /// 最後のステップでの圧力補正反復回数
    int getLastPressureIterations() const { return lastPressureIter_; }

private:
    int lastPressureIter_ = 0;

    /// 圧力補正場
    std::vector<std::vector<double>> p_prime_;

    /// 第2中間速度場（second correctorで使用）
    std::vector<std::vector<double>> u_double_star_;
    std::vector<std::vector<double>> v_double_star_;

    //==========================================================================
    // PISO法固有のメソッド
    //==========================================================================

    /**
     * @brief 圧力補正方程式を解く
     * @param grid 格子
     * @param bc 境界条件
     * @param u_field 速度場（u_starまたはu_double_star）
     * @param v_field 速度場（v_starまたはv_double_star）
     * @return 収束に要した反復回数（-1: 非収束）
     */
    int solvePressureCorrection(Grid& grid, const BoundaryCondition& bc,
                                 const std::vector<std::vector<double>>& u_field,
                                 const std::vector<std::vector<double>>& v_field);

    /// 2次元配列を初期化
    void ensureArraySize(std::vector<std::vector<double>>& arr, int rows, int cols);

    /**
     * @brief 圧力補正p'による速度補正（PISO法専用）
     * @param grid 格子
     * @param u_field 入力速度場（補正前）
     * @param v_field 入力速度場（補正前）
     *
     * u^{n+1} = u_field - (dt/ρ) * ∇p'
     */
    void correctVelocityWithPressureCorrection(Grid& grid,
                                                const std::vector<std::vector<double>>& u_field,
                                                const std::vector<std::vector<double>>& v_field);
};

} // namespace fluid
