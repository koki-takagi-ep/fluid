/**
 * @file SolverBase.hpp
 * @brief 非圧縮性Navier-Stokes方程式ソルバーの基底クラス
 *
 * Projection法とSIMPLE法の共通機能を提供する抽象基底クラス。
 */

#pragma once

#include "Grid.hpp"
#include "BoundaryCondition.hpp"
#include <functional>

namespace fluid {

/**
 * @class SolverBase
 * @brief Navier-Stokesソルバーの抽象基底クラス
 *
 * 共通機能:
 * - 時間刻みの計算（CFL条件・粘性条件）
 * - 移流項・拡散項の計算
 * - シミュレーション実行ループ
 */
class SolverBase {
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
    // コンストラクタ・デストラクタ
    //==========================================================================

    /**
     * @brief ソルバーを初期化
     * @param rho 密度 [kg/m³]（デフォルト: 1000 = 水）
     * @param nu 動粘性係数 [m²/s]（デフォルト: 1e-6 = 水）
     * @param dt 時間刻み [s]
     */
    SolverBase(double rho = 1000.0, double nu = 1.0e-6, double dt = 0.001);

    /// 仮想デストラクタ
    virtual ~SolverBase() = default;

    //==========================================================================
    // 時間発展（純粋仮想関数）
    //==========================================================================

    /**
     * @brief 1タイムステップ進める
     * @param grid 格子
     * @param bc 境界条件
     * @return ソルバー固有の反復回数（-1: 非収束）
     */
    virtual int step(Grid& grid, const BoundaryCondition& bc) = 0;

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
    double getTime() const { return time_; }

    /// 現在のステップ数
    int getStepCount() const { return stepCount_; }

protected:
    double time_ = 0.0;    ///< 現在時刻 [s]
    int stepCount_ = 0;    ///< ステップカウンタ

    //==========================================================================
    // 共通メソッド
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

    /**
     * @brief 中間速度場を計算（予測ステップ共通部分）
     *
     * 圧力項を除いた運動量方程式を解く
     */
    void computeIntermediateVelocity(Grid& grid, const BoundaryCondition& bc);

    /**
     * @brief 圧力勾配による速度補正
     *
     * u^{n+1} = u* - (dt/ρ) * ∇p
     */
    void correctVelocity(Grid& grid);
};

} // namespace fluid
