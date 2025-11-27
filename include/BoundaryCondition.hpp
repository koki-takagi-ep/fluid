/**
 * @file BoundaryCondition.hpp
 * @brief 境界条件の定義と適用
 *
 * 速度場と圧力場に対する各種境界条件を実装。
 */

#pragma once

#include "Grid.hpp"

namespace fluid {

/**
 * @enum BCType
 * @brief 境界条件の種類
 */
enum class BCType {
    NoSlip,      ///< 滑りなし壁（Dirichlet条件）: u = u_wall
    Slip,        ///< 滑り壁（Neumann条件）: ∂u_t/∂n = 0, u_n = 0
    Inflow,      ///< 流入境界: u = u_in（指定値）
    Outflow,     ///< 流出境界: ∂u/∂n = 0（Neumann条件）
    Periodic     ///< 周期境界
};

/**
 * @class BoundaryCondition
 * @brief 境界条件クラス
 *
 * 4つの境界（左、右、下、上）に対する条件を管理。
 *
 * ## 実装されている境界条件
 *
 * ### NoSlip（滑りなし壁）
 * 壁面で速度がゼロ（または壁面速度）になる条件。
 * - ゴーストセルを使用した反射条件で実装
 *
 * ### Inflow（流入）
 * 流入速度を指定するDirichlet条件。
 * - u = u_in, v = 0
 *
 * ### Outflow（流出）
 * 流出境界でのNeumann条件。
 * - ∂u/∂n = 0（速度勾配ゼロ）
 *
 * ## 使用例
 * @code
 * // キャビティ流れ（上壁が移動）
 * auto bc = BoundaryCondition::cavityFlow(1.0);  // U_lid = 1.0 m/s
 *
 * // チャネル流れ（左から流入）
 * auto bc = BoundaryCondition::channelFlow(0.01);  // U_in = 0.01 m/s
 * @endcode
 */
class BoundaryCondition {
public:
    //==========================================================================
    // 境界条件タイプ
    //==========================================================================

    BCType left   = BCType::NoSlip;   ///< 左境界 (x = 0)
    BCType right  = BCType::NoSlip;   ///< 右境界 (x = Lx)
    BCType bottom = BCType::NoSlip;   ///< 下境界 (y = 0)
    BCType top    = BCType::NoSlip;   ///< 上境界 (y = Ly)

    //==========================================================================
    // 境界での速度値（Dirichlet条件用）
    //==========================================================================

    double u_left = 0.0;    ///< 左境界でのu速度 [m/s]
    double u_right = 0.0;   ///< 右境界でのu速度 [m/s]
    double u_bottom = 0.0;  ///< 下境界でのu速度 [m/s]
    double u_top = 0.0;     ///< 上境界でのu速度 [m/s]

    double v_left = 0.0;    ///< 左境界でのv速度 [m/s]
    double v_right = 0.0;   ///< 右境界でのv速度 [m/s]
    double v_bottom = 0.0;  ///< 下境界でのv速度 [m/s]
    double v_top = 0.0;     ///< 上境界でのv速度 [m/s]

    //==========================================================================
    // コンストラクタ
    //==========================================================================

    /// デフォルトコンストラクタ（全境界NoSlip）
    BoundaryCondition() = default;

    //==========================================================================
    // 境界条件の適用
    //==========================================================================

    /**
     * @brief 速度場に境界条件を適用
     * @param grid 格子（速度場が更新される）
     */
    void apply(Grid& grid) const;

    /**
     * @brief 圧力場に境界条件を適用
     *
     * 全境界でNeumann条件（∂p/∂n = 0）を適用。
     *
     * @param grid 格子（圧力場が更新される）
     */
    void applyPressureBC(Grid& grid) const;

    //==========================================================================
    // ファクトリメソッド（典型的な境界条件設定）
    //==========================================================================

    /**
     * @brief Lid-driven cavity flow の境界条件
     *
     * - 上壁: NoSlip with u = lidVelocity
     * - 他の壁: NoSlip (u = v = 0)
     *
     * @param lidVelocity 上壁の速度 [m/s]
     * @return 境界条件オブジェクト
     */
    static BoundaryCondition cavityFlow(double lidVelocity = 1.0);

    /**
     * @brief Channel flow (Poiseuille flow) の境界条件
     *
     * - 左: Inflow (u = inflowVelocity)
     * - 右: Outflow (∂u/∂x = 0)
     * - 上下: NoSlip (u = v = 0)
     *
     * @param inflowVelocity 流入速度 [m/s]
     * @return 境界条件オブジェクト
     */
    static BoundaryCondition channelFlow(double inflowVelocity = 1.0);

private:
    /// 左境界条件を適用
    void applyLeftBC(Grid& grid) const;

    /// 右境界条件を適用
    void applyRightBC(Grid& grid) const;

    /// 下境界条件を適用
    void applyBottomBC(Grid& grid) const;

    /// 上境界条件を適用
    void applyTopBC(Grid& grid) const;
};

} // namespace fluid
