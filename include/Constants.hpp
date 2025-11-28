/**
 * @file Constants.hpp
 * @brief 数値計算用の定数定義
 *
 * マジックナンバーを避けるための名前付き定数を定義。
 */

#pragma once

namespace fluid {
namespace constants {

//==============================================================================
// 物理定数
//==============================================================================

/// 大気圧 [Pa]
constexpr double ATMOSPHERIC_PRESSURE = 101325.0;

/// 水の密度 [kg/m³]
constexpr double WATER_DENSITY = 1000.0;

/// 水の動粘性係数 [m²/s]
constexpr double WATER_KINEMATIC_VISCOSITY = 1.0e-6;

//==============================================================================
// 数値計算パラメータ
//==============================================================================

/// 速度がゼロとみなされる閾値
constexpr double VELOCITY_ZERO_THRESHOLD = 1e-10;

/// 時間刻みの初期最大値
constexpr double TIMESTEP_MAX_INITIAL = 1e10;

/// 粘性CFL条件の係数
constexpr double VISCOUS_CFL_FACTOR = 0.5;

/// デフォルトのCFL数
constexpr double DEFAULT_CFL = 0.5;

/// デフォルトの時間刻み [s]
constexpr double DEFAULT_TIMESTEP = 0.001;

//==============================================================================
// 圧力ソルバーパラメータ
//==============================================================================

/// SORの緩和係数（デフォルト）
constexpr double DEFAULT_SOR_OMEGA = 1.8;

/// 圧力ソルバーの最大反復回数（デフォルト）
constexpr int DEFAULT_PRESSURE_MAX_ITERATIONS = 10000;

/// 圧力ソルバーの収束判定閾値（デフォルト）
constexpr double DEFAULT_PRESSURE_TOLERANCE = 1e-5;

//==============================================================================
// 差分スキームの係数
//==============================================================================

/// ラプラシアン中心係数
constexpr double LAPLACIAN_CENTER_COEFF = 2.0;

/// 4点平均の係数
constexpr double FOUR_POINT_AVERAGE_COEFF = 0.25;

/// 2点補間の係数
constexpr double INTERPOLATION_COEFF = 0.5;

} // namespace constants
} // namespace fluid
