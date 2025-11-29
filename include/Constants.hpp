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

//==============================================================================
// SIMPLE法パラメータ
//==============================================================================

/// 速度の緩和係数（デフォルト）
constexpr double DEFAULT_VELOCITY_RELAXATION = 0.7;

/// 圧力の緩和係数（デフォルト）
constexpr double DEFAULT_PRESSURE_RELAXATION = 0.3;

/// SIMPLE法の最大外部反復回数
constexpr int DEFAULT_SIMPLE_MAX_ITERATIONS = 100;

/// SIMPLE法の収束判定閾値
constexpr double DEFAULT_SIMPLE_CONVERGENCE_TOL = 1e-4;

//==============================================================================
// 物理パラメータ（流体力学）
//==============================================================================

/// 水の動粘度 [Pa·s]
constexpr double WATER_DYNAMIC_VISCOSITY = 1.0e-3;

/// Hagen-Poiseuille流れの平均速度/最大速度比（2/3）
constexpr double POISEUILLE_MEAN_MAX_RATIO = 2.0 / 3.0;

//==============================================================================
// デフォルトシミュレーションパラメータ
//==============================================================================

/// デフォルトのキャビティ流れ終了時間 [s]
constexpr double DEFAULT_CAVITY_END_TIME = 10.0;

/// デフォルトのチャネル流れ終了時間 [s]
constexpr double DEFAULT_CHANNEL_END_TIME = 5.0;

} // namespace constants
} // namespace fluid
