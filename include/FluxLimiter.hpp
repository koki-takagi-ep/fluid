/**
 * @file FluxLimiter.hpp
 * @brief TVD（Total Variation Diminishing）スキーム用フラックスリミッター
 *
 * 高次精度スキームの数値振動を抑制するためのリミッター関数を提供する。
 * 各リミッターは勾配比 r に対して制限係数 φ(r) を返す。
 */

#pragma once

#include <cmath>
#include <algorithm>
#include <string>

namespace fluid {

/**
 * @enum LimiterType
 * @brief 利用可能なフラックスリミッターの種類
 */
enum class LimiterType {
    None,       ///< リミッターなし（1次風上）
    Minmod,     ///< Minmod（最も拡散的、最も安定）
    Superbee,   ///< Superbee（最も圧縮的、急勾配に強い）
    VanLeer,    ///< Van Leer（滑らかで対称的）
    MC          ///< Monotonized Central（中心差分に近い）
};

/**
 * @brief リミッタータイプを文字列に変換
 */
inline std::string limiterTypeToString(LimiterType type) {
    switch (type) {
        case LimiterType::None:     return "None (1st-order upwind)";
        case LimiterType::Minmod:   return "Minmod";
        case LimiterType::Superbee: return "Superbee";
        case LimiterType::VanLeer:  return "Van Leer";
        case LimiterType::MC:       return "MC (Monotonized Central)";
        default:                    return "Unknown";
    }
}

/**
 * @brief 文字列からリミッタータイプに変換
 */
inline LimiterType stringToLimiterType(const std::string& str) {
    if (str == "none" || str == "upwind") return LimiterType::None;
    if (str == "minmod") return LimiterType::Minmod;
    if (str == "superbee") return LimiterType::Superbee;
    if (str == "vanleer" || str == "van_leer") return LimiterType::VanLeer;
    if (str == "mc") return LimiterType::MC;
    return LimiterType::None;
}

/**
 * @namespace limiter
 * @brief フラックスリミッター関数群
 *
 * 各リミッター関数は勾配比 r = (u_i - u_{i-1}) / (u_{i+1} - u_i) に対して
 * 制限係数 φ(r) を返す。TVD条件を満たすため、φ(r) は特定の領域内に収まる必要がある。
 */
namespace limiter {

/**
 * @brief Minmodリミッター
 *
 * 最も拡散的だが最も安定なリミッター。
 * φ(r) = max(0, min(1, r))
 *
 * 特徴:
 * - TVD領域の下限に位置
 * - 数値拡散は大きいが、振動は完全に抑制
 * - 滑らかな解には不向き
 */
inline double minmod(double r) {
    if (r <= 0.0) return 0.0;
    return std::min(1.0, r);
}

/**
 * @brief Superbeeリミッター
 *
 * 最も圧縮的なリミッター。急勾配の保存に優れる。
 * φ(r) = max(0, min(2r, 1), min(r, 2))
 *
 * 特徴:
 * - TVD領域の上限に位置
 * - 数値拡散は最小
 * - 衝撃波や不連続面の捕捉に優れる
 * - 滑らかな解では過度に急峻になることがある
 */
inline double superbee(double r) {
    if (r <= 0.0) return 0.0;
    return std::max(std::min(2.0 * r, 1.0), std::min(r, 2.0));
}

/**
 * @brief Van Leerリミッター
 *
 * 滑らかで対称的なリミッター。
 * φ(r) = (r + |r|) / (1 + |r|)
 *
 * 特徴:
 * - MinmodとSuperbeeの中間的な特性
 * - 滑らかな遷移を持つ
 * - 対称性: φ(r) / r = φ(1/r)
 * - 汎用的で安定
 */
inline double vanLeer(double r) {
    if (r <= 0.0) return 0.0;
    return (r + std::abs(r)) / (1.0 + std::abs(r));
}

/**
 * @brief MC（Monotonized Central）リミッター
 *
 * 中心差分に近い高精度リミッター。
 * φ(r) = max(0, min(2r, (1+r)/2, 2))
 *
 * 特徴:
 * - 2次精度中心差分の性質を保持
 * - Superbeeより滑らか
 * - 滑らかな解に対して高精度
 */
inline double mc(double r) {
    if (r <= 0.0) return 0.0;
    return std::max(0.0, std::min({2.0 * r, 0.5 * (1.0 + r), 2.0}));
}

/**
 * @brief 指定されたタイプのリミッターを適用
 * @param type リミッタータイプ
 * @param r 勾配比
 * @return 制限係数 φ(r)
 */
inline double apply(LimiterType type, double r) {
    switch (type) {
        case LimiterType::Minmod:   return minmod(r);
        case LimiterType::Superbee: return superbee(r);
        case LimiterType::VanLeer:  return vanLeer(r);
        case LimiterType::MC:       return mc(r);
        case LimiterType::None:
        default:                    return 0.0;  // 1次風上に相当
    }
}

/**
 * @brief 符号関数
 */
inline double sign(double x) {
    if (x > 0.0) return 1.0;
    if (x < 0.0) return -1.0;
    return 0.0;
}

/**
 * @brief 2値のminmod関数
 */
inline double minmod2(double a, double b) {
    if (a * b <= 0.0) return 0.0;
    if (std::abs(a) < std::abs(b)) return a;
    return b;
}

} // namespace limiter

} // namespace fluid
