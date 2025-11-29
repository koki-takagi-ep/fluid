# TVDスキーム（Total Variation Diminishing）

本ドキュメントでは、移流項の高次精度離散化に用いるTVDスキームについて解説する。

## 目次

1. [TVDスキームの原理](#1-tvdスキームの原理)
2. [フラックスリミッター](#2-フラックスリミッター)
3. [実装されているリミッター](#3-実装されているリミッター)
4. [使用方法](#4-使用方法)
5. [ベンチマーク結果](#5-ベンチマーク結果)

---

## 1. TVDスキームの原理

### 1.1 問題の背景

移流方程式の数値解法において：

- **1次風上差分**: 安定だが数値拡散が大きい（$O(\Delta x)$）
- **2次中心差分**: 高精度だが不連続付近で数値振動が発生

TVD（Total Variation Diminishing）スキームは、この問題を解決するために開発された。

### 1.2 Total Variation

離散解 $\{u_i\}$ のTotal Variation（全変動）は次のように定義される：

$$
TV(u) = \sum_i |u_{i+1} - u_i|
$$

物理的に正しい解では、時間発展に伴いTVは増加しない（むしろ減少する）。数値振動が発生すると、TVは人工的に増加する。

### 1.3 TVD条件

スキームがTVDであるとは、次の条件を満たすことを意味する：

$$
TV(u^{n+1}) \leq TV(u^n)
$$

この条件を満たすスキームは、新たな極値を生成せず、数値振動を抑制できる。

### 1.4 MUSCL型再構築

本実装では、MUSCL（Monotone Upstream-centered Schemes for Conservation Laws）型の再構築を用いる：

$$
u_{i+1/2}^L = u_i + \frac{1}{2}\psi(r_i)(u_{i+1} - u_i)
$$

ここで：
- $u_{i+1/2}^L$: セル境界での左側からの外挿値
- $\psi(r)$: リミッター関数
- $r_i = \frac{u_i - u_{i-1}}{u_{i+1} - u_i}$: 勾配比

---

## 2. フラックスリミッター

### 2.1 リミッター関数の役割

リミッター関数 $\psi(r)$ は：
- $r \approx 1$（滑らかな領域）: $\psi \approx 1$（2次精度）
- $r < 0$（極値付近）: $\psi = 0$（1次風上に退化）
- $r \gg 1$ または $r \ll 1$（急勾配）: 適切に制限

### 2.2 TVD領域（Sweby図）

リミッターがTVD条件を満たすには、$(r, \psi)$ 平面で以下の領域内に収まる必要がある：

```
ψ
 ^
2|     ___________
 |    /          |
 |   /           |
1|--+            |
 | /             |
 |/______________|______> r
 0   1           2
```

- 下限: $\psi = 0$ と $\psi = r$
- 上限: $\psi = 2r$ と $\psi = 2$

---

## 3. 実装されているリミッター

### 3.1 Minmod

最も拡散的だが最も安定なリミッター。

$$
\psi(r) = \max(0, \min(1, r))
$$

**特徴**:
- TVD領域の下限に位置
- 数値拡散は大きいが、振動は完全に抑制
- 滑らかな解には不向き

### 3.2 Superbee

最も圧縮的なリミッター。

$$
\psi(r) = \max(0, \min(2r, 1), \min(r, 2))
$$

**特徴**:
- TVD領域の上限に位置
- 数値拡散は最小
- 衝撃波や不連続面の捕捉に優れる
- 滑らかな解では過度に急峻になることがある

### 3.3 Van Leer

滑らかで対称的なリミッター。

$$
\psi(r) = \frac{r + |r|}{1 + |r|}
$$

**特徴**:
- MinmodとSuperbeeの中間的な特性
- 滑らかな遷移を持つ
- 対称性: $\psi(r)/r = \psi(1/r)$
- 汎用的で安定

### 3.4 MC（Monotonized Central）

中心差分に近い高精度リミッター。

$$
\psi(r) = \max\left(0, \min\left(2r, \frac{1+r}{2}, 2\right)\right)
$$

**特徴**:
- 2次精度中心差分の性質を保持
- Superbeeより滑らか
- 滑らかな解に対して高精度

### 3.5 リミッターの比較

| リミッター | 数値拡散 | 振動抑制 | 滑らかさ | 推奨用途 |
|-----------|---------|---------|---------|---------|
| Minmod | 大 | 最強 | 低 | 衝撃波、不連続 |
| Superbee | 最小 | 強 | 低 | 急勾配の保存 |
| Van Leer | 中 | 強 | 高 | 汎用 |
| MC | 中〜小 | 中 | 高 | 滑らかな流れ |

---

## 4. 使用方法

### 4.1 C++コードでの設定

```cpp
#include "Solver.hpp"
#include "FluxLimiter.hpp"

// ソルバーの作成
fluid::Solver solver(rho, nu);

// リミッターの設定
solver.setLimiter(fluid::LimiterType::VanLeer);

// 利用可能なリミッター:
// - fluid::LimiterType::None      (1次風上)
// - fluid::LimiterType::Minmod
// - fluid::LimiterType::Superbee
// - fluid::LimiterType::VanLeer
// - fluid::LimiterType::MC
```

### 4.2 ベンチマークプログラムの実行

```bash
cd build
./tvd_benchmark [nx] [U_lid] [endTime]

# 例: 64x64格子、10mm/s、20秒
./tvd_benchmark 64 0.01 20.0
```

---

## 5. ベンチマーク結果

### 5.1 計算時間比較（64x64格子、Re=100）

| スキーム | 計算時間 [s] | ステップ数 | 最大発散 |
|---------|-------------|-----------|---------|
| 1次風上 | 6.5 | 3277 | 9.3e-7 |
| Minmod | 4.8 | 3277 | 5.4e-7 |
| Superbee | 4.9 | 3277 | 4.2e-7 |
| Van Leer | 4.8 | 3277 | 5.2e-7 |
| MC | 4.8 | 3277 | 5.4e-7 |

TVDスキームは1次風上より高速で収束が良い結果が得られている。

### 5.2 精度に関する注意

TVDスキームは「高次精度」を目指すものだが、リミッターにより局所的に1次精度に落ちる。そのため：
- 滑らかな定常流れでは、1次風上と同程度の精度となることがある
- 非定常流れや急勾配のある流れでは、TVDの利点が顕著になる
- 格子を細かくすると、TVDスキームの精度向上が明確になる

---

## 6. 理論的背景

### 6.1 Godunov定理

**定理**: 線形の保存スキームで2次以上の精度を持つものは、TVDにはなり得ない。

この定理から、TVDスキームは本質的に非線形でなければならない。リミッター関数による勾配の制限が、この非線形性を導入している。

### 6.2 精度次数

TVDスキームの精度は：
- 滑らかな領域: $O(\Delta x^2)$（2次精度）
- 極値付近: $O(\Delta x)$（1次精度に退化）

これは「精度の落ちる」ことを意味するが、振動を抑制するためには必要なトレードオフである。

---

## 参考文献

- Sweby, P. K. (1984). High Resolution Schemes Using Flux Limiters for Hyperbolic Conservation Laws. *SIAM Journal on Numerical Analysis*, 21(5), 995-1011.
- van Leer, B. (1979). Towards the Ultimate Conservative Difference Scheme. V. A Second-Order Sequel to Godunov's Method. *Journal of Computational Physics*, 32(1), 101-136.
- Harten, A. (1983). High Resolution Schemes for Hyperbolic Conservation Laws. *Journal of Computational Physics*, 49(3), 357-393.
- LeVeque, R. J. (2002). *Finite Volume Methods for Hyperbolic Problems*. Cambridge University Press.
