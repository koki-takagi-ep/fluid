# 離散化スキームの詳細

本ドキュメントでは、非圧縮性Navier-Stokes方程式の空間・時間離散化について詳細に解説する。

> **Note**: 各離散化公式の打ち切り誤差（精度次数）は、テイラー展開により導出される。詳細な導出過程は [TAYLOR_EXPANSION.md](TAYLOR_EXPANSION.md) を参照。

## 目次

1. [スタガード格子（MAC格子）](#1-スタガード格子mac格子)
2. [Projection法の離散化](#2-projection法の離散化)
3. [SIMPLE法の離散化](#3-simple法の離散化)
4. [移流項の離散化](#4-移流項の離散化)
5. [拡散項の離散化](#5-拡散項の離散化)
6. [圧力Poisson方程式](#6-圧力poisson方程式)
7. [境界条件の実装](#7-境界条件の実装)
8. [離散化精度のまとめ](#8-離散化精度のまとめ)

---

## 1. スタガード格子（MAC格子）

### 1.1 格子配置

Harlow & Welch (1965) により提案されたMAC法（Marker-and-Cell）では、圧力と速度を異なる位置に配置する**スタガード格子**を使用する。

```
        v[i,j]
          |
    ------+------
    |            |
u[i-1,j]  p[i,j]  u[i,j]
    |            |
    ------+------
        v[i,j-1]
```

- **圧力 $p_{i,j}$**: セル $(i,j)$ の中心
- **x方向速度 $u_{i,j}$**: セル $(i,j)$ と $(i+1,j)$ の境界（右面）
- **y方向速度 $v_{i,j}$**: セル $(i,j)$ と $(i,j+1)$ の境界（上面）

### 1.2 スタガード格子の利点

1. **圧力振動の抑制**: 同一格子点に圧力と速度を配置すると、チェッカーボード状の圧力振動が発生しやすい。スタガード配置ではこれを自然に回避できる。

2. **連続の式の正確な離散化**: セル境界に速度があるため、質量保存（連続の式）を正確に離散化できる。

3. **圧力勾配の自然な評価**: 速度点に隣接する2つの圧力点から直接勾配を計算できる。

### 1.3 配列サイズ

格子サイズを $n_x \times n_y$ とすると：

| 変数 | 配列サイズ | 説明 |
|------|-----------|------|
| $p$ | $(n_x+2) \times (n_y+2)$ | ゴーストセル含む |
| $u$ | $(n_x+1) \times (n_y+2)$ | x方向面、ゴーストセル含む |
| $v$ | $(n_x+2) \times (n_y+1)$ | y方向面、ゴーストセル含む |

---

## 2. Projection法の離散化

### 2.1 概要

Projection法（Chorin, 1968）は、速度場と圧力場を分離して解く分離解法である。

### 2.2 Step 1: 中間速度場の計算

圧力項を除いた運動量方程式を離散化：

$$
\frac{u^*_{i,j} - u^n_{i,j}}{\Delta t} = -[(\mathbf{u} \cdot \nabla) u]^n_{i,j} + \nu [\nabla^2 u]^n_{i,j}
$$

$$
\frac{v^*_{i,j} - v^n_{i,j}}{\Delta t} = -[(\mathbf{u} \cdot \nabla) v]^n_{i,j} + \nu [\nabla^2 v]^n_{i,j}
$$

**時間離散化精度**: $O(\Delta t)$（1次精度、前進Euler法）

### 2.3 Step 2: 圧力Poisson方程式

連続の式 $\nabla \cdot \mathbf{u}^{n+1} = 0$ を満たすための圧力を求める：

$$
\frac{p_{i+1,j} - 2p_{i,j} + p_{i-1,j}}{\Delta x^2} + \frac{p_{i,j+1} - 2p_{i,j} + p_{i,j-1}}{\Delta y^2} = \frac{\rho}{\Delta t}\left(\frac{u^*_{i,j} - u^*_{i-1,j}}{\Delta x} + \frac{v^*_{i,j} - v^*_{i,j-1}}{\Delta y}\right)
$$

**空間離散化精度**:
- 左辺（ラプラシアン）: $O(\Delta x^2, \Delta y^2)$（2次精度中心差分）
- 右辺（発散）: $O(\Delta x, \Delta y)$（1次精度）

> ラプラシアンの2次精度は、前進・後退展開を足し合わせることで1次誤差項が相殺されることから導かれる。[詳細な導出](TAYLOR_EXPANSION.md#3-2階微分の差分公式)

### 2.4 Step 3: 速度の補正

圧力勾配で中間速度を補正：

$$
u^{n+1}_{i,j} = u^*_{i,j} - \frac{\Delta t}{\rho} \frac{p_{i+1,j} - p_{i,j}}{\Delta x} + O(\Delta x)
$$

$$
v^{n+1}_{i,j} = v^*_{i,j} - \frac{\Delta t}{\rho} \frac{p_{i,j+1} - p_{i,j}}{\Delta y} + O(\Delta y)
$$

**空間離散化精度**: $O(\Delta x)$（1次精度前進差分）

---

## 3. SIMPLE法の離散化

### 3.1 概要

SIMPLE法（Semi-Implicit Method for Pressure-Linked Equations, Patankar & Spalding, 1972）は、定常または非定常流れに適用できる圧力-速度連成解法である。

### 3.2 アルゴリズム

1. **仮の速度場 $\mathbf{u}^*$ を計算**（仮の圧力場 $p^*$ を使用）
2. **圧力補正方程式を解く**: $p' = p - p^*$
3. **速度を補正**: $\mathbf{u} = \mathbf{u}^* + \mathbf{u}'$
4. **圧力を更新**: $p = p^* + \alpha_p p'$（緩和係数 $\alpha_p$ を適用）
5. 収束するまで繰り返し

### 3.3 緩和係数

収束を安定化するため、以下の緩和係数を使用：

- **速度緩和係数** $\alpha_u$: 典型値 0.5〜0.8
- **圧力緩和係数** $\alpha_p$: 典型値 0.2〜0.3

推奨: $\alpha_u + \alpha_p \approx 1.0$

### 3.4 Projection法との違い

| 特徴 | Projection法 | SIMPLE法 |
|------|-------------|----------|
| 時間積分 | 陽的 | 反復的 |
| 適用 | 非定常流れ | 定常/非定常流れ |
| 収束性 | 時間ステップで制御 | 緩和係数で制御 |
| 計算コスト | 各ステップ軽い | 反復が必要 |

---

## 4. 移流項の離散化

### 4.1 1次風上差分（Upwind）

移流項 $(\mathbf{u} \cdot \nabla)u$ を風上差分で離散化：

$$
u \frac{\partial u}{\partial x} \approx u_{i,j} \times \begin{cases}
\frac{u_{i,j} - u_{i-1,j}}{\Delta x} + O(\Delta x) & (u_{i,j} \geq 0) \\[8pt]
\frac{u_{i+1,j} - u_{i,j}}{\Delta x} + O(\Delta x) & (u_{i,j} < 0)
\end{cases}
$$

**空間離散化精度**: $O(\Delta x)$（1次精度）

> 風上差分は流れの上流側から情報を取得することで数値安定性を確保する。テイラー展開により、打ち切り誤差が2階微分に比例する**数値拡散**を含むことが示される。[詳細な導出](TAYLOR_EXPANSION.md#4-風上差分の導出)

### 4.2 速度の補間

スタガード格子では、v速度をu点に補間する必要がある（逆も同様）：

$$
v|_{u点} = \frac{1}{4}(v_{i,j-1} + v_{i+1,j-1} + v_{i,j} + v_{i+1,j}) + O(\Delta x^2, \Delta y^2)
$$

**補間精度**: $O(\Delta x^2, \Delta y^2)$（2次精度、4点平均）

### 4.3 数値拡散

1次風上差分のテイラー展開による誤差解析：

$$
\frac{u_i - u_{i-1}}{\Delta x} = \frac{\partial u}{\partial x} - \frac{\Delta x}{2} \frac{\partial^2 u}{\partial x^2} + O(\Delta x^2)
$$

この $\frac{\Delta x}{2} \frac{\partial^2 u}{\partial x^2}$ 項が**数値拡散**を引き起こす。物理的な拡散ではないが、解を人工的に平滑化する効果がある。より高次のスキーム（QUICK、TVDなど）で数値拡散を低減可能。

---

## 5. 拡散項の離散化

### 5.1 2次中心差分

拡散項 $\nu \nabla^2 u$ を中心差分で離散化：

$$
\nu \nabla^2 u \approx \nu \left(\frac{u_{i+1,j} - 2u_{i,j} + u_{i-1,j}}{\Delta x^2} + \frac{u_{i,j+1} - 2u_{i,j} + u_{i,j-1}}{\Delta y^2}\right) + O(\Delta x^2, \Delta y^2)
$$

**空間離散化精度**: $O(\Delta x^2, \Delta y^2)$（2次精度中心差分）

> 2次精度は、$f(x+h)$ と $f(x-h)$ のテイラー展開を足し合わせることで1次の誤差項が相殺され、主誤差項が $\frac{h^2}{12}f''''(x)$ となることから導かれる。[詳細な導出](TAYLOR_EXPANSION.md#31-2次中心差分)

### 5.2 安定性条件

陽的スキームでは、拡散による安定性条件（von Neumann解析より）：

$$
\Delta t \leq \frac{1}{2\nu}\left(\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2}\right)^{-1}
$$

---

## 6. 圧力Poisson方程式

### 6.1 SOR法（Successive Over-Relaxation）

圧力Poisson方程式を反復的に解く：

$$
p^{(k+1)}_{i,j} = (1-\omega)p^{(k)}_{i,j} + \frac{\omega}{a_{i,j}}\left(b_{i,j} - \sum_{m \neq (i,j)} a_m p^{(k)}_m\right)
$$

ここで：
- $\omega$: 緩和係数（$1 < \omega < 2$、典型値: 1.8）
- $a_{i,j} = -2/\Delta x^2 - 2/\Delta y^2$: 対角成分

### 6.2 Red-Black SOR

並列化のため、Red-Black順序で更新：

1. **Red点**（$i+j$ が偶数）を更新
2. **Black点**（$i+j$ が奇数）を更新

### 6.3 収束判定

残差ベースの収束判定を使用：

$$
\max_{i,j} |(\nabla^2 p)_{i,j} - \text{RHS}_{i,j}| < \epsilon
$$

ここで $\epsilon$ は発散の閾値（典型値: $10^{-6}$）。

---

## 7. 境界条件の実装

### 7.1 滑りなし壁（No-slip）

壁面で速度がゼロ：

- **接線速度**: ゴーストセルを使用して $u_{ghost} = -u_{interior}$
- **法線速度**: 壁面で直接 $u_n = 0$

### 7.2 流入境界（Inflow）

指定した速度分布を設定：

$$
u|_{inlet} = U_{in}, \quad v|_{inlet} = 0
$$

### 7.3 流出境界（Outflow）

ゼロ勾配条件（自然境界条件）：

$$
\frac{\partial u}{\partial n} = 0, \quad \frac{\partial v}{\partial n} = 0
$$

### 7.4 圧力境界条件

- **壁面**: Neumann条件 $\frac{\partial p}{\partial n} = 0$
- **流出**: 参照圧力を設定（例: $p = p_{atm}$）

---

## 8. 離散化精度のまとめ

本ソルバーで使用する離散化スキームの精度一覧：

| 項目 | 離散化手法 | 空間精度 | 時間精度 |
|------|-----------|---------|---------|
| 移流項 | 1次風上差分 | $O(\Delta x)$ | - |
| 拡散項 | 2次中心差分 | $O(\Delta x^2)$ | - |
| 圧力ラプラシアン | 2次中心差分 | $O(\Delta x^2)$ | - |
| 圧力勾配 | 1次前進差分 | $O(\Delta x)$ | - |
| 速度発散 | 1次差分 | $O(\Delta x)$ | - |
| 時間積分 | 前進Euler | - | $O(\Delta t)$ |

> 各離散化公式の精度は、テイラー展開により厳密に導出される。詳細は [差分公式のテイラー展開による導出](TAYLOR_EXPANSION.md) を参照。

---

## 参考文献

- Harlow, F. H., & Welch, J. E. (1965). Numerical calculation of time-dependent viscous incompressible flow of fluid with free surface. *Physics of Fluids*, 8(12), 2182-2189.
- Chorin, A. J. (1968). Numerical solution of the Navier-Stokes equations. *Mathematics of Computation*, 22(104), 745-762.
- Patankar, S. V., & Spalding, D. B. (1972). A calculation procedure for heat, mass and momentum transfer in three-dimensional parabolic flows. *International Journal of Heat and Mass Transfer*, 15(10), 1787-1806.
- Ferziger, J. H., & Peric, M. (2002). *Computational Methods for Fluid Dynamics* (3rd ed.). Springer.
- LeVeque, R. J. (2007). *Finite Difference Methods for Ordinary and Partial Differential Equations*. SIAM.

詳細な参考文献リストは [REFERENCES.md](REFERENCES.md) を参照。
