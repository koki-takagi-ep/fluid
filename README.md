# Incompressible Navier-Stokes Solver

2次元非圧縮性Navier-Stokes方程式のための数値解法ライブラリ。
MAC法（Marker-and-Cell）による空間離散化と、**Projection法**および**SIMPLE法**による時間積分を実装。

## プロジェクト構成

```
fluid/
├── include/                  # ヘッダファイル
│   ├── Grid.hpp              # 格子・速度場・圧力場
│   ├── SolverBase.hpp        # ソルバー基底クラス
│   ├── Solver.hpp            # Projection法ソルバー
│   ├── SimpleSolver.hpp      # SIMPLE法ソルバー
│   ├── PressureSolver.hpp    # 圧力Poisson方程式ソルバー
│   ├── BoundaryCondition.hpp # 境界条件
│   └── CSVWriter.hpp         # データ出力
├── src/                      # ソースファイル
├── examples/                 # 計算例
│   ├── cavity_flow.cpp       # キャビティ流れ（Projection法）
│   ├── cavity_flow_simple.cpp # キャビティ流れ（SIMPLE法）
│   ├── channel_flow.cpp      # チャネル流れ（Projection法）
│   └── channel_flow_simple.cpp # チャネル流れ（SIMPLE法）
├── scripts/                  # 可視化・解析スクリプト
│   ├── visualize.py          # 結果の可視化
│   ├── validation.py         # ベンチマーク検証
│   └── convergence.py        # 収束解析
└── docs/                     # ドキュメント
```

## ドキュメント

| ドキュメント | 内容 |
|-------------|------|
| [使用方法](docs/USAGE.md) | ビルド、シミュレーション実行、可視化 |
| [計算結果](docs/RESULTS.md) | シミュレーション結果と検証 |
| [離散化スキーム](docs/DISCRETIZATION.md) | 空間・時間離散化の詳細 |
| [参考文献](docs/REFERENCES.md) | 参考文献一覧 |

## クイックスタート

```bash
# ビルド
mkdir build && cd build
cmake .. && make -j4

# キャビティ流れを実行
./cavity_flow 64 0.01 10.0

# 結果を可視化
python ../scripts/visualize.py output/cavity_projection --plot-final
```

## 支配方程式

### Navier-Stokes方程式

$$
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2 \mathbf{u}
$$

$$
\nabla \cdot \mathbf{u} = 0
$$

ここで：
- $\mathbf{u} = (u, v)$: 速度ベクトル (m/s)
- $p$: 圧力 (Pa)
- $\rho$: 密度 (kg/m³)
- $\nu$: 動粘性係数 (m²/s)

## 数値解法

### Projection法（時間積分）

Projection法（Chorin, 1968）は、非圧縮性Navier-Stokes方程式の時間積分において最も広く使われる手法の一つである。この手法の核心は、**速度場と圧力場を分離して解く**ことで、計算効率を大幅に向上させる点にある。

#### なぜProjection法を使うのか？

非圧縮性流れでは、速度と圧力が連続の式 $\nabla \cdot \mathbf{u} = 0$ を通じて強く結合している。これを直接解くことは計算コストが高い。Projection法は**Helmholtz-Hodge分解**の原理に基づき、任意のベクトル場を発散なし成分と勾配成分に分解できることを利用する：

$$
\mathbf{w} = \mathbf{u} + \nabla \phi
$$

ここで $\nabla \cdot \mathbf{u} = 0$ である。この性質を使い、まず圧力を無視した「仮の速度場」を計算し、その後圧力によって連続の式を満たすように補正する。

#### アルゴリズム

**Step 1. 予測ステップ（Predictor）**

圧力項を除いた運動量方程式を解き、中間速度場 $\mathbf{u}^*$ を求める：

$$
\frac{\mathbf{u}^* - \mathbf{u}^n}{\Delta t} = -(\mathbf{u}^n \cdot \nabla)\mathbf{u}^n + \nu \nabla^2 \mathbf{u}^n
$$

この段階では連続の式 $\nabla \cdot \mathbf{u}^* = 0$ は一般に満たされない。

**Step 2. 圧力Poisson方程式**

次の時刻で連続の式を満たす圧力場を求める。修正ステップの式の発散をとり $\nabla \cdot \mathbf{u}^{n+1} = 0$ を課すと：

$$
\nabla^2 p^{n+1} = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*
$$

この楕円型偏微分方程式をSOR法（逐次過緩和法）で反復的に解く。

**Step 3. 修正ステップ（Corrector）**

求めた圧力勾配で中間速度場を補正し、発散なしの速度場を得る：

$$
\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla p^{n+1}
$$

### SIMPLE法（時間積分）

SIMPLE法（Semi-Implicit Method for Pressure-Linked Equations, Patankar & Spalding, 1972）は、定常・非定常流れの両方に適用できる圧力-速度連成解法である。

#### Projection法との比較

| 特徴 | Projection法 | SIMPLE法 |
|------|-------------|----------|
| 時間積分 | 陽的 | 反復的 |
| 適用 | 非定常流れ | 定常/非定常流れ |
| 収束制御 | 時間ステップ | 緩和係数 |

### 空間離散化

#### スタガード格子（MAC格子）

Harlow & Welch (1965) により提案されたMAC法では、圧力と速度を異なる位置に配置する：

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

スタガード格子の利点：
1. チェッカーボード状の圧力振動を自然に抑制
2. 質量保存（連続の式）を正確に離散化
3. 圧力勾配を直接評価可能

#### 移流項（1次風上差分）

$$
u \frac{\partial u}{\partial x} \approx u_{i,j} \times \begin{cases}
\frac{u_{i,j} - u_{i-1,j}}{\Delta x} & (u_{i,j} \geq 0) \\
\frac{u_{i+1,j} - u_{i,j}}{\Delta x} & (u_{i,j} < 0)
\end{cases}
$$

#### 拡散項（2次中心差分）

$$
\nu \nabla^2 u \approx \nu \left(\frac{u_{i+1,j} - 2u_{i,j} + u_{i-1,j}}{\Delta x^2} + \frac{u_{i,j+1} - 2u_{i,j} + u_{i,j-1}}{\Delta y^2}\right)
$$

#### 圧力Poisson方程式（SOR法）

$$
p^{(k+1)}_{i,j} = (1-\omega)p^{(k)}_{i,j} + \frac{\omega}{a_{i,j}}\left(b_{i,j} - \sum_{m \neq (i,j)} a_m p^{(k)}_m\right)
$$

- $\omega$: 緩和係数（$1 < \omega < 2$、典型値: 1.8）

### 安定性条件

CFLおよび粘性条件による適応的時間刻み：

$$
\Delta t \leq \min\left( \frac{\Delta x}{|u|_{\max}}, \frac{\Delta y}{|v|_{\max}}, \frac{\Delta x^2}{4\nu}, \frac{\Delta y^2}{4\nu} \right)
$$

詳細は [離散化スキーム](docs/DISCRETIZATION.md) を参照。

## 計算結果

### キャビティ流れ（Re = 100）

![Cavity Flow Result](docs/images/cavity_projection_result.svg)

### チャネル流れ（Re = 30）

**Projection法**

![Channel Flow - Projection](docs/images/channel_projection_result.svg)

**SIMPLE法**

![Channel Flow - SIMPLE](docs/images/channel_simple_result.svg)

### ベンチマーク検証

Ghia, Ghia & Shin (1982) のベンチマークデータとの比較：

![Validation](docs/images/cavity_validation_Re100.svg)

詳細は [計算結果](docs/RESULTS.md) を参照。

## ライセンス

MIT License
