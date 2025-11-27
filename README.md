# Incompressible Navier-Stokes Solver

2次元非圧縮性Navier-Stokes方程式のための数値解法ライブラリ。
MAC法（Marker-and-Cell）とProjection法に基づく実装。

## 支配方程式

### Navier-Stokes方程式

$$
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu \nabla^2 \mathbf{u}
$$

$$
\nabla \cdot \mathbf{u} = 0
$$

ここで：
- $\mathbf{u} = (u, v)$: 速度ベクトル [m/s]
- $p$: 圧力 [Pa]
- $\rho$: 密度 [kg/m³]
- $\nu$: 動粘性係数 [m²/s]

## 数値解法

### 時間積分（Projection法）

1. **予測ステップ**: 圧力項を除いた運動量方程式を解く
   $$
   \frac{\mathbf{u}^* - \mathbf{u}^n}{\Delta t} = -(\mathbf{u}^n \cdot \nabla)\mathbf{u}^n + \nu \nabla^2 \mathbf{u}^n
   $$

2. **圧力Poisson方程式**: 連続の式を満たす圧力場を求める
   $$
   \nabla^2 p^{n+1} = \frac{\rho}{\Delta t} \nabla \cdot \mathbf{u}^*
   $$

3. **修正ステップ**: 速度場を更新
   $$
   \mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho} \nabla p^{n+1}
   $$

### 空間離散化

- **格子**: スタガード格子（MAC格子）
  - 圧力 $p$: セル中心
  - 速度 $u$: セル面（x方向）
  - 速度 $v$: セル面（y方向）

- **移流項**: 1次風上差分
- **拡散項**: 2次中心差分
- **圧力Poisson方程式**: SOR法（Successive Over-Relaxation）

### 安定性条件

CFLおよび粘性条件による適応的時間刻み：

$$
\Delta t \leq \min\left( \frac{\Delta x}{|u|_{\max}}, \frac{\Delta y}{|v|_{\max}}, \frac{\Delta x^2}{4\nu}, \frac{\Delta y^2}{4\nu} \right)
$$

## プロジェクト構成

```
fluid/
├── include/              # ヘッダファイル
│   ├── Grid.hpp          # 格子・速度場・圧力場
│   ├── Solver.hpp        # 時間発展ソルバー
│   ├── PressureSolver.hpp # 圧力Poisson方程式ソルバー
│   ├── BoundaryCondition.hpp # 境界条件
│   └── CSVWriter.hpp     # データ出力
├── src/                  # ソースファイル
├── examples/             # 計算例
│   ├── cavity_flow.cpp   # キャビティ流れ（Lid-driven cavity）
│   └── channel_flow.cpp  # チャネル流れ（Poiseuille flow）
├── scripts/              # 可視化・解析スクリプト
│   ├── visualize.py      # 結果の可視化
│   ├── validation.py     # ベンチマーク検証
│   └── convergence.py    # 収束解析
└── docs/                 # ドキュメント
    └── REFERENCES.md     # 参考文献
```

## ビルド

```bash
mkdir build && cd build
cmake ..
make -j4
```

## 使用例

### キャビティ流れ（Lid-driven cavity）

```bash
./cavity_flow [nx] [Re] [end_time]
# 例: ./cavity_flow 128 100 20.0
```

上壁が一定速度で移動する正方形キャビティ内の流れ。
Ghia et al. (1982) のベンチマークデータとの比較検証が可能。

### チャネル流れ（Poiseuille flow）

```bash
./channel_flow [nx] [ny] [U_in] [end_time]
# 例: ./channel_flow 256 32 0.01 5.0
```

2枚の平行平板間の圧力駆動流れ。
理論解（放物線速度分布）との比較が可能。

## 可視化

```bash
cd build

# 最終状態の可視化
python ../scripts/visualize.py output_cavity --plot-final

# アニメーション生成
python ../scripts/visualize.py output_channel --animation --save animation.gif

# ベンチマーク検証（キャビティ流れ）
python ../scripts/validation.py output_cavity

# 収束解析
python ../scripts/convergence.py output_channel
```

## 検証

### キャビティ流れ

Ghia, Ghia & Shin (1982) のベンチマークデータとの比較により検証：
- Re = 100, 400, 1000 で良好な一致を確認
- 中心線速度分布の比較

### チャネル流れ

Poiseuille流れの理論解との比較：
- 最大速度: $U_{\max} = \frac{3}{2} U_{\text{mean}}$
- 速度分布: $u(y) = U_{\max} \left(1 - \frac{4y^2}{H^2}\right)$

## 参考文献

詳細は [docs/REFERENCES.md](docs/REFERENCES.md) を参照。

## ライセンス

MIT License
