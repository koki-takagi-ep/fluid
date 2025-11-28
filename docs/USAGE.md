# 使用方法

## ビルド

```bash
mkdir build && cd build
cmake ..
make -j4
```

## シミュレーション実行

### キャビティ流れ（Lid-driven cavity）

```bash
./cavity_flow [nx] [U_lid] [end_time]
# 例: ./cavity_flow 64 0.01 10.0
```

上壁が一定速度で移動する正方形キャビティ内の流れ。
Ghia et al. (1982) のベンチマークデータとの比較検証が可能。

### チャネル流れ - Projection法

```bash
./channel_flow [nx] [ny] [U_in] [end_time]
# 例: ./channel_flow 128 32 0.01 2.0
```

### チャネル流れ - SIMPLE法

```bash
./channel_flow_simple [nx] [ny] [U_in] [end_time]
# 例: ./channel_flow_simple 128 32 0.01 1.0
```

2枚の平行平板間の圧力駆動流れ。
理論解（放物線速度分布）との比較が可能。

## 可視化

```bash
cd build

# 最終状態の可視化（Projection法）
python ../scripts/visualize.py output/cavity_projection --plot-final
python ../scripts/visualize.py output/channel_projection --plot-final

# SIMPLE法の結果
python ../scripts/visualize.py output/channel_simple --plot-final

# アニメーション生成
python ../scripts/visualize.py output/channel_projection --animation

# ベンチマーク検証（キャビティ流れ）
python ../scripts/validation.py output/cavity_projection

# 複数ケースの比較検証
python ../scripts/validation.py output/cavity_64x64 output/cavity_128x128 \
  --labels "64×64" "128×128" --Re 100

# 収束解析
python ../scripts/convergence.py output/channel_projection
```

## 出力ファイル

シミュレーション結果は以下の構造で出力されます：

```
output/
├── cavity_projection/      # Cavity flow (Projection method)
├── channel_projection/     # Channel flow (Projection method)
└── channel_simple/         # Channel flow (SIMPLE method)
    ├── data/
    │   ├── field_000000.csv    # x, y, u, v, p, magnitude
    │   ├── field_000001.csv
    │   └── metadata.csv        # Grid parameters
    ├── figures/
    │   ├── result.svg          # SVG format
    │   ├── result.pdf          # PDF format
    │   └── result.png          # PNG format (600 dpi)
    └── simulation.log          # Computation time and settings
```

## Python依存関係

```bash
pip install numpy matplotlib pandas scipy
```
