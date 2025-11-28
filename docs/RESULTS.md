# 計算結果

## キャビティ流れ（Lid-driven cavity, Re = 100）

上壁が一定速度で移動する正方形キャビティ内の流れのシミュレーション結果（Projection法）：

![Cavity Flow Result](images/cavity_projection_result.svg)

▶︎ 上から順に：速度場（ベクトル＋カラーマップ）、流線、圧力場、中心線速度分布

▶︎ 主渦が右上に形成され、左下に二次渦が発生している様子が確認できる。

## チャネル流れ（Poiseuille flow, Re = 30）

2枚の平行平板間の圧力駆動流れ。入口で一様流入、出口で自然流出境界条件を設定。

### Projection法による結果

![Channel Flow - Projection](images/channel_projection_result.svg)

▶︎ 上から順に：速度場、流線、圧力場、中心線速度分布

### SIMPLE法による結果

![Channel Flow - SIMPLE](images/channel_simple_result.svg)

▶︎ 上から順に：速度場、流線、圧力場、中心線速度分布

▶︎ 両手法とも、発達した流れでは放物線状の速度分布（Poiseuille流れの理論解）に近づいていることが確認できる。

## ベンチマーク検証

Ghia, Ghia & Shin (1982) のベンチマークデータとの比較により、数値解の妥当性を検証。

### 格子収束性

![Validation Re=100](images/cavity_validation_Re100.svg)

▶︎ 異なる格子解像度（64×64, 128×128）での数値解（ライン）とGhiaの参照データ（ドット）を比較。格子を細かくするにつれて参照解に収束していることが確認できる。

### Projection法 vs SIMPLE法

![Projection vs SIMPLE](images/cavity_validation_comparison.svg)

▶︎ 同一格子（64×64）でのProjection法とSIMPLE法の比較。両手法とも同等の精度でGhiaの参照解に一致している。

## 計算性能

MacBook Air (Apple M2, 24GB RAM) での計算時間：

| ケース | 格子サイズ | 手法 | 総ステップ数 | 計算時間 | 1ステップあたり |
|--------|-----------|------|-------------|---------|----------------|
| Cavity Flow (Re=100) | 64×64 | Projection | 1,639 | 7.6 s | 4.6 ms |
| Channel Flow (Re=30) | 128×32 | Projection | 1,056 | 84 s | 79 ms |
| Channel Flow (Re=30) | 128×32 | SIMPLE | 264 | 113 s | 430 ms |

▶︎ Projection法はSIMPLE法に比べて1ステップあたり約10倍高速（圧力ソルバーの収束特性の違いによる）

## 検証

### キャビティ流れ

Ghia, Ghia & Shin (1982) のベンチマークデータとの比較により検証：
- Re = 100, 400, 1000 で良好な一致を確認
- 中心線速度分布の比較

### チャネル流れ

Poiseuille流れの理論解との比較：

$$
U_{\max} = \frac{3}{2} U_{\text{mean}}
$$

$$
u(y) = U_{\max} \left(1 - \frac{4y^2}{H^2}\right)
$$
