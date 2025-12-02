# Literature Survey: 電気浸透流とイオン輸送のベンチマーク問題

本ドキュメントでは、電解液・イオン液体における電気浸透流（Electroosmotic Flow, EOF）およびイオン輸送シミュレーションで使用されるベンチマーク問題について文献調査を行う。

## 目次

1. [背景と動機](#1-背景と動機)
2. [支配方程式](#2-支配方程式)
3. [古典的ベンチマーク問題](#3-古典的ベンチマーク問題)
4. [数値検証に用いられる問題](#4-数値検証に用いられる問題)
5. [オープンソースソルバー](#5-オープンソースソルバー)
6. [ライセンスに関する考察](#6-ライセンスに関する考察)
7. [参考文献](#7-参考文献)

---

## 1. 背景と動機

### 1.1 なぜ電気浸透流が重要か

電気浸透流は、マイクロ・ナノ流体デバイス、電気化学システム、電池・キャパシタ、生体医療デバイスなど幅広い分野で重要な役割を果たす。特に：

- **Lab-on-a-chip**: DNA分離、タンパク質分析
- **電池・スーパーキャパシタ**: 電解液中のイオン輸送
- **脱塩・水処理**: 電気透析、容量性脱塩
- **燃料電池**: プロトン交換膜内のイオン輸送

### 1.2 Hagen-Poiseuille流れとの類似性

本ソルバーで検証に用いているHagen-Poiseuille流れ（チャネル流れ）は、電気浸透流との類似性がある：

| 特徴 | Hagen-Poiseuille | 電気浸透流 |
|------|-----------------|-----------|
| 駆動力 | 圧力勾配 $\nabla p$ | 電場 $\mathbf{E}$ |
| 速度分布（完全発達） | 放物線 | プラグ流（EDL外）|
| 支配方程式 | Navier-Stokes | Navier-Stokes + Poisson-Nernst-Planck |

電気浸透流では、電気二重層（EDL: Electric Double Layer）の存在が本質的であり、EDLの厚さ（Debye長）と流路スケールの比が重要なパラメータとなる。

---

## 2. 支配方程式

### 2.1 Poisson-Nernst-Planck-Navier-Stokes (PNPNS) 方程式

電気浸透流の完全な記述には、以下の連成方程式系が必要：

**Poisson方程式**（電位分布）:
$$
\nabla^2 \phi = -\frac{\rho_e}{\varepsilon}
$$

**Nernst-Planck方程式**（イオン濃度）:
$$
\frac{\partial c_i}{\partial t} + \nabla \cdot (c_i \mathbf{u}) = D_i \nabla^2 c_i + \frac{z_i F D_i}{RT} \nabla \cdot (c_i \nabla \phi)
$$

**Navier-Stokes方程式**（流体流れ）:
$$
\rho \left( \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u} \right) = -\nabla p + \mu \nabla^2 \mathbf{u} + \rho_e \mathbf{E}
$$

ここで：
- $\phi$: 電位
- $c_i$: イオン種 $i$ の濃度
- $\rho_e = F \sum_i z_i c_i$: 電荷密度
- $z_i$: イオン価数
- $D_i$: 拡散係数
- $F$: Faraday定数
- $\mathbf{E} = -\nabla \phi$: 電場

### 2.2 簡略化モデル

多くのベンチマーク問題では、以下の簡略化が用いられる：

1. **Poisson-Boltzmann近似**: 平衡状態を仮定し、Nernst-Planck方程式をBoltzmann分布で置換
2. **Debye-Hückel近似**: 低電位（$ze\phi \ll k_B T$）での線形化
3. **Helmholtz-Smoluchowski近似**: 薄いEDL（$\kappa a \gg 1$）での速度スリップ境界条件

---

## 3. 古典的ベンチマーク問題

### 3.1 Helmholtz-Smoluchowski速度

**最も基本的なベンチマーク**

薄いEDL極限での電気浸透流速度（理論解）：

$$
u_{HS} = -\frac{\varepsilon \zeta E}{\mu}
$$

ここで $\zeta$ はゼータ電位、$E$ は印加電場。

**検証方法**:
- 直線チャネル内のEOF速度を計算し、Helmholtz-Smoluchowski速度と比較
- EDL厚さを十分小さくした極限で一致を確認

**参考文献**:
- Helmholtz (1879), Smoluchowski (1903)
- Park & Lee (2008) "Extension of the Helmholtz-Smoluchowski velocity to hydrophobic microchannels"

### 3.2 有限EDL厚さでの解析解

**2次元直線チャネル**

Debye-Hückel近似下での解析解（Rice & Whitehead, 1965）：

$$
u(y) = u_{HS} \left( 1 - \frac{\cosh(\kappa y)}{\cosh(\kappa H/2)} \right)
$$

ここで $\kappa^{-1}$ はDebye長、$H$ はチャネル幅。

**検証方法**:
- EDL厚さ/チャネル幅の比を変えて速度分布を比較
- 薄いEDL極限でプラグ流、厚いEDL極限で放物線分布に近づくことを確認

### 3.3 混合EOF/圧力駆動流

電気浸透流と圧力駆動流の重ね合わせ（Santiago, 2001）：

$$
u(y) = u_{HS} \left( 1 - \frac{\cosh(\kappa y)}{\cosh(\kappa H/2)} \right) + \frac{H^2}{8\mu} \frac{dp}{dx} \left( 1 - \frac{4y^2}{H^2} \right)
$$

**検証方法**:
- 圧力勾配と電場を独立に変化させ、理論解との一致を確認

---

## 4. 数値検証に用いられる問題

### 4.1 1次元イオン輸送（平板間）

**問題設定**:
- 2枚の平行平板間のイオン濃度・電位分布
- 定常状態での解析解が存在

**支配方程式**: 1D Poisson-Nernst-Planck

**検証項目**:
- 濃度分布の解析解との比較
- 電位分布の解析解との比較
- 電流-電圧特性

**参考文献**: Bazant et al. (2004) "Diffuse-charge dynamics in electrochemical systems"

### 4.2 導電性円柱周りの誘起電荷電気浸透（ICEO）

**問題設定** (Squires & Bazant, 2004):
- 一様電場中に置かれた導電性円柱
- 円柱表面での誘起電荷による渦流れ

**特徴**:
- 非線形現象（誘起電荷は印加電場の2乗に比例）
- 定常・非定常の両方で検証可能
- 解析的な漸近解が存在

**検証項目**:
- 円柱周りの流れ構造（4つの渦）
- 流れの強度の電場依存性

### 4.3 イオン選択性膜での電気対流

**問題設定** (Rubinstein & Zaltzman, 2000):
- イオン選択性膜表面での電気対流不安定性
- 限界電流以上での渦構造形成

**特徴**:
- 電気水力学的不安定性
- 乱流的な流れパターン
- 実験との比較が可能

**参考文献**:
- Pham et al. (2012) "Direct numerical simulation of electroconvective instability"
- Druzgalski et al. (2013) "Direct numerical simulation of electroconvective instability and hydrodynamic chaos"

### 4.4 ナノポア内のイオン輸送

**問題設定**:
- 単一ナノポア（円筒形または円錐形）を通るイオン輸送
- 電流整流効果（円錐形ポア）

**検証項目**:
- 電流-電圧特性
- イオン濃度分布
- 整流比（円錐形ポア）

**参考文献**:
- Liu et al. (2007) "Electrokinetic effects in nanochannels"
- Vlassiouk & Siwy (2007) "Nanofluidic diode"

### 4.5 イオン液体の電気動力学

**問題設定** (Bazant et al., 2011):
- イオン液体中の電気二重層（crowding効果）
- 修正Poisson-Nernst-Planck方程式

**特徴**:
- 古典的なDebye-Hückel理論が適用不可
- イオンサイズ効果、短距離相関が重要
- 実験データとの比較が可能

**参考文献**:
- Bazant et al. (2011) "Double layer in ionic liquids: Overscreening versus crowding"
- Fedorov & Kornyshev (2014) "Ionic liquids at electrified interfaces"

---

## 5. オープンソースソルバー

### 5.1 OpenFOAM ベース

**rheoTool** (Pimenta & Alves, 2017):
- 電気動力学的流れに対応したOpenFOAMソルバー
- Poisson-Nernst-Planck + Navier-Stokes
- GPL v3ライセンス

**pnpFoam** (Schaurer et al., 2023):
- ポアスケールでのStokes-Poisson-Nernst-Planck
- OpenFOAMベース
- オープンソース

### 5.2 その他のツール

**COMSOL Multiphysics**:
- 商用ソフトウェア
- 広範なベンチマーク検証済み
- 多くの論文で参照される「準標準」

**FEniCS**:
- 有限要素法ライブラリ
- PNPNS方程式の実装例あり
- LGPL/GPLライセンス

---

## 6. ライセンスに関する考察

### 6.1 学術研究でのライセンス選択

**目的別の推奨ライセンス**:

| 目的 | 推奨ライセンス | 理由 |
|------|--------------|------|
| 最大限の普及 | MIT, Apache 2.0 | 制限が少なく、商用利用も可能 |
| 派生物のオープン性維持 | GPL v3 | コピーレフト条項により派生物もオープンに |
| 特許保護 | Apache 2.0 | 明示的な特許条項あり |
| 学術引用の確保 | MIT + 引用要請 | ライセンス + READMEでの引用要請 |

### 6.2 CC BY-NC-NDとソフトウェアライセンス

CC BY-NC-ND 4.0（学術文献でよく使用）は**ソフトウェアには推奨されない**：

1. **非互換性**: Open Source Initiative (OSI) の定義を満たさない
2. **改変禁止**: ソフトウェア開発に本質的な「fork & improve」を禁止
3. **商用制限**: 企業での利用を制限し、普及を妨げる可能性

**代替案**:
- **公開前**: 非公開リポジトリで開発
- **論文発表後**: MIT, Apache 2.0, GPL v3 などで公開
- **引用要請**: LICENSE とは別に CITATION.cff や README で引用を求める

### 6.3 推奨アプローチ

本プロジェクトの目的（オープンな共同開発、アウトリーチ）には：

1. **コアソルバー**: MIT または Apache 2.0
   - フォーク・プルリクエストを最大限促進
   - 企業・研究機関での利用障壁を最小化

2. **引用の確保**:
   - CITATION.cff ファイルで引用情報を明示
   - README に引用のお願いを記載
   - 論文発表時に DOI を付与

3. **コア研究部分**:
   - 論文発表前は別リポジトリ（プライベート）で管理
   - 発表後にマージまたは公開

---

## 7. 参考文献

### 基礎理論
- Hunter, R. J. (2001). *Foundations of Colloid Science* (2nd ed.). Oxford University Press.
- Probstein, R. F. (2003). *Physicochemical Hydrodynamics* (2nd ed.). Wiley.

### 電気浸透流の理論
- Helmholtz, H. (1879). Studien über electrische Grenzschichten. *Annalen der Physik*, 243(7), 337-382.
- Smoluchowski, M. (1903). Contribution à la théorie de l'endosmose électrique et de quelques phénomènes corrélatifs. *Bulletin International de l'Académie des Sciences de Cracovie*, 8, 182-200.
- Rice, C. L., & Whitehead, R. (1965). Electrokinetic flow in a narrow cylindrical capillary. *The Journal of Physical Chemistry*, 69(11), 4017-4024.

### 数値解法
- Patankar, N. A., & Hu, H. H. (1998). Numerical simulation of electroosmotic flow. *Analytical Chemistry*, 70(9), 1870-1881.
- Park, H. M., & Lee, W. M. (2008). Helmholtz–Smoluchowski velocity for viscoelastic electroosmotic flows. *Journal of Colloid and Interface Science*, 317(2), 631-636.

### ベンチマーク問題
- Squires, T. M., & Bazant, M. Z. (2004). Induced-charge electro-osmosis. *Journal of Fluid Mechanics*, 509, 217-252.
- Bazant, M. Z., et al. (2004). Diffuse-charge dynamics in electrochemical systems. *Physical Review E*, 70(2), 021506.
- Rubinstein, I., & Zaltzman, B. (2000). Electro-osmotically induced convection at a permselective membrane. *Physical Review E*, 62(2), 2238.

### イオン液体
- Bazant, M. Z., et al. (2011). Double layer in ionic liquids: Overscreening versus crowding. *Physical Review Letters*, 106(4), 046102.
- Fedorov, M. V., & Kornyshev, A. A. (2014). Ionic liquids at electrified interfaces. *Chemical Reviews*, 114(5), 2978-3036.

### OpenFOAMソルバー
- Roache, A., et al. (2018). Numerical simulation of electrically-driven flows using OpenFOAM. *arXiv:1802.02843*.
- Pimenta, F., & Alves, M. A. (2017). Stabilization of an open-source finite-volume solver for viscoelastic fluid flows. *Journal of Non-Newtonian Fluid Mechanics*, 239, 85-104.
- Schaurer, C., et al. (2023). Electrochemical transport modelling and open-source simulation of pore-scale solid–liquid systems. *Engineering with Computers*.

---

## 付録: 次のステップ

本ソルバーを電気浸透流シミュレーションに拡張する場合の検討事項：

1. **Poisson方程式ソルバーの追加**: 電位場の計算
2. **Nernst-Planck方程式ソルバーの追加**: イオン濃度場の計算
3. **電気体積力項の追加**: Navier-Stokes方程式への $\rho_e \mathbf{E}$ 項
4. **ベンチマーク検証**: Helmholtz-Smoluchowski速度、解析解との比較
