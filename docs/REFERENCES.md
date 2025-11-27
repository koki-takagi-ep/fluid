# References

## 数値解法

### Projection法（分離法）

- Chorin, A. J. (1968). "Numerical solution of the Navier-Stokes equations." *Mathematics of Computation*, 22(104), 745-762.
  - Projection法の原論文

- Temam, R. (1969). "Sur l'approximation de la solution des équations de Navier-Stokes par la méthode des pas fractionnaires." *Archive for Rational Mechanics and Analysis*, 33(5), 377-385.

### MAC法（Marker-and-Cell）

- Harlow, F. H., & Welch, J. E. (1965). "Numerical calculation of time-dependent viscous incompressible flow of fluid with free surface." *Physics of Fluids*, 8(12), 2182-2189.
  - スタガード格子の原論文

### SOR法

- Young, D. M. (1954). "Iterative methods for solving partial difference equations of elliptic type." *Transactions of the American Mathematical Society*, 76(1), 92-111.

## ベンチマーク

### Lid-driven cavity flow

- **Ghia, U., Ghia, K. N., & Shin, C. T. (1982).** "High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method." *Journal of Computational Physics*, 48(3), 387-411.
  - 本コードの検証に使用したベンチマークデータ
  - Re = 100, 400, 1000, 3200, 5000, 7500, 10000 での中心線速度分布

- Botella, O., & Peyret, R. (1998). "Benchmark spectral results on the lid-driven cavity flow." *Computers & Fluids*, 27(4), 421-433.
  - 高精度スペクトル法による参照解

### Channel flow (Poiseuille flow)

- White, F. M. (2006). *Viscous Fluid Flow* (3rd ed.). McGraw-Hill.
  - 第3章: Poiseuille流れの理論解

## 教科書

### 数値流体力学

- Ferziger, J. H., Perić, M., & Street, R. L. (2020). *Computational Methods for Fluid Dynamics* (4th ed.). Springer.
  - 第7章: 非圧縮性流れの解法

- Versteeg, H. K., & Malalasekera, W. (2007). *An Introduction to Computational Fluid Dynamics: The Finite Volume Method* (2nd ed.). Pearson.

### 流体力学

- Kundu, P. K., Cohen, I. M., & Dowling, D. R. (2015). *Fluid Mechanics* (6th ed.). Academic Press.

- Batchelor, G. K. (2000). *An Introduction to Fluid Dynamics*. Cambridge University Press.

## オンラインリソース

- CFD-Online: https://www.cfd-online.com/
  - CFDコミュニティのフォーラム・リソース

- NASA CFD Tutorial: https://www.grc.nasa.gov/www/wind/valid/tutorial/tutorial.html
  - CFD検証・妥当性確認のチュートリアル
