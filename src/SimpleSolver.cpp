/**
 * @file SimpleSolver.cpp
 * @brief SIMPLE法ソルバーの実装
 */

#include "SimpleSolver.hpp"
#include "Constants.hpp"
#include <cmath>
#include <algorithm>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace fluid {

SimpleSolver::SimpleSolver(double rho, double nu, double dt)
    : SolverBase(rho, nu, dt),
      alpha_u(0.7),
      alpha_p(0.3),
      omega(constants::DEFAULT_SOR_OMEGA),
      maxPressureIter(constants::DEFAULT_PRESSURE_MAX_ITERATIONS),
      pressureTol(constants::DEFAULT_PRESSURE_TOLERANCE),
      maxOuterIter(100),
      convergenceTol(1e-4)
{
}

int SimpleSolver::solvePressureCorrection(Grid& grid, const BoundaryCondition& bc) {
    int nx = grid.nx;
    int ny = grid.ny;
    double dx = grid.dx;
    double dy = grid.dy;

    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double factor = constants::LAPLACIAN_CENTER_COEFF * (1.0 / dx2 + 1.0 / dy2);

    // RHS（右辺）を事前計算
    std::vector<std::vector<double>> rhs_array(nx + 2, std::vector<double>(ny + 2, 0.0));
    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            double dudx = (grid.u_star[i][j] - grid.u_star[i - 1][j]) / dx;
            double dvdy = (grid.v_star[i][j] - grid.v_star[i][j - 1]) / dy;
            rhs_array[i][j] = rho * (dudx + dvdy) / dt;
        }
    }

    // 非定常SIMPLE法: 前回のタイムステップの圧力場を初期値として使用
    // （Projection法と同様に、grid.pを直接更新する）
    // これにより2次精度の格子収束性を達成できる

    // SOR反復で圧力Poisson方程式を解く
    for (int iter = 0; iter < maxPressureIter; ++iter) {
        // Red-Black SOR
        for (int color = 0; color < 2; ++color) {
            #ifdef USE_OPENMP
            #pragma omp parallel for collapse(2)
            #endif
            for (int i = 1; i <= nx; ++i) {
                for (int j = 1; j <= ny; ++j) {
                    if ((i + j) % 2 != color) continue;

                    double p_new = (
                        (grid.p[i + 1][j] + grid.p[i - 1][j]) / dx2 +
                        (grid.p[i][j + 1] + grid.p[i][j - 1]) / dy2 -
                        rhs_array[i][j]
                    ) / factor;

                    grid.p[i][j] = grid.p[i][j] + omega * (p_new - grid.p[i][j]);
                }
            }
        }

        // 境界条件を適用
        bc.applyPressureBC(grid);

        // Poisson方程式の残差で収束判定
        double maxResidual = 0.0;
        for (int i = 1; i <= nx; ++i) {
            for (int j = 1; j <= ny; ++j) {
                double laplacian_p = (grid.p[i + 1][j] - 2.0 * grid.p[i][j] + grid.p[i - 1][j]) / dx2
                                   + (grid.p[i][j + 1] - 2.0 * grid.p[i][j] + grid.p[i][j - 1]) / dy2;
                double residual = std::abs(laplacian_p - rhs_array[i][j]);
                maxResidual = std::max(maxResidual, residual);
            }
        }

        if (maxResidual < pressureTol * rho / dt) {
            return iter + 1;
        }
    }

    return -1;  // 収束しなかった
}

void SimpleSolver::correctVelocityAndPressure(Grid& grid) {
    // 圧力はsolvePressureCorrectionで既に更新済み
    // 速度の補正（基底クラスのメソッドを使用）
    correctVelocity(grid);
}

int SimpleSolver::step(Grid& grid, const BoundaryCondition& bc) {
    // 適応時間刻み
    if (autoTimeStep) {
        dt = computeTimeStep(grid);
    }

    // Step 1: 推定速度場の計算（基底クラスのメソッド）
    computeIntermediateVelocity(grid, bc);

    // Step 2: 圧力Poisson方程式を解く
    int pressureIter = solvePressureCorrection(grid, bc);

    // Step 3: 速度と圧力の補正
    correctVelocityAndPressure(grid);

    // 境界条件を適用
    bc.apply(grid);

    lastOuterIter_ = pressureIter;
    time_ += dt;
    stepCount_++;

    return pressureIter;
}

} // namespace fluid
