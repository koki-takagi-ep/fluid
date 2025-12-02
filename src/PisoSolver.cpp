/**
 * @file PisoSolver.cpp
 * @brief PISO法ソルバーの実装
 */

#include "PisoSolver.hpp"
#include "Constants.hpp"
#include <cmath>
#include <algorithm>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace fluid {

PisoSolver::PisoSolver(double rho, double nu, double dt, int nCorrectors)
    : SolverBase(rho, nu, dt),
      omega(constants::DEFAULT_SOR_OMEGA),
      maxPressureIter(constants::DEFAULT_PRESSURE_MAX_ITERATIONS),
      pressureTol(constants::DEFAULT_PRESSURE_TOLERANCE),
      nCorrectors(nCorrectors)
{
}

void PisoSolver::ensureArraySize(std::vector<std::vector<double>>& arr, int rows, int cols) {
    if (arr.empty() || static_cast<int>(arr.size()) != rows ||
        static_cast<int>(arr[0].size()) != cols) {
        arr.resize(rows, std::vector<double>(cols, 0.0));
    }
}

int PisoSolver::solvePressureCorrection(Grid& grid, const BoundaryCondition& bc,
                                         const std::vector<std::vector<double>>& u_field,
                                         const std::vector<std::vector<double>>& v_field) {
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
            double dudx = (u_field[i][j] - u_field[i - 1][j]) / dx;
            double dvdy = (v_field[i][j] - v_field[i][j - 1]) / dy;
            rhs_array[i][j] = rho * (dudx + dvdy) / dt;
        }
    }

    // 非定常PISO法: 前回のタイムステップの圧力場を初期値として使用
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

int PisoSolver::step(Grid& grid, const BoundaryCondition& bc) {
    int nx = grid.nx;
    int ny = grid.ny;

    // 適応時間刻み
    if (autoTimeStep) {
        dt = computeTimeStep(grid);
    }

    // 配列のサイズを確保
    ensureArraySize(u_double_star_, nx + 1, ny + 2);
    ensureArraySize(v_double_star_, nx + 2, ny + 1);

    // =========================================================================
    // Predictor: 中間速度場 u* を計算
    // =========================================================================
    computeIntermediateVelocity(grid, bc);

    int totalPressureIter = 0;

    // =========================================================================
    // First Corrector: 圧力Poisson方程式を解き、速度を補正
    // =========================================================================
    // 圧力はsolvePressureCorrectionで直接grid.pに更新される
    int iter1 = solvePressureCorrection(grid, bc, grid.u_star, grid.v_star);
    if (iter1 > 0) {
        totalPressureIter += iter1;
    } else {
        totalPressureIter += maxPressureIter;
    }

    // 速度を補正: u** = u* - (dt/ρ) * ∇p
    correctVelocityWithPressure(grid, grid.u_star, grid.v_star);

    // 補正後の速度をu_double_star_にコピー
    for (int i = 0; i < nx + 1; ++i) {
        for (int j = 0; j < ny + 2; ++j) {
            u_double_star_[i][j] = grid.u[i][j];
        }
    }
    for (int i = 0; i < nx + 2; ++i) {
        for (int j = 0; j < ny + 1; ++j) {
            v_double_star_[i][j] = grid.v[i][j];
        }
    }

    // 追加の補正ステップ（nCorrectors >= 2 の場合）
    for (int corrector = 2; corrector <= nCorrectors; ++corrector) {
        // =========================================================================
        // Additional Corrector: 圧力を再度解き、速度を再補正
        // =========================================================================
        int iterN = solvePressureCorrection(grid, bc, u_double_star_, v_double_star_);
        if (iterN > 0) {
            totalPressureIter += iterN;
        } else {
            totalPressureIter += maxPressureIter;
        }

        // 速度を補正: u = u** - (dt/ρ) * ∇p
        correctVelocityWithPressure(grid, u_double_star_, v_double_star_);

        // 次の補正ステップのためにu_double_star_を更新
        if (corrector < nCorrectors) {
            for (int i = 0; i < nx + 1; ++i) {
                for (int j = 0; j < ny + 2; ++j) {
                    u_double_star_[i][j] = grid.u[i][j];
                }
            }
            for (int i = 0; i < nx + 2; ++i) {
                for (int j = 0; j < ny + 1; ++j) {
                    v_double_star_[i][j] = grid.v[i][j];
                }
            }
        }
    }

    // 境界条件を適用
    bc.apply(grid);

    lastPressureIter_ = totalPressureIter;
    time_ += dt;
    stepCount_++;

    return totalPressureIter;
}

void PisoSolver::correctVelocityWithPressure(
    Grid& grid,
    const std::vector<std::vector<double>>& u_field,
    const std::vector<std::vector<double>>& v_field) {

    int nx = grid.nx;
    int ny = grid.ny;
    double dx = grid.dx;
    double dy = grid.dy;

    // u^{n+1} = u_field - (dt/ρ) * ∂p/∂x
    #ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int i = 1; i < nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            grid.u[i][j] = u_field[i][j] - (dt / rho) * (grid.p[i + 1][j] - grid.p[i][j]) / dx;
        }
    }

    // v^{n+1} = v_field - (dt/ρ) * ∂p/∂y
    #ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j < ny; ++j) {
            grid.v[i][j] = v_field[i][j] - (dt / rho) * (grid.p[i][j + 1] - grid.p[i][j]) / dy;
        }
    }
}

} // namespace fluid
