/**
 * @file SimpleSolver.cpp
 * @brief SIMPLE法の実装
 */

#include "SimpleSolver.hpp"
#include "Constants.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace fluid {

SimpleSolver::SimpleSolver(double rho, double nu, double dt)
    : rho(rho), nu(nu), dt(dt),
      cfl(constants::DEFAULT_CFL),
      autoTimeStep(true),
      alpha_u(0.7),           // 速度緩和係数
      alpha_p(0.3),           // 圧力緩和係数
      omega(constants::DEFAULT_SOR_OMEGA),
      maxPressureIter(constants::DEFAULT_PRESSURE_MAX_ITERATIONS),
      pressureTol(constants::DEFAULT_PRESSURE_TOLERANCE),
      maxOuterIter(100),      // SIMPLE外部反復
      convergenceTol(1e-4)
{
}

double SimpleSolver::computeTimeStep(const Grid& grid) const {
    double maxU = 0.0, maxV = 0.0;

    for (int i = 0; i < grid.nx + 1; ++i) {
        for (int j = 0; j < grid.ny + 2; ++j) {
            maxU = std::max(maxU, std::abs(grid.u[i][j]));
        }
    }

    for (int i = 0; i < grid.nx + 2; ++i) {
        for (int j = 0; j < grid.ny + 1; ++j) {
            maxV = std::max(maxV, std::abs(grid.v[i][j]));
        }
    }

    // 移流CFL条件
    double dt_conv = constants::TIMESTEP_MAX_INITIAL;
    if (maxU > constants::VELOCITY_ZERO_THRESHOLD) {
        dt_conv = std::min(dt_conv, cfl * grid.dx / maxU);
    }
    if (maxV > constants::VELOCITY_ZERO_THRESHOLD) {
        dt_conv = std::min(dt_conv, cfl * grid.dy / maxV);
    }

    // 拡散CFL条件
    double dt_diff = cfl * constants::VISCOUS_CFL_FACTOR
                   / (nu * (1.0 / (grid.dx * grid.dx) + 1.0 / (grid.dy * grid.dy)));

    return std::min(dt_conv, dt_diff);
}

double SimpleSolver::convectionU(const Grid& grid, int i, int j) const {
    const double dx = grid.dx;
    const double dy = grid.dy;
    const double u_here = grid.u[i][j];

    // ∂u/∂x: 風上差分
    const double dudx = (u_here >= 0)
        ? (u_here - grid.u[i - 1][j]) / dx
        : (grid.u[i + 1][j] - u_here) / dx;

    // v速度をu点に補間
    const double v_here = constants::FOUR_POINT_AVERAGE_COEFF
        * (grid.v[i][j - 1] + grid.v[i + 1][j - 1] + grid.v[i][j] + grid.v[i + 1][j]);

    // ∂u/∂y: 風上差分
    const double dudy = (v_here >= 0)
        ? (u_here - grid.u[i][j - 1]) / dy
        : (grid.u[i][j + 1] - u_here) / dy;

    return u_here * dudx + v_here * dudy;
}

double SimpleSolver::convectionV(const Grid& grid, int i, int j) const {
    const double dx = grid.dx;
    const double dy = grid.dy;
    const double v_here = grid.v[i][j];

    // u速度をv点に補間
    const double u_here = constants::FOUR_POINT_AVERAGE_COEFF
        * (grid.u[i - 1][j] + grid.u[i][j] + grid.u[i - 1][j + 1] + grid.u[i][j + 1]);

    // ∂v/∂x: 風上差分
    const double dvdx = (u_here >= 0)
        ? (v_here - grid.v[i - 1][j]) / dx
        : (grid.v[i + 1][j] - v_here) / dx;

    // ∂v/∂y: 風上差分
    const double dvdy = (v_here >= 0)
        ? (v_here - grid.v[i][j - 1]) / dy
        : (grid.v[i][j + 1] - v_here) / dy;

    return u_here * dvdx + v_here * dvdy;
}

double SimpleSolver::diffusionU(const Grid& grid, int i, int j) const {
    const double dx2 = grid.dx * grid.dx;
    const double dy2 = grid.dy * grid.dy;
    constexpr double c = constants::LAPLACIAN_CENTER_COEFF;

    return (grid.u[i + 1][j] - c * grid.u[i][j] + grid.u[i - 1][j]) / dx2 +
           (grid.u[i][j + 1] - c * grid.u[i][j] + grid.u[i][j - 1]) / dy2;
}

double SimpleSolver::diffusionV(const Grid& grid, int i, int j) const {
    const double dx2 = grid.dx * grid.dx;
    const double dy2 = grid.dy * grid.dy;
    constexpr double c = constants::LAPLACIAN_CENTER_COEFF;

    return (grid.v[i + 1][j] - c * grid.v[i][j] + grid.v[i - 1][j]) / dx2 +
           (grid.v[i][j + 1] - c * grid.v[i][j] + grid.v[i][j - 1]) / dy2;
}

void SimpleSolver::computeStarVelocity(Grid& grid, const BoundaryCondition& /*bc*/) {
    int nx = grid.nx;
    int ny = grid.ny;

    // u* の計算 (SIMPLE法: 圧力勾配なしで予測)
    // Projection法と同様のアプローチ
    #ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int i = 1; i < nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            double conv = convectionU(grid, i, j);
            double diff = diffusionU(grid, i, j);

            // 中間速度（圧力勾配なし）
            grid.u_star[i][j] = grid.u[i][j] + dt * (-conv + nu * diff);
        }
    }

    // v* の計算
    #ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j < ny; ++j) {
            double conv = convectionV(grid, i, j);
            double diff = diffusionV(grid, i, j);

            // 中間速度（圧力勾配なし）
            grid.v_star[i][j] = grid.v[i][j] + dt * (-conv + nu * diff);
        }
    }

    // 境界条件を中間速度場にも適用
    for (int j = 0; j < ny + 2; ++j) {
        grid.u_star[0][j] = grid.u[0][j];
        grid.u_star[nx][j] = grid.u[nx][j];
    }

    for (int i = 0; i < nx + 2; ++i) {
        grid.v_star[i][0] = grid.v[i][0];
        grid.v_star[i][ny] = grid.v[i][ny];
    }
}

int SimpleSolver::solvePressureCorrection(Grid& grid, const BoundaryCondition& bc) {
    int nx = grid.nx;
    int ny = grid.ny;
    double dx = grid.dx;
    double dy = grid.dy;

    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double factor = constants::LAPLACIAN_CENTER_COEFF * (1.0 / dx2 + 1.0 / dy2);

    // 圧力補正場を初期化（現在の圧力をコピー）
    if (p_prime.empty() || static_cast<int>(p_prime.size()) != nx + 2) {
        p_prime.resize(nx + 2, std::vector<double>(ny + 2, 0.0));
    }

    // 初期値として現在の圧力場をコピー
    for (int i = 0; i < nx + 2; ++i) {
        for (int j = 0; j < ny + 2; ++j) {
            p_prime[i][j] = grid.p[i][j];
        }
    }

    // SOR反復で圧力Poisson方程式を解く
    // ∇²p = (ρ/Δt) * ∇・u*
    for (int iter = 0; iter < maxPressureIter; ++iter) {
        double maxResidual = 0.0;

        // Red-Black SOR
        for (int color = 0; color < 2; ++color) {
            #ifdef USE_OPENMP
            #pragma omp parallel for reduction(max:maxResidual)
            #endif
            for (int i = 1; i <= nx; ++i) {
                for (int j = 1; j <= ny; ++j) {
                    if ((i + j) % 2 != color) continue;

                    // 右辺: ∇・u*
                    double dudx = (grid.u_star[i][j] - grid.u_star[i - 1][j]) / dx;
                    double dvdy = (grid.v_star[i][j] - grid.v_star[i][j - 1]) / dy;
                    double rhs = rho * (dudx + dvdy) / dt;

                    double p_new = (
                        (p_prime[i + 1][j] + p_prime[i - 1][j]) / dx2 +
                        (p_prime[i][j + 1] + p_prime[i][j - 1]) / dy2 -
                        rhs
                    ) / factor;

                    double residual = std::abs(p_new - p_prime[i][j]);
                    maxResidual = std::max(maxResidual, residual);

                    p_prime[i][j] = p_prime[i][j] + omega * (p_new - p_prime[i][j]);
                }
            }
        }

        // 境界条件を適用（元のapplyPressureBCと同様）
        // 左境界: Neumann
        for (int j = 1; j <= ny; ++j) {
            p_prime[0][j] = p_prime[1][j];
        }
        // 右境界: Dirichlet（流出境界でゲージ圧=0）
        if (bc.right == BCType::Outflow) {
            for (int j = 1; j <= ny; ++j) {
                p_prime[nx + 1][j] = 0.0;
            }
        } else {
            for (int j = 1; j <= ny; ++j) {
                p_prime[nx + 1][j] = p_prime[nx][j];
            }
        }
        // 下境界: Neumann
        for (int i = 1; i <= nx; ++i) {
            p_prime[i][0] = p_prime[i][1];
        }
        // 上境界: Neumann
        for (int i = 1; i <= nx; ++i) {
            p_prime[i][ny + 1] = p_prime[i][ny];
        }
        // 角
        p_prime[0][0] = 0.5 * (p_prime[1][0] + p_prime[0][1]);
        p_prime[nx + 1][0] = 0.5 * (p_prime[nx][0] + p_prime[nx + 1][1]);
        p_prime[0][ny + 1] = 0.5 * (p_prime[1][ny + 1] + p_prime[0][ny]);
        p_prime[nx + 1][ny + 1] = 0.5 * (p_prime[nx][ny + 1] + p_prime[nx + 1][ny]);

        if (maxResidual < pressureTol) {
            return iter + 1;
        }
    }

    return -1;  // 収束しなかった
}

void SimpleSolver::correctVelocityAndPressure(Grid& grid) {
    int nx = grid.nx;
    int ny = grid.ny;
    double dx = grid.dx;
    double dy = grid.dy;

    // 圧力の更新: p = p' (Projection法と同様、p'が圧力そのもの)
    #ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            grid.p[i][j] = p_prime[i][j];
        }
    }

    // 速度の補正: u = u* - (Δt/ρ) * ∂p/∂x
    #ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int i = 1; i < nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            grid.u[i][j] = grid.u_star[i][j]
                         - (dt / rho) * (grid.p[i + 1][j] - grid.p[i][j]) / dx;
        }
    }

    // v = v* - (Δt/ρ) * ∂p/∂y
    #ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j < ny; ++j) {
            grid.v[i][j] = grid.v_star[i][j]
                         - (dt / rho) * (grid.p[i][j + 1] - grid.p[i][j]) / dy;
        }
    }
}

double SimpleSolver::computeResidual(const Grid& grid) const {
    // 連続の式の残差（質量保存）
    double maxDiv = 0.0;
    for (int i = 0; i < grid.nx; ++i) {
        for (int j = 0; j < grid.ny; ++j) {
            maxDiv = std::max(maxDiv, std::abs(grid.computeDivergence(i, j)));
        }
    }
    return maxDiv;
}

int SimpleSolver::step(Grid& grid, const BoundaryCondition& bc) {
    // 適応時間刻み
    if (autoTimeStep) {
        dt = computeTimeStep(grid);
    }

    // Step 1: 推定速度場の計算（圧力勾配なし）
    computeStarVelocity(grid, bc);

    // Step 2: 圧力Poisson方程式を解く
    int pressureIter = solvePressureCorrection(grid, bc);

    // Step 3: 速度と圧力の補正
    correctVelocityAndPressure(grid);

    // 境界条件を適用
    bc.apply(grid);

    lastOuterIter = pressureIter;
    time += dt;
    stepCount++;

    return pressureIter;
}

void SimpleSolver::run(Grid& grid, const BoundaryCondition& bc, double endTime,
                       int outputInterval,
                       std::function<void(int, double, const Grid&)> outputCallback) {

    // 初期境界条件
    bc.apply(grid);

    if (outputCallback) {
        outputCallback(0, 0.0, grid);
    }

    while (time < endTime) {
        int pressureIter = step(grid, bc);

        if (stepCount % outputInterval == 0) {
            double maxDiv = grid.maxDivergence();
            std::cout << "Step: " << stepCount
                      << ", Time: " << time
                      << ", dt: " << dt
                      << ", Pressure iter: " << pressureIter
                      << ", Max divergence: " << maxDiv
                      << std::endl;

            if (outputCallback) {
                outputCallback(stepCount, time, grid);
            }
        }
    }

    // 最終結果を出力
    if (outputCallback) {
        outputCallback(stepCount, time, grid);
    }
}

} // namespace fluid
