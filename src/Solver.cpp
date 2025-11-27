#include "Solver.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace fluid {

Solver::Solver(double rho, double nu, double dt)
    : rho(rho), nu(nu), dt(dt), cfl(0.5), autoTimeStep(true)
{
}

double Solver::computeTimeStep(const Grid& grid) const {
    // CFL条件から時間刻みを計算
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
    double dt_conv = 1e10;
    if (maxU > 1e-10) dt_conv = std::min(dt_conv, cfl * grid.dx / maxU);
    if (maxV > 1e-10) dt_conv = std::min(dt_conv, cfl * grid.dy / maxV);

    // 拡散CFL条件
    double dt_diff = cfl * 0.5 / (nu * (1.0 / (grid.dx * grid.dx) + 1.0 / (grid.dy * grid.dy)));

    return std::min(dt_conv, dt_diff);
}

double Solver::convectionU(const Grid& grid, int i, int j) const {
    // u速度の移流項: (u・∇)u at u[i][j]
    double dx = grid.dx;
    double dy = grid.dy;

    // u速度の補間
    double u_here = grid.u[i][j];

    // 風上差分
    double u_e = 0.5 * (grid.u[i][j] + grid.u[i + 1][j]);
    double u_w = 0.5 * (grid.u[i - 1][j] + grid.u[i][j]);

    double dudx;
    if (u_here >= 0) {
        dudx = (u_here - grid.u[i - 1][j]) / dx;
    } else {
        dudx = (grid.u[i + 1][j] - u_here) / dx;
    }

    // v速度をu点に補間
    double v_here = 0.25 * (grid.v[i][j - 1] + grid.v[i + 1][j - 1] +
                           grid.v[i][j] + grid.v[i + 1][j]);

    double dudy;
    if (v_here >= 0) {
        dudy = (u_here - grid.u[i][j - 1]) / dy;
    } else {
        dudy = (grid.u[i][j + 1] - u_here) / dy;
    }

    return u_here * dudx + v_here * dudy;
}

double Solver::convectionV(const Grid& grid, int i, int j) const {
    // v速度の移流項: (u・∇)v at v[i][j]
    double dx = grid.dx;
    double dy = grid.dy;

    double v_here = grid.v[i][j];

    // u速度をv点に補間
    double u_here = 0.25 * (grid.u[i - 1][j] + grid.u[i][j] +
                           grid.u[i - 1][j + 1] + grid.u[i][j + 1]);

    double dvdx;
    if (u_here >= 0) {
        dvdx = (v_here - grid.v[i - 1][j]) / dx;
    } else {
        dvdx = (grid.v[i + 1][j] - v_here) / dx;
    }

    double dvdy;
    if (v_here >= 0) {
        dvdy = (v_here - grid.v[i][j - 1]) / dy;
    } else {
        dvdy = (grid.v[i][j + 1] - v_here) / dy;
    }

    return u_here * dvdx + v_here * dvdy;
}

double Solver::diffusionU(const Grid& grid, int i, int j) const {
    double dx2 = grid.dx * grid.dx;
    double dy2 = grid.dy * grid.dy;

    return (grid.u[i + 1][j] - 2.0 * grid.u[i][j] + grid.u[i - 1][j]) / dx2 +
           (grid.u[i][j + 1] - 2.0 * grid.u[i][j] + grid.u[i][j - 1]) / dy2;
}

double Solver::diffusionV(const Grid& grid, int i, int j) const {
    double dx2 = grid.dx * grid.dx;
    double dy2 = grid.dy * grid.dy;

    return (grid.v[i + 1][j] - 2.0 * grid.v[i][j] + grid.v[i - 1][j]) / dx2 +
           (grid.v[i][j + 1] - 2.0 * grid.v[i][j] + grid.v[i][j - 1]) / dy2;
}

void Solver::computeIntermediateVelocity(Grid& grid, const BoundaryCondition& bc) {
    int nx = grid.nx;
    int ny = grid.ny;

    // u* の計算
    #ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int i = 1; i < nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            double conv = convectionU(grid, i, j);
            double diff = diffusionU(grid, i, j);
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
            grid.v_star[i][j] = grid.v[i][j] + dt * (-conv + nu * diff);
        }
    }

    // 境界条件を中間速度場にも適用
    // 左右境界のu_star
    for (int j = 0; j < ny + 2; ++j) {
        grid.u_star[0][j] = grid.u[0][j];
        grid.u_star[nx][j] = grid.u[nx][j];
    }

    // 上下境界のv_star
    for (int i = 0; i < nx + 2; ++i) {
        grid.v_star[i][0] = grid.v[i][0];
        grid.v_star[i][ny] = grid.v[i][ny];
    }
}

void Solver::correctVelocity(Grid& grid) {
    int nx = grid.nx;
    int ny = grid.ny;
    double dx = grid.dx;
    double dy = grid.dy;

    // u^{n+1} = u* - (dt/ρ) * ∂p/∂x
    #ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int i = 1; i < nx; ++i) {
        for (int j = 1; j <= ny; ++j) {
            grid.u[i][j] = grid.u_star[i][j] - (dt / rho) * (grid.p[i + 1][j] - grid.p[i][j]) / dx;
        }
    }

    // v^{n+1} = v* - (dt/ρ) * ∂p/∂y
    #ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
    #endif
    for (int i = 1; i <= nx; ++i) {
        for (int j = 1; j < ny; ++j) {
            grid.v[i][j] = grid.v_star[i][j] - (dt / rho) * (grid.p[i][j + 1] - grid.p[i][j]) / dy;
        }
    }
}

int Solver::step(Grid& grid, const BoundaryCondition& bc) {
    // 適応時間刻み
    if (autoTimeStep) {
        dt = computeTimeStep(grid);
    }

    // Step 1: 中間速度場の計算
    computeIntermediateVelocity(grid, bc);

    // Step 2: 圧力ポアソン方程式を解く
    int pressureIter = pressureSolver.solve(grid, bc, dt, rho);

    // Step 3: 速度の補正
    correctVelocity(grid);

    // 境界条件を適用
    bc.apply(grid);

    time += dt;
    stepCount++;

    return pressureIter;
}

void Solver::run(Grid& grid, const BoundaryCondition& bc, double endTime,
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
