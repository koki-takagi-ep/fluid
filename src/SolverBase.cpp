/**
 * @file SolverBase.cpp
 * @brief SolverBase クラスの実装
 */

#include "SolverBase.hpp"
#include "Constants.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace fluid {

SolverBase::SolverBase(double rho, double nu, double dt)
    : rho(rho), nu(nu), dt(dt), cfl(constants::DEFAULT_CFL), autoTimeStep(true),
      limiterType(LimiterType::None)
{
}

double SolverBase::computeTimeStep(const Grid& grid) const {
    // 最大速度を計算
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

double SolverBase::tvdConvection1D(double vel, double phi_mm, double phi_m,
                                    double phi_c, double phi_p, double dx) const {
    // TVDスキームによる移流項計算
    //
    // Flux Corrected Transport (FCT) / TVD形式:
    //   F = F_low + ψ(r) * (F_high - F_low)
    //
    // F_low  = 1次風上差分（安定だが拡散的）
    // F_high = 2次中心差分（高精度だが振動）
    // ψ(r)   = リミッター関数（0〜1の範囲）

    // 1次風上差分
    double upwind;
    if (vel >= 0.0) {
        upwind = vel * (phi_c - phi_m) / dx;
    } else {
        upwind = vel * (phi_p - phi_c) / dx;
    }

    if (limiterType == LimiterType::None) {
        return upwind;
    }

    // 2次中心差分
    double central = vel * (phi_p - phi_m) / (2.0 * dx);

    // antidiffusive flux = central - upwind
    double adf = central - upwind;

    // 勾配比を計算してリミッターを適用
    const double eps = 1.0e-10;
    double r;

    if (vel >= 0.0) {
        // 正の速度: 上流側（左）の勾配比
        // r = (φ_c - φ_m) / (φ_p - φ_c) を使う場合もあるが、
        // TVD条件では上流の連続勾配比が標準的
        double delta_c = phi_c - phi_m;   // 現在点での後方差分
        double delta_d = phi_p - phi_c;   // 現在点での前方差分

        r = (std::abs(delta_d) > eps) ? delta_c / delta_d : 1.0;
    } else {
        // 負の速度: 上流側（右）の勾配比
        double delta_c = phi_c - phi_p;   // 現在点での後方差分（右から見て）
        double delta_d = phi_m - phi_c;   // 現在点での前方差分（右から見て）

        r = (std::abs(delta_d) > eps) ? delta_c / delta_d : 1.0;
    }

    double psi = limiter::apply(limiterType, r);

    // TVDスキーム: 1次風上 + 制限された高次補正
    return upwind + psi * adf;
}

double SolverBase::convectionU(const Grid& grid, int i, int j) const {
    // u速度の移流項: (u・∇)u at u[i][j]
    const double dx = grid.dx;
    const double dy = grid.dy;
    const double u_here = grid.u[i][j];

    // v速度をu点に補間（4点平均）
    const double v_here = constants::FOUR_POINT_AVERAGE_COEFF
        * (grid.v[i][j - 1] + grid.v[i + 1][j - 1] + grid.v[i][j] + grid.v[i + 1][j]);

    double dudx, dudy;

    if (limiterType == LimiterType::None) {
        // 1次精度風上差分スキーム
        dudx = (u_here >= 0)
            ? (u_here - grid.u[i - 1][j]) / dx
            : (grid.u[i + 1][j] - u_here) / dx;

        dudy = (v_here >= 0)
            ? (u_here - grid.u[i][j - 1]) / dy
            : (grid.u[i][j + 1] - u_here) / dy;

        return u_here * dudx + v_here * dudy;
    }

    // TVDスキーム
    // x方向: 境界チェック付きでステンシルを取得
    int nx = grid.nx;
    double u_mm = (i >= 2) ? grid.u[i - 2][j] : grid.u[i - 1][j];
    double u_m  = grid.u[i - 1][j];
    double u_c  = grid.u[i][j];
    double u_p  = (i < nx - 1) ? grid.u[i + 1][j] : grid.u[i][j];
    double u_pp = (i < nx - 2) ? grid.u[i + 2][j] : u_p;

    // y方向
    int ny = grid.ny;
    double u_jmm = (j >= 2) ? grid.u[i][j - 2] : grid.u[i][j - 1];
    double u_jm  = grid.u[i][j - 1];
    double u_jc  = grid.u[i][j];
    double u_jp  = (j < ny) ? grid.u[i][j + 1] : grid.u[i][j];
    double u_jpp = (j < ny - 1) ? grid.u[i][j + 2] : u_jp;

    // TVD移流項: 速度方向に関係なく同じステンシルを渡す
    // tvdConvection1D内で速度の符号に基づいて処理を切り替え
    double conv_x = tvdConvection1D(u_here, u_mm, u_m, u_c, u_p, dx);
    double conv_y = tvdConvection1D(v_here, u_jmm, u_jm, u_jc, u_jp, dy);

    return conv_x + conv_y;
}

double SolverBase::convectionV(const Grid& grid, int i, int j) const {
    // v速度の移流項: (u・∇)v at v[i][j]
    const double dx = grid.dx;
    const double dy = grid.dy;
    const double v_here = grid.v[i][j];

    // u速度をv点に補間（4点平均）
    const double u_here = constants::FOUR_POINT_AVERAGE_COEFF
        * (grid.u[i - 1][j] + grid.u[i][j] + grid.u[i - 1][j + 1] + grid.u[i][j + 1]);

    double dvdx, dvdy;

    if (limiterType == LimiterType::None) {
        // 1次精度風上差分スキーム
        dvdx = (u_here >= 0)
            ? (v_here - grid.v[i - 1][j]) / dx
            : (grid.v[i + 1][j] - v_here) / dx;

        dvdy = (v_here >= 0)
            ? (v_here - grid.v[i][j - 1]) / dy
            : (grid.v[i][j + 1] - v_here) / dy;

        return u_here * dvdx + v_here * dvdy;
    }

    // TVDスキーム
    // x方向: 境界チェック付きでステンシルを取得
    int nx = grid.nx;
    double v_mm = (i >= 2) ? grid.v[i - 2][j] : grid.v[i - 1][j];
    double v_m  = grid.v[i - 1][j];
    double v_c  = grid.v[i][j];
    double v_p  = (i < nx) ? grid.v[i + 1][j] : grid.v[i][j];
    double v_pp = (i < nx - 1) ? grid.v[i + 2][j] : v_p;

    // y方向
    int ny = grid.ny;
    double v_jmm = (j >= 2) ? grid.v[i][j - 2] : grid.v[i][j - 1];
    double v_jm  = grid.v[i][j - 1];
    double v_jc  = grid.v[i][j];
    double v_jp  = (j < ny - 1) ? grid.v[i][j + 1] : grid.v[i][j];
    double v_jpp = (j < ny - 2) ? grid.v[i][j + 2] : v_jp;

    // TVD移流項: 速度方向に関係なく同じステンシルを渡す
    // tvdConvection1D内で速度の符号に基づいて処理を切り替え
    double conv_x = tvdConvection1D(u_here, v_mm, v_m, v_c, v_p, dx);
    double conv_y = tvdConvection1D(v_here, v_jmm, v_jm, v_jc, v_jp, dy);

    return conv_x + conv_y;
}

double SolverBase::diffusionU(const Grid& grid, int i, int j) const {
    // 2次精度中心差分によるラプラシアン: ∇²u
    const double dx2 = grid.dx * grid.dx;
    const double dy2 = grid.dy * grid.dy;
    constexpr double c = constants::LAPLACIAN_CENTER_COEFF;

    return (grid.u[i + 1][j] - c * grid.u[i][j] + grid.u[i - 1][j]) / dx2 +
           (grid.u[i][j + 1] - c * grid.u[i][j] + grid.u[i][j - 1]) / dy2;
}

double SolverBase::diffusionV(const Grid& grid, int i, int j) const {
    // 2次精度中心差分によるラプラシアン: ∇²v
    const double dx2 = grid.dx * grid.dx;
    const double dy2 = grid.dy * grid.dy;
    constexpr double c = constants::LAPLACIAN_CENTER_COEFF;

    return (grid.v[i + 1][j] - c * grid.v[i][j] + grid.v[i - 1][j]) / dx2 +
           (grid.v[i][j + 1] - c * grid.v[i][j] + grid.v[i][j - 1]) / dy2;
}

void SolverBase::computeIntermediateVelocity(Grid& grid, const BoundaryCondition& /*bc*/) {
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

void SolverBase::correctVelocity(Grid& grid) {
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

void SolverBase::run(Grid& grid, const BoundaryCondition& bc, double endTime,
                     int outputInterval,
                     std::function<void(int, double, const Grid&)> outputCallback) {

    // 初期境界条件
    bc.apply(grid);

    if (outputCallback) {
        outputCallback(0, 0.0, grid);
    }

    while (time_ < endTime) {
        int iterCount = step(grid, bc);

        if (stepCount_ % outputInterval == 0) {
            double maxDiv = grid.maxDivergence();
            std::cout << "Step: " << stepCount_
                      << ", Time: " << time_
                      << ", dt: " << dt
                      << ", Pressure iter: " << iterCount
                      << ", Max divergence: " << maxDiv
                      << std::endl;

            if (outputCallback) {
                outputCallback(stepCount_, time_, grid);
            }
        }
    }

    // 最終結果を出力
    if (outputCallback) {
        outputCallback(stepCount_, time_, grid);
    }
}

} // namespace fluid
