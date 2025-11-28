#include "PressureSolver.hpp"
#include <cmath>
#include <algorithm>

#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace fluid {

PressureSolver::PressureSolver(double omega, int maxIter, double tol)
    : omega(omega), maxIterations(maxIter), tolerance(tol)
{
}

double PressureSolver::computeRHS(const Grid& grid, int i, int j, double dt, double rho) const {
    // 中間速度場の発散を計算
    // i, j は内部セルインデックス (1 <= i <= nx, 1 <= j <= ny)
    double dudx = (grid.u_star[i][j] - grid.u_star[i - 1][j]) / grid.dx;
    double dvdy = (grid.v_star[i][j] - grid.v_star[i][j - 1]) / grid.dy;
    return rho * (dudx + dvdy) / dt;
}

int PressureSolver::solve(Grid& grid, const BoundaryCondition& bc, double dt, double rho) {
    int nx = grid.nx;
    int ny = grid.ny;
    double dx = grid.dx;
    double dy = grid.dy;

    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double factor = 2.0 * (1.0 / dx2 + 1.0 / dy2);

    // 初期圧力を0にリセット（ゲージ圧として計算）
    for (int i = 0; i < nx + 2; ++i) {
        for (int j = 0; j < ny + 2; ++j) {
            grid.p[i][j] = 0.0;
        }
    }

    // SOR反復
    for (int iter = 0; iter < maxIterations; ++iter) {
        double maxResidual = 0.0;

        // Red-Black SOR (並列化可能)
        for (int color = 0; color < 2; ++color) {
            #ifdef USE_OPENMP
            #pragma omp parallel for reduction(max:maxResidual)
            #endif
            for (int i = 1; i <= nx; ++i) {
                for (int j = 1; j <= ny; ++j) {
                    if ((i + j) % 2 != color) continue;

                    double rhs = computeRHS(grid, i, j, dt, rho);

                    double p_new = (
                        (grid.p[i + 1][j] + grid.p[i - 1][j]) / dx2 +
                        (grid.p[i][j + 1] + grid.p[i][j - 1]) / dy2 -
                        rhs
                    ) / factor;

                    double residual = std::abs(p_new - grid.p[i][j]);
                    maxResidual = std::max(maxResidual, residual);

                    grid.p[i][j] = grid.p[i][j] + omega * (p_new - grid.p[i][j]);
                }
            }
        }

        // 境界条件を適用
        bc.applyPressureBC(grid);

        lastResidual = maxResidual;

        if (maxResidual < tolerance) {
            return iter + 1;
        }
    }

    return -1;  // 収束しなかった
}

} // namespace fluid
