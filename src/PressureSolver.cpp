#include "PressureSolver.hpp"
#include "Constants.hpp"
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
    double factor = constants::LAPLACIAN_CENTER_COEFF * (1.0 / dx2 + 1.0 / dy2);

    // 圧力を0にリセットせず、前回の解を初期値として使用（収束を早める）
    // Grid側で初期化されているので、ここではリセットしない

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

                    // Poisson方程式の残差: ∇²p - rhs
                    double laplacian_p = (grid.p[i + 1][j] - 2.0 * grid.p[i][j] + grid.p[i - 1][j]) / dx2 +
                                         (grid.p[i][j + 1] - 2.0 * grid.p[i][j] + grid.p[i][j - 1]) / dy2;
                    double residual = std::abs(laplacian_p - rhs);
                    maxResidual = std::max(maxResidual, residual);

                    double p_new = (
                        (grid.p[i + 1][j] + grid.p[i - 1][j]) / dx2 +
                        (grid.p[i][j + 1] + grid.p[i][j - 1]) / dy2 -
                        rhs
                    ) / factor;

                    grid.p[i][j] = grid.p[i][j] + omega * (p_new - grid.p[i][j]);
                }
            }
        }

        // 境界条件を適用
        bc.applyPressureBC(grid);

        lastResidual = maxResidual;

        // 残差はρ/dtでスケールされているので、発散の閾値に変換
        // tolerance は発散の閾値として解釈（例: 1e-6 なら ∇・u < 1e-6）
        if (maxResidual < tolerance * rho / dt) {
            return iter + 1;
        }
    }

    return -1;  // 収束しなかった
}

} // namespace fluid
