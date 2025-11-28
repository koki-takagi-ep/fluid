/**
 * @file Solver.cpp
 * @brief Projection法ソルバーの実装
 */

#include "Solver.hpp"

namespace fluid {

Solver::Solver(double rho, double nu, double dt)
    : SolverBase(rho, nu, dt)
{
}

int Solver::step(Grid& grid, const BoundaryCondition& bc) {
    // 適応時間刻み
    if (autoTimeStep) {
        dt = computeTimeStep(grid);
    }

    // Step 1: 中間速度場の計算（基底クラスのメソッド）
    computeIntermediateVelocity(grid, bc);

    // Step 2: 圧力ポアソン方程式を解く
    int pressureIter = pressureSolver.solve(grid, bc, dt, rho);

    // Step 3: 速度の補正（基底クラスのメソッド）
    correctVelocity(grid);

    // 境界条件を適用
    bc.apply(grid);

    time_ += dt;
    stepCount_++;

    return pressureIter;
}

} // namespace fluid
