#include "BoundaryCondition.hpp"

namespace fluid {

void BoundaryCondition::apply(Grid& grid) const {
    applyLeftBC(grid);
    applyRightBC(grid);
    applyBottomBC(grid);
    applyTopBC(grid);
}

void BoundaryCondition::applyLeftBC(Grid& grid) const {
    int ny = grid.ny;

    switch (left) {
        case BCType::NoSlip:
            // u = 0 on left wall (i=0)
            for (int j = 0; j < ny + 2; ++j) {
                grid.u[0][j] = u_left;
            }
            // v: 反射条件
            for (int j = 0; j < ny + 1; ++j) {
                grid.v[0][j] = 2.0 * v_left - grid.v[1][j];
            }
            break;

        case BCType::Slip:
            for (int j = 0; j < ny + 2; ++j) {
                grid.u[0][j] = 0.0;
            }
            for (int j = 0; j < ny + 1; ++j) {
                grid.v[0][j] = grid.v[1][j];  // ∂v/∂x = 0
            }
            break;

        case BCType::Inflow:
            for (int j = 0; j < ny + 2; ++j) {
                grid.u[0][j] = u_left;
            }
            for (int j = 0; j < ny + 1; ++j) {
                grid.v[0][j] = 2.0 * v_left - grid.v[1][j];
            }
            break;

        case BCType::Outflow:
            for (int j = 0; j < ny + 2; ++j) {
                grid.u[0][j] = grid.u[1][j];  // ∂u/∂x = 0
            }
            for (int j = 0; j < ny + 1; ++j) {
                grid.v[0][j] = grid.v[1][j];
            }
            break;

        default:
            break;
    }
}

void BoundaryCondition::applyRightBC(Grid& grid) const {
    int nx = grid.nx;
    int ny = grid.ny;

    switch (right) {
        case BCType::NoSlip:
            for (int j = 0; j < ny + 2; ++j) {
                grid.u[nx][j] = u_right;
            }
            for (int j = 0; j < ny + 1; ++j) {
                grid.v[nx + 1][j] = 2.0 * v_right - grid.v[nx][j];
            }
            break;

        case BCType::Slip:
            for (int j = 0; j < ny + 2; ++j) {
                grid.u[nx][j] = 0.0;
            }
            for (int j = 0; j < ny + 1; ++j) {
                grid.v[nx + 1][j] = grid.v[nx][j];
            }
            break;

        case BCType::Outflow:
            for (int j = 0; j < ny + 2; ++j) {
                grid.u[nx][j] = grid.u[nx - 1][j];
            }
            for (int j = 0; j < ny + 1; ++j) {
                grid.v[nx + 1][j] = grid.v[nx][j];
            }
            break;

        default:
            break;
    }
}

void BoundaryCondition::applyBottomBC(Grid& grid) const {
    int nx = grid.nx;

    switch (bottom) {
        case BCType::NoSlip:
            for (int i = 0; i < nx + 2; ++i) {
                grid.v[i][0] = v_bottom;
            }
            for (int i = 0; i < nx + 1; ++i) {
                grid.u[i][0] = 2.0 * u_bottom - grid.u[i][1];
            }
            break;

        case BCType::Slip:
            for (int i = 0; i < nx + 2; ++i) {
                grid.v[i][0] = 0.0;
            }
            for (int i = 0; i < nx + 1; ++i) {
                grid.u[i][0] = grid.u[i][1];
            }
            break;

        default:
            break;
    }
}

void BoundaryCondition::applyTopBC(Grid& grid) const {
    int nx = grid.nx;
    int ny = grid.ny;

    switch (top) {
        case BCType::NoSlip:
            for (int i = 0; i < nx + 2; ++i) {
                grid.v[i][ny] = v_top;
            }
            for (int i = 0; i < nx + 1; ++i) {
                grid.u[i][ny + 1] = 2.0 * u_top - grid.u[i][ny];
            }
            break;

        case BCType::Slip:
            for (int i = 0; i < nx + 2; ++i) {
                grid.v[i][ny] = 0.0;
            }
            for (int i = 0; i < nx + 1; ++i) {
                grid.u[i][ny + 1] = grid.u[i][ny];
            }
            break;

        default:
            break;
    }
}

void BoundaryCondition::applyPressureBC(Grid& grid) const {
    int nx = grid.nx;
    int ny = grid.ny;

    // 注: 圧力はゲージ圧（大気圧基準の相対圧力）として計算する
    // 絶対圧力 = ゲージ圧 + p_ref（大気圧）
    // 流出境界でゲージ圧 = 0（絶対圧力 = 大気圧）となるように設定

    // 左境界: ノイマン条件（∂p/∂n = 0）
    for (int j = 1; j <= ny; ++j) {
        grid.p[0][j] = grid.p[1][j];
    }

    // 右境界
    if (right == BCType::Outflow) {
        // 流出境界: ディリクレ条件（ゲージ圧 = 0、つまり絶対圧力 = 大気圧）
        for (int j = 1; j <= ny; ++j) {
            grid.p[nx + 1][j] = 0.0;
        }
    } else {
        // ノイマン条件（∂p/∂n = 0）
        for (int j = 1; j <= ny; ++j) {
            grid.p[nx + 1][j] = grid.p[nx][j];
        }
    }

    // 下境界: ノイマン条件
    for (int i = 1; i <= nx; ++i) {
        grid.p[i][0] = grid.p[i][1];
    }

    // 上境界: ノイマン条件
    for (int i = 1; i <= nx; ++i) {
        grid.p[i][ny + 1] = grid.p[i][ny];
    }

    // 角
    grid.p[0][0] = 0.5 * (grid.p[1][0] + grid.p[0][1]);
    grid.p[nx + 1][0] = 0.5 * (grid.p[nx][0] + grid.p[nx + 1][1]);
    grid.p[0][ny + 1] = 0.5 * (grid.p[1][ny + 1] + grid.p[0][ny]);
    grid.p[nx + 1][ny + 1] = 0.5 * (grid.p[nx][ny + 1] + grid.p[nx + 1][ny]);
}

void BoundaryCondition::applyPressureBC(const Grid& grid, std::vector<std::vector<double>>& p_field) const {
    int nx = grid.nx;
    int ny = grid.ny;

    // 圧力補正場p'に対する境界条件
    // 全境界でNeumann条件（∂p'/∂n = 0）を適用

    // 左境界
    for (int j = 1; j <= ny; ++j) {
        p_field[0][j] = p_field[1][j];
    }

    // 右境界
    for (int j = 1; j <= ny; ++j) {
        p_field[nx + 1][j] = p_field[nx][j];
    }

    // 下境界
    for (int i = 1; i <= nx; ++i) {
        p_field[i][0] = p_field[i][1];
    }

    // 上境界
    for (int i = 1; i <= nx; ++i) {
        p_field[i][ny + 1] = p_field[i][ny];
    }

    // 角
    p_field[0][0] = 0.5 * (p_field[1][0] + p_field[0][1]);
    p_field[nx + 1][0] = 0.5 * (p_field[nx][0] + p_field[nx + 1][1]);
    p_field[0][ny + 1] = 0.5 * (p_field[1][ny + 1] + p_field[0][ny]);
    p_field[nx + 1][ny + 1] = 0.5 * (p_field[nx][ny + 1] + p_field[nx + 1][ny]);
}

BoundaryCondition BoundaryCondition::cavityFlow(double lidVelocity) {
    BoundaryCondition bc;
    bc.left = BCType::NoSlip;
    bc.right = BCType::NoSlip;
    bc.bottom = BCType::NoSlip;
    bc.top = BCType::NoSlip;

    // 上壁のみ動く
    bc.u_top = lidVelocity;

    return bc;
}

BoundaryCondition BoundaryCondition::channelFlow(double inflowVelocity) {
    BoundaryCondition bc;
    bc.left = BCType::Inflow;
    bc.right = BCType::Outflow;
    bc.bottom = BCType::NoSlip;
    bc.top = BCType::NoSlip;

    bc.u_left = inflowVelocity;

    return bc;
}

} // namespace fluid
