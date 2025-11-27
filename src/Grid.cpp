#include "Grid.hpp"
#include <algorithm>

namespace fluid {

Grid::Grid(int nx, int ny, double lx, double ly)
    : nx(nx), ny(ny), lx(lx), ly(ly)
{
    dx = lx / nx;
    dy = ly / ny;

    // 配列の初期化
    p = make2DArray(nx + 2, ny + 2);
    u = make2DArray(nx + 1, ny + 2);
    v = make2DArray(nx + 2, ny + 1);
    u_star = make2DArray(nx + 1, ny + 2);
    v_star = make2DArray(nx + 2, ny + 1);
}

std::vector<std::vector<double>> Grid::make2DArray(int rows, int cols) {
    return std::vector<std::vector<double>>(rows, std::vector<double>(cols, 0.0));
}

double Grid::computeDivergence(int i, int j) const {
    // i, j は内部セルのインデックス (0 <= i < nx, 0 <= j < ny)
    // 配列インデックスは +1 のオフセットが必要
    return (u[i + 1][j + 1] - u[i][j + 1]) / dx +
           (v[i + 1][j + 1] - v[i + 1][j]) / dy;
}

double Grid::maxDivergence() const {
    double maxDiv = 0.0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            maxDiv = std::max(maxDiv, std::abs(computeDivergence(i, j)));
        }
    }
    return maxDiv;
}

double Grid::uAtCellCenter(int i, int j) const {
    // セル(i,j)の中心でのu速度を補間
    return 0.5 * (u[i][j + 1] + u[i + 1][j + 1]);
}

double Grid::vAtCellCenter(int i, int j) const {
    // セル(i,j)の中心でのv速度を補間
    return 0.5 * (v[i + 1][j] + v[i + 1][j + 1]);
}

double Grid::velocityMagnitude(int i, int j) const {
    double uc = uAtCellCenter(i, j);
    double vc = vAtCellCenter(i, j);
    return std::sqrt(uc * uc + vc * vc);
}

} // namespace fluid
