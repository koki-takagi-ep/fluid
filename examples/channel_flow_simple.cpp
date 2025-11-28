/**
 * チャネル流れ（SIMPLE法）
 *
 * SIMPLE法を使用したチャネル流れのシミュレーション
 * Projection法との比較用
 */

#include "Grid.hpp"
#include "SimpleSolver.hpp"
#include "BoundaryCondition.hpp"
#include "CSVWriter.hpp"
#include <iostream>
#include <cstdlib>

int main(int argc, char* argv[]) {
    // 物理パラメータ（SI単位系）
    double L = 0.003;          // チャネル高さ [m] (3 mm)
    double Lx = 0.03;          // チャネル長さ [m] (30 mm)
    double rho = 1000.0;       // 密度 [kg/m³] (水)
    double mu = 1.0e-3;        // 粘性係数 [Pa·s] (水)
    double nu = mu / rho;      // 動粘性係数 [m²/s]
    double U_in = 0.01;        // 流入速度 [m/s] (10 mm/s)

    // レイノルズ数
    double Re = U_in * L / nu;

    // 計算パラメータ
    int nx = 128;
    int ny = 32;
    double endTime = 10.0;

    // コマンドライン引数の処理
    if (argc > 1) nx = std::atoi(argv[1]);
    if (argc > 2) ny = std::atoi(argv[2]);
    if (argc > 3) U_in = std::atof(argv[3]);
    if (argc > 4) endTime = std::atof(argv[4]);

    Re = U_in * L / nu;

    std::cout << "=== Channel Flow (SIMPLE Method) ===" << std::endl;
    std::cout << "Physical parameters:" << std::endl;
    std::cout << "  Channel height: " << L * 1000 << " mm" << std::endl;
    std::cout << "  Channel length: " << Lx * 1000 << " mm" << std::endl;
    std::cout << "  Density: " << rho << " kg/m^3" << std::endl;
    std::cout << "  Viscosity: " << mu << " Pa.s" << std::endl;
    std::cout << "  Inflow velocity: " << U_in * 1000 << " mm/s" << std::endl;
    std::cout << "  Reynolds number: " << Re << std::endl;
    std::cout << "Numerical parameters:" << std::endl;
    std::cout << "  Grid: " << nx << " x " << ny << std::endl;
    std::cout << "  End time: " << endTime << " s" << std::endl;
    std::cout << "=====================================" << std::endl;

    // 格子の作成
    fluid::Grid grid(nx, ny, Lx, L);

    // SIMPLE法ソルバーの作成
    fluid::SimpleSolver solver(rho, nu);
    solver.autoTimeStep = true;
    solver.cfl = 0.5;
    solver.alpha_u = 0.7;  // 速度緩和係数
    solver.alpha_p = 0.3;  // 圧力緩和係数

    // 境界条件
    fluid::BoundaryCondition bc = fluid::BoundaryCondition::channelFlow(U_in);

    // CSV出力（SIMPLE法）
    fluid::CSVWriter writer("output/channel_simple");
    writer.p_ref = bc.p_ref;
    writer.createOutputDirectory();
    writer.writeMetadata(grid);

    // 出力コールバック
    int outputCount = 0;
    auto outputCallback = [&](int /*step*/, double time, const fluid::Grid& g) {
        writer.writeAll(g, outputCount, time);
        outputCount++;
    };

    // シミュレーション実行
    solver.run(grid, bc, endTime, 100, outputCallback);

    std::cout << "\nSimulation completed!" << std::endl;
    std::cout << "Total steps: " << solver.getStepCount() << std::endl;
    std::cout << "Final time: " << solver.getTime() << " s" << std::endl;
    std::cout << "Output files: " << outputCount << std::endl;

    return 0;
}
