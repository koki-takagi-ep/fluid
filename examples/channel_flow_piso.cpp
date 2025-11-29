/**
 * チャネル流れ（PISO法 / Hagen-Poiseuille Flow）
 *
 * PISO法を使用したチャネル流れのシミュレーション
 * 放物線流入プロファイルでHagen-Poiseuille理論解と比較
 */

#include "Grid.hpp"
#include "PisoSolver.hpp"
#include "BoundaryCondition.hpp"
#include "CSVWriter.hpp"
#include "Constants.hpp"
#include <iostream>
#include <cstdlib>
#include <chrono>

using namespace fluid::constants;

int main(int argc, char* argv[]) {
    // 物理パラメータ（SI単位系）
    double L = 0.003;          // チャネル高さ [m] (3 mm)
    double Lx = 0.03;          // チャネル長さ [m] (30 mm)
    double rho = WATER_DENSITY;       // 密度 [kg/m³] (水)
    double mu = WATER_DYNAMIC_VISCOSITY;  // 粘性係数 [Pa·s] (水)
    double nu = mu / rho;      // 動粘性係数 [m²/s]
    double U_max = 0.015;      // 最大流速（中心速度） [m/s] (15 mm/s)
    double U_mean = POISEUILLE_MEAN_MAX_RATIO * U_max;  // 平均流速 [m/s]

    // レイノルズ数
    double Re = U_mean * L / nu;

    // 計算パラメータ
    int nx = 128;
    int ny = 32;
    double endTime = 10.0;
    int nCorrectors = 2;  // PISO補正ステップ数

    // コマンドライン引数の処理
    if (argc > 1) nx = std::atoi(argv[1]);
    if (argc > 2) ny = std::atoi(argv[2]);
    if (argc > 3) U_max = std::atof(argv[3]);
    if (argc > 4) endTime = std::atof(argv[4]);
    if (argc > 5) nCorrectors = std::atoi(argv[5]);

    U_mean = POISEUILLE_MEAN_MAX_RATIO * U_max;
    Re = U_mean * L / nu;

    std::cout << "=== Channel Flow - PISO (Hagen-Poiseuille) ===" << std::endl;
    std::cout << "Physical parameters:" << std::endl;
    std::cout << "  Channel height: " << L * 1000 << " mm" << std::endl;
    std::cout << "  Channel length: " << Lx * 1000 << " mm" << std::endl;
    std::cout << "  Density: " << rho << " kg/m^3" << std::endl;
    std::cout << "  Viscosity: " << mu << " Pa.s" << std::endl;
    std::cout << "  Max velocity (center): " << U_max * 1000 << " mm/s" << std::endl;
    std::cout << "  Mean velocity: " << U_mean * 1000 << " mm/s" << std::endl;
    std::cout << "  Reynolds number: " << Re << std::endl;
    std::cout << "  Inflow profile: Parabolic (Hagen-Poiseuille)" << std::endl;
    std::cout << "Numerical parameters:" << std::endl;
    std::cout << "  Grid: " << nx << " x " << ny << std::endl;
    std::cout << "  End time: " << endTime << " s" << std::endl;
    std::cout << "  PISO correctors: " << nCorrectors << std::endl;
    std::cout << "===============================================" << std::endl;

    // 格子の作成
    fluid::Grid grid(nx, ny, Lx, L);

    // PISO法ソルバーの作成
    fluid::PisoSolver solver(rho, nu, 0.001, nCorrectors);
    solver.autoTimeStep = true;
    solver.cfl = DEFAULT_CFL;

    // 境界条件（放物線流入 = Hagen-Poiseuille）
    fluid::BoundaryCondition bc = fluid::BoundaryCondition::channelFlowParabolic(U_max);

    // CSV出力（PISO法）
    fluid::CSVWriter writer("output/channel_piso");
    writer.p_ref = bc.p_ref;
    writer.createOutputDirectory();
    writer.writeMetadata(grid);

    // 出力コールバック
    int outputCount = 0;
    auto outputCallback = [&](int /*step*/, double time, const fluid::Grid& g) {
        writer.writeAll(g, outputCount, time);
        outputCount++;
    };

    // 計算時間の計測開始
    auto startTime = std::chrono::high_resolution_clock::now();

    // シミュレーション実行
    solver.run(grid, bc, endTime, 100, outputCallback);

    // 計算時間の計測終了
    auto endTimePoint = std::chrono::high_resolution_clock::now();
    double wallTime = std::chrono::duration<double>(endTimePoint - startTime).count();

    std::cout << "\nSimulation completed!" << std::endl;
    std::cout << "Total steps: " << solver.getStepCount() << std::endl;
    std::cout << "Final time: " << solver.getTime() << " s" << std::endl;
    std::cout << "Output files: " << outputCount << std::endl;
    std::cout << "Wall time: " << wallTime << " s" << std::endl;

    // シミュレーションログを出力
    writer.writeSimulationLog("PISO", grid, rho, nu, endTime,
                              solver.getStepCount(), wallTime, Re);

    return 0;
}
