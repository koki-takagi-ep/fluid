/**
 * キャビティ流れ（SIMPLE法）
 *
 * SIMPLE法を使用したキャビティ流れのシミュレーション
 * Projection法との比較用
 */

#include "Grid.hpp"
#include "SimpleSolver.hpp"
#include "BoundaryCondition.hpp"
#include "CSVWriter.hpp"
#include "Constants.hpp"
#include <iostream>
#include <cstdlib>
#include <chrono>

using namespace fluid::constants;

int main(int argc, char* argv[]) {
    // 物理パラメータ（SI単位系）
    double L = 0.01;           // キャビティサイズ [m] (10 mm)
    double rho = WATER_DENSITY;       // 密度 [kg/m³] (水)
    double mu = WATER_DYNAMIC_VISCOSITY;  // 粘性係数 [Pa·s] (水)
    double nu = mu / rho;      // 動粘性係数 [m²/s]
    double U_lid = 0.01;       // 上壁速度 [m/s] (10 mm/s)

    // レイノルズ数: Re = ρUL/μ = UL/ν
    double Re = U_lid * L / nu;

    // 計算パラメータ
    int nx = 64;               // 格子数
    int ny = 64;
    double endTime = DEFAULT_CAVITY_END_TIME;  // 終了時刻 [s]

    // コマンドライン引数の処理
    if (argc > 1) nx = ny = std::atoi(argv[1]);
    if (argc > 2) U_lid = std::atof(argv[2]);
    if (argc > 3) endTime = std::atof(argv[3]);

    // レイノルズ数を再計算
    Re = U_lid * L / nu;

    std::cout << "=== Lid-Driven Cavity Flow (SIMPLE Method) ===" << std::endl;
    std::cout << "Physical parameters:" << std::endl;
    std::cout << "  Cavity size: " << L * 1000 << " mm" << std::endl;
    std::cout << "  Density: " << rho << " kg/m^3" << std::endl;
    std::cout << "  Viscosity: " << mu << " Pa.s" << std::endl;
    std::cout << "  Lid velocity: " << U_lid * 1000 << " mm/s" << std::endl;
    std::cout << "  Reynolds number: " << Re << std::endl;
    std::cout << "Numerical parameters:" << std::endl;
    std::cout << "  Grid: " << nx << " x " << ny << std::endl;
    std::cout << "  End time: " << endTime << " s" << std::endl;
    std::cout << "===============================================" << std::endl;

    // 格子の作成 [m]
    fluid::Grid grid(nx, ny, L, L);

    // SIMPLE法ソルバーの作成（物理パラメータを指定）
    fluid::SimpleSolver solver(rho, nu);
    solver.autoTimeStep = true;
    solver.cfl = DEFAULT_CFL;
    solver.alpha_u = DEFAULT_VELOCITY_RELAXATION;  // 速度緩和係数
    solver.alpha_p = DEFAULT_PRESSURE_RELAXATION;  // 圧力緩和係数

    // 境界条件（キャビティ流れ）
    fluid::BoundaryCondition bc = fluid::BoundaryCondition::cavityFlow(U_lid);

    // CSV出力（SIMPLE法）
    fluid::CSVWriter writer("output/cavity_simple");
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
    writer.writeSimulationLog("SIMPLE", grid, rho, nu, endTime,
                              solver.getStepCount(), wallTime, Re);

    return 0;
}
