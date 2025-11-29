/**
 * チャネル流れ（Channel Flow / Hagen-Poiseuille Flow）
 *
 * 左から放物線プロファイルで流入、右から流出、上下は滑りなし壁
 * 2枚の平行平板間の完全発達層流（Hagen-Poiseuille流れ）
 *
 * 理論解との比較用：
 *   u(y) = U_max * (1 - (2y/H - 1)^2)
 *   U_mean = (2/3) * U_max
 *
 * 物理量（SI単位系）で計算
 */

#include "Grid.hpp"
#include "Solver.hpp"
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
    double Lx = 0.03;          // チャネル長さ [m] (30 mm) = 10L
    double rho = WATER_DENSITY;       // 密度 [kg/m³] (水)
    double mu = WATER_DYNAMIC_VISCOSITY;  // 粘性係数 [Pa·s] (水)
    double nu = mu / rho;      // 動粘性係数 [m²/s]
    double U_max = 0.015;      // 最大流速（中心速度） [m/s] (15 mm/s)
    double U_mean = POISEUILLE_MEAN_MAX_RATIO * U_max;  // 平均流速 [m/s] (Hagen-Poiseuille: U_mean = 2/3 * U_max)

    // レイノルズ数: Re = U_mean * L / ν（水力直径ではなくチャネル高さで定義）
    double Re = U_mean * L / nu;

    // 計算パラメータ
    int nx = 128;              // 格子数（流れ方向）
    int ny = 32;               // 格子数（壁面垂直方向）
    double endTime = 10.0;     // 終了時刻 [s]

    // コマンドライン引数の処理
    if (argc > 1) nx = std::atoi(argv[1]);
    if (argc > 2) ny = std::atoi(argv[2]);
    if (argc > 3) U_max = std::atof(argv[3]);
    if (argc > 4) endTime = std::atof(argv[4]);

    // 平均流速とレイノルズ数を再計算
    U_mean = POISEUILLE_MEAN_MAX_RATIO * U_max;
    Re = U_mean * L / nu;

    std::cout << "=== Channel Flow (Hagen-Poiseuille Flow) ===" << std::endl;
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
    std::cout << "=======================================" << std::endl;

    // 格子の作成 [m]
    fluid::Grid grid(nx, ny, Lx, L);

    // ソルバーの作成（物理パラメータを指定）
    fluid::Solver solver(rho, nu);
    solver.autoTimeStep = true;
    solver.cfl = DEFAULT_CFL;

    // 境界条件（放物線流入 = Hagen-Poiseuille、右流出、上下滑りなし壁）
    fluid::BoundaryCondition bc = fluid::BoundaryCondition::channelFlowParabolic(U_max);

    // CSV出力（Projection法）
    fluid::CSVWriter writer("output/channel_projection");
    writer.p_ref = bc.p_ref;  // 参照圧力を設定（絶対圧力で出力するため）
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
    writer.writeSimulationLog("Projection", grid, rho, nu, endTime,
                              solver.getStepCount(), wallTime, Re);

    return 0;
}
