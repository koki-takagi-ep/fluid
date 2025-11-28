/**
 * チャネル流れ（Channel Flow / Poiseuille Flow）
 *
 * 左から流入、右から流出、上下は滑りなし壁
 * 2枚の平行平板間の流れ
 *
 * 物理量（SI単位系）で計算
 */

#include "Grid.hpp"
#include "Solver.hpp"
#include "BoundaryCondition.hpp"
#include "CSVWriter.hpp"
#include <iostream>
#include <cstdlib>

int main(int argc, char* argv[]) {
    // 物理パラメータ（SI単位系）
    double L = 0.003;          // チャネル高さ [m] (3 mm)
    double Lx = 0.03;          // チャネル長さ [m] (30 mm) = 10L
    double rho = 1000.0;       // 密度 [kg/m³] (水)
    double mu = 1.0e-3;        // 粘性係数 [Pa·s] (水)
    double nu = mu / rho;      // 動粘性係数 [m²/s]
    double U_in = 0.01;        // 流入速度 [m/s] (10 mm/s)

    // レイノルズ数: Re = ρUL/μ = UL/ν
    double Re = U_in * L / nu;

    // 計算パラメータ
    int nx = 128;              // 格子数（流れ方向）
    int ny = 32;               // 格子数（壁面垂直方向）
    double endTime = 10.0;     // 終了時刻 [s]

    // コマンドライン引数の処理
    if (argc > 1) nx = std::atoi(argv[1]);
    if (argc > 2) ny = std::atoi(argv[2]);
    if (argc > 3) U_in = std::atof(argv[3]);
    if (argc > 4) endTime = std::atof(argv[4]);

    // レイノルズ数を再計算
    Re = U_in * L / nu;

    std::cout << "=== Channel Flow (Poiseuille Flow) ===" << std::endl;
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
    std::cout << "=======================================" << std::endl;

    // 格子の作成 [m]
    fluid::Grid grid(nx, ny, Lx, L);

    // ソルバーの作成（物理パラメータを指定）
    fluid::Solver solver(rho, nu);
    solver.autoTimeStep = true;
    solver.cfl = 0.5;

    // 境界条件（チャネル流れ: 左流入、右流出、上下滑りなし壁）
    fluid::BoundaryCondition bc = fluid::BoundaryCondition::channelFlow(U_in);

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

    // シミュレーション実行
    solver.run(grid, bc, endTime, 100, outputCallback);

    std::cout << "\nSimulation completed!" << std::endl;
    std::cout << "Total steps: " << solver.getStepCount() << std::endl;
    std::cout << "Final time: " << solver.getTime() << " s" << std::endl;
    std::cout << "Output files: " << outputCount << std::endl;

    return 0;
}
