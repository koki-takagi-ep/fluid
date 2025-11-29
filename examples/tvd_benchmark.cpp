/**
 * TVDスキームのベンチマーク
 *
 * 各リミッター（None, Minmod, Superbee, Van Leer, MC）を使用して
 * キャビティ流れを計算し、精度と計算時間を比較する。
 */

#include "Grid.hpp"
#include "Solver.hpp"
#include "BoundaryCondition.hpp"
#include "CSVWriter.hpp"
#include "FluxLimiter.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>

struct BenchmarkResult {
    std::string name;
    double wallTime;
    int steps;
    double maxDivergence;
    double centerlineU;  // y=0.5でのu速度
};

BenchmarkResult runBenchmark(fluid::LimiterType limiterType, const std::string& name,
                              int nx, double U_lid, double endTime) {
    // 物理パラメータ
    double L = 0.01;           // キャビティサイズ [m]
    double rho = 1000.0;       // 密度 [kg/m³]
    double mu = 1.0e-3;        // 粘性係数 [Pa·s]
    double nu = mu / rho;      // 動粘性係数 [m²/s]

    // 格子の作成
    fluid::Grid grid(nx, nx, L, L);

    // ソルバーの作成
    fluid::Solver solver(rho, nu);
    solver.autoTimeStep = true;
    solver.cfl = 0.5;
    solver.setLimiter(limiterType);

    // 境界条件
    fluid::BoundaryCondition bc = fluid::BoundaryCondition::cavityFlow(U_lid);

    // CSV出力
    std::string outputDir = "output/tvd_" + name;
    fluid::CSVWriter writer(outputDir);
    writer.createOutputDirectory();
    writer.writeMetadata(grid);

    int outputCount = 0;
    auto outputCallback = [&](int /*step*/, double time, const fluid::Grid& g) {
        writer.writeAll(g, outputCount, time);
        outputCount++;
    };

    // 計算時間の計測
    auto startTime = std::chrono::high_resolution_clock::now();

    // サイレント実行（出力を抑制）
    bc.apply(grid);
    outputCallback(0, 0.0, grid);

    double time = 0.0;
    int stepCount = 0;
    while (time < endTime) {
        solver.step(grid, bc);
        time = solver.getTime();
        stepCount = solver.getStepCount();

        if (stepCount % 500 == 0) {
            outputCallback(stepCount, time, grid);
        }
    }
    outputCallback(stepCount, time, grid);

    auto endTimePoint = std::chrono::high_resolution_clock::now();
    double wallTime = std::chrono::duration<double>(endTimePoint - startTime).count();

    // 結果の収集
    BenchmarkResult result;
    result.name = name;
    result.wallTime = wallTime;
    result.steps = stepCount;
    result.maxDivergence = grid.maxDivergence();

    // 中心線のu速度を取得
    int i_center = nx / 2;
    int j_center = nx / 2;
    result.centerlineU = grid.u[i_center][j_center] / U_lid;

    return result;
}

int main(int argc, char* argv[]) {
    // パラメータ
    int nx = 64;
    double U_lid = 0.01;
    double endTime = 10.0;

    if (argc > 1) nx = std::atoi(argv[1]);
    if (argc > 2) U_lid = std::atof(argv[2]);
    if (argc > 3) endTime = std::atof(argv[3]);

    std::cout << "=== TVD Scheme Benchmark ===" << std::endl;
    std::cout << "Grid: " << nx << " x " << nx << std::endl;
    std::cout << "Lid velocity: " << U_lid * 1000 << " mm/s" << std::endl;
    std::cout << "End time: " << endTime << " s" << std::endl;
    std::cout << "============================\n" << std::endl;

    // テストするリミッター
    std::vector<std::pair<fluid::LimiterType, std::string>> limiters = {
        {fluid::LimiterType::None, "upwind"},
        {fluid::LimiterType::Minmod, "minmod"},
        {fluid::LimiterType::Superbee, "superbee"},
        {fluid::LimiterType::VanLeer, "vanleer"},
        {fluid::LimiterType::MC, "mc"}
    };

    std::vector<BenchmarkResult> results;

    for (const auto& [limiterType, name] : limiters) {
        std::cout << "Running: " << fluid::limiterTypeToString(limiterType) << "... " << std::flush;
        BenchmarkResult result = runBenchmark(limiterType, name, nx, U_lid, endTime);
        results.push_back(result);
        std::cout << "done (" << std::fixed << std::setprecision(2)
                  << result.wallTime << " s)" << std::endl;
    }

    // 結果の表示
    std::cout << "\n=== Benchmark Results ===" << std::endl;
    std::cout << std::left << std::setw(25) << "Limiter"
              << std::right << std::setw(12) << "Time [s]"
              << std::setw(10) << "Steps"
              << std::setw(15) << "Max Div"
              << std::setw(12) << "u(0.5,0.5)" << std::endl;
    std::cout << std::string(74, '-') << std::endl;

    for (const auto& r : results) {
        std::cout << std::left << std::setw(25) << r.name
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(12) << r.wallTime
                  << std::setw(10) << r.steps
                  << std::scientific << std::setprecision(2)
                  << std::setw(15) << r.maxDivergence
                  << std::fixed << std::setprecision(4)
                  << std::setw(12) << r.centerlineU << std::endl;
    }

    std::cout << "\nOutput directories created: output/tvd_*" << std::endl;

    return 0;
}
