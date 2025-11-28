#include "CSVWriter.hpp"
#include "Constants.hpp"
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <filesystem>

namespace fluid {

CSVWriter::CSVWriter(const std::string& outputDir)
    : outputDir(outputDir), p_ref(constants::ATMOSPHERIC_PRESSURE)
{
}

void CSVWriter::createOutputDirectory() const {
    // std::filesystem を使用してポータブルなディレクトリ作成
    std::filesystem::create_directories(getDataDir());
    std::filesystem::create_directories(getFiguresDir());
}

std::string CSVWriter::getDataDir() const {
    return outputDir + "/data";
}

std::string CSVWriter::getFiguresDir() const {
    return outputDir + "/figures";
}

std::string CSVWriter::getFilename(const std::string& prefix, int step) const {
    std::ostringstream oss;
    oss << getDataDir() << "/" << prefix << "_"
        << std::setfill('0') << std::setw(6) << step << ".csv";
    return oss.str();
}

void CSVWriter::writeVelocity(const Grid& grid, int step, double time) const {
    const std::string filename = getFilename("velocity", step);
    std::ofstream ofs(filename);

    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    // ヘッダ
    ofs << "# time=" << time << "\n";
    ofs << "x,y,u,v,magnitude\n";

    // セル中心での速度を出力
    for (int i = 0; i < grid.nx; ++i) {
        for (int j = 0; j < grid.ny; ++j) {
            const double x = grid.cellCenterX(i);
            const double y = grid.cellCenterY(j);
            const double u = grid.uAtCellCenter(i, j);
            const double v = grid.vAtCellCenter(i, j);
            const double mag = grid.velocityMagnitude(i, j);

            ofs << std::fixed << std::setprecision(8)
                << x << "," << y << "," << u << "," << v << "," << mag << "\n";
        }
    }
}

void CSVWriter::writePressure(const Grid& grid, int step, double time) const {
    const std::string filename = getFilename("pressure", step);
    std::ofstream ofs(filename);

    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    // ヘッダ
    ofs << "# time=" << time << "\n";
    ofs << "x,y,p\n";

    // セル中心での圧力を出力（ゲージ圧 + 参照圧力 = 絶対圧力）
    for (int i = 0; i < grid.nx; ++i) {
        for (int j = 0; j < grid.ny; ++j) {
            const double x = grid.cellCenterX(i);
            const double y = grid.cellCenterY(j);
            const double p = grid.p[i + 1][j + 1] + p_ref;

            ofs << std::fixed << std::setprecision(8)
                << x << "," << y << "," << p << "\n";
        }
    }
}

void CSVWriter::writeAll(const Grid& grid, int step, double time) const {
    const std::string filename = getFilename("field", step);
    std::ofstream ofs(filename);

    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    // ヘッダ
    ofs << "# time=" << time << "\n";
    ofs << "x,y,u,v,p,magnitude\n";

    // 全フィールドを出力（圧力は絶対圧力で出力）
    for (int i = 0; i < grid.nx; ++i) {
        for (int j = 0; j < grid.ny; ++j) {
            const double x = grid.cellCenterX(i);
            const double y = grid.cellCenterY(j);
            const double u = grid.uAtCellCenter(i, j);
            const double v = grid.vAtCellCenter(i, j);
            const double p = grid.p[i + 1][j + 1] + p_ref;
            const double mag = grid.velocityMagnitude(i, j);

            ofs << std::fixed << std::setprecision(8)
                << x << "," << y << "," << u << "," << v << "," << p << "," << mag << "\n";
        }
    }
}

void CSVWriter::writeMetadata(const Grid& grid) const {
    const std::string filename = getDataDir() + "/metadata.csv";
    std::ofstream ofs(filename);

    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    ofs << "parameter,value\n";
    ofs << "nx," << grid.nx << "\n";
    ofs << "ny," << grid.ny << "\n";
    ofs << "lx," << grid.lx << "\n";
    ofs << "ly," << grid.ly << "\n";
    ofs << "dx," << grid.dx << "\n";
    ofs << "dy," << grid.dy << "\n";
}

void CSVWriter::writeSimulationLog(const std::string& solverType,
                                    const Grid& grid,
                                    double rho, double nu,
                                    double endTime, int totalSteps,
                                    double wallTime, double Re) const {
    const std::string filename = outputDir + "/simulation.log";
    std::ofstream ofs(filename);

    if (!ofs.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    // 現在時刻を取得
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);

    ofs << "========================================\n";
    ofs << "  Simulation Log\n";
    ofs << "========================================\n";
    ofs << "Date: " << std::ctime(&time_t_now);
    ofs << "\n";

    ofs << "[Solver]\n";
    ofs << "  Type: " << solverType << "\n";
    ofs << "  Scheme: 1st-order upwind (convection), 2nd-order central (diffusion)\n";
    ofs << "  Pressure solver: SOR (omega=1.8, tol=1e-6)\n";
    ofs << "\n";

    ofs << "[Grid]\n";
    ofs << "  Resolution: " << grid.nx << " x " << grid.ny << "\n";
    ofs << "  Domain size: " << grid.lx * 1000 << " mm x " << grid.ly * 1000 << " mm\n";
    ofs << "  Cell size: dx=" << grid.dx * 1000 << " mm, dy=" << grid.dy * 1000 << " mm\n";
    ofs << "\n";

    ofs << "[Physical Parameters]\n";
    ofs << "  Density (rho): " << rho << " kg/m^3\n";
    ofs << "  Kinematic viscosity (nu): " << nu << " m^2/s\n";
    ofs << "  Reynolds number: " << Re << "\n";
    ofs << "\n";

    ofs << "[Simulation]\n";
    ofs << "  End time: " << endTime << " s\n";
    ofs << "  Total steps: " << totalSteps << "\n";
    ofs << "\n";

    ofs << "[Performance]\n";
    ofs << "  Wall time: " << std::fixed << std::setprecision(2) << wallTime << " s\n";
    if (totalSteps > 0) {
        ofs << "  Time per step: " << std::fixed << std::setprecision(4)
            << wallTime / totalSteps * 1000 << " ms\n";
    }
    ofs << "========================================\n";
}

} // namespace fluid
