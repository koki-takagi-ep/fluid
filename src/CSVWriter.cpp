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

} // namespace fluid
