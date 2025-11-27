#include "CSVWriter.hpp"
#include <iomanip>
#include <sstream>
#include <sys/stat.h>

namespace fluid {

CSVWriter::CSVWriter(const std::string& outputDir)
    : outputDir(outputDir)
{
}

void CSVWriter::createOutputDirectory() const {
    mkdir(outputDir.c_str(), 0755);
    mkdir(getDataDir().c_str(), 0755);
    mkdir(getFiguresDir().c_str(), 0755);
}

std::string CSVWriter::getDataDir() const {
    return outputDir + "/data";
}

std::string CSVWriter::getFiguresDir() const {
    return outputDir + "/figures";
}

std::string CSVWriter::getFilename(const std::string& prefix, int step) const {
    std::ostringstream oss;
    oss << getDataDir() << "/" << prefix << "_" << std::setfill('0') << std::setw(6) << step << ".csv";
    return oss.str();
}

void CSVWriter::writeVelocity(const Grid& grid, int step, double time) const {
    std::string filename = getFilename("velocity", step);
    std::ofstream ofs(filename);

    if (!ofs.is_open()) {
        return;
    }

    // ヘッダ
    ofs << "# time=" << time << "\n";
    ofs << "x,y,u,v,magnitude\n";

    // セル中心での速度を出力
    for (int i = 0; i < grid.nx; ++i) {
        for (int j = 0; j < grid.ny; ++j) {
            double x = grid.cellCenterX(i);
            double y = grid.cellCenterY(j);
            double u = grid.uAtCellCenter(i, j);
            double v = grid.vAtCellCenter(i, j);
            double mag = grid.velocityMagnitude(i, j);

            ofs << std::fixed << std::setprecision(8)
                << x << "," << y << "," << u << "," << v << "," << mag << "\n";
        }
    }

    ofs.close();
}

void CSVWriter::writePressure(const Grid& grid, int step, double time) const {
    std::string filename = getFilename("pressure", step);
    std::ofstream ofs(filename);

    if (!ofs.is_open()) {
        return;
    }

    // ヘッダ
    ofs << "# time=" << time << "\n";
    ofs << "x,y,p\n";

    // セル中心での圧力を出力
    for (int i = 0; i < grid.nx; ++i) {
        for (int j = 0; j < grid.ny; ++j) {
            double x = grid.cellCenterX(i);
            double y = grid.cellCenterY(j);
            double p = grid.p[i + 1][j + 1];

            ofs << std::fixed << std::setprecision(8)
                << x << "," << y << "," << p << "\n";
        }
    }

    ofs.close();
}

void CSVWriter::writeAll(const Grid& grid, int step, double time) const {
    std::string filename = getFilename("field", step);
    std::ofstream ofs(filename);

    if (!ofs.is_open()) {
        return;
    }

    // ヘッダ
    ofs << "# time=" << time << "\n";
    ofs << "x,y,u,v,p,magnitude\n";

    // 全フィールドを出力
    for (int i = 0; i < grid.nx; ++i) {
        for (int j = 0; j < grid.ny; ++j) {
            double x = grid.cellCenterX(i);
            double y = grid.cellCenterY(j);
            double u = grid.uAtCellCenter(i, j);
            double v = grid.vAtCellCenter(i, j);
            double p = grid.p[i + 1][j + 1];
            double mag = grid.velocityMagnitude(i, j);

            ofs << std::fixed << std::setprecision(8)
                << x << "," << y << "," << u << "," << v << "," << p << "," << mag << "\n";
        }
    }

    ofs.close();
}

void CSVWriter::writeMetadata(const Grid& grid) const {
    std::string filename = getDataDir() + "/metadata.csv";
    std::ofstream ofs(filename);

    if (!ofs.is_open()) {
        return;
    }

    ofs << "parameter,value\n";
    ofs << "nx," << grid.nx << "\n";
    ofs << "ny," << grid.ny << "\n";
    ofs << "lx," << grid.lx << "\n";
    ofs << "ly," << grid.ly << "\n";
    ofs << "dx," << grid.dx << "\n";
    ofs << "dy," << grid.dy << "\n";

    ofs.close();
}

} // namespace fluid
