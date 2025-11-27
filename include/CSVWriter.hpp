#pragma once

#include "Grid.hpp"
#include <string>
#include <fstream>

namespace fluid {

/**
 * CSV出力クラス
 *
 * 出力ディレクトリ構造:
 *   outputDir/
 *   ├── data/           # CSVファイル
 *   │   ├── metadata.csv
 *   │   └── field_XXXXXX.csv
 *   └── figures/        # 可視化結果（Pythonスクリプトが出力）
 */
class CSVWriter {
public:
    std::string outputDir;

    CSVWriter(const std::string& outputDir = "output");

    /**
     * 出力ディレクトリを作成（data/, figures/ サブディレクトリも作成）
     */
    void createOutputDirectory() const;

    /**
     * dataディレクトリのパスを取得
     */
    std::string getDataDir() const;

    /**
     * figuresディレクトリのパスを取得
     */
    std::string getFiguresDir() const;

    /**
     * 速度場をCSV出力（セル中心値）
     * @param grid 格子
     * @param step ステップ番号
     * @param time 時刻
     */
    void writeVelocity(const Grid& grid, int step, double time) const;

    /**
     * 圧力場をCSV出力
     */
    void writePressure(const Grid& grid, int step, double time) const;

    /**
     * 全場をまとめてCSV出力
     */
    void writeAll(const Grid& grid, int step, double time) const;

    /**
     * メタデータ（格子情報）を出力
     */
    void writeMetadata(const Grid& grid) const;

private:
    std::string getFilename(const std::string& prefix, int step) const;
};

} // namespace fluid
