#pragma once

#include "Grid.hpp"
#include <string>
#include <fstream>
#include <chrono>

namespace fluid {

/**
 * CSV出力クラス
 *
 * 出力ディレクトリ構造:
 *   outputDir/
 *   ├── data/           # CSVファイル
 *   │   ├── metadata.csv
 *   │   └── field_XXXXXX.csv
 *   ├── figures/        # 可視化結果（Pythonスクリプトが出力）
 *   └── simulation.log  # 計算ログ
 */
class CSVWriter {
public:
    std::string outputDir;
    double p_ref = 101325.0;  ///< 参照圧力（大気圧） [Pa]

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

    /**
     * シミュレーションログを出力
     * @param solverType ソルバーの種類 ("Projection", "SIMPLE")
     * @param grid 格子情報
     * @param rho 密度 [kg/m³]
     * @param nu 動粘性係数 [m²/s]
     * @param endTime シミュレーション終了時刻 [s]
     * @param totalSteps 総ステップ数
     * @param wallTime 実行時間 [s]
     * @param Re レイノルズ数
     */
    void writeSimulationLog(const std::string& solverType,
                            const Grid& grid,
                            double rho, double nu,
                            double endTime, int totalSteps,
                            double wallTime, double Re) const;

private:
    std::string getFilename(const std::string& prefix, int step) const;
};

} // namespace fluid
