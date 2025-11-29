/**
 * Channel flow convergence study with multiple solvers and TVD limiters
 *
 * Command-line arguments:
 *   ./channel_flow_convergence <solver> <limiter> <nx> <ny> <U_max> <end_time> <output_dir>
 *
 * Solver options: projection, simple, piso
 * Limiter options: none, minmod, superbee, vanleer, mc
 */

#include "Grid.hpp"
#include "Solver.hpp"
#include "SimpleSolver.hpp"
#include "PisoSolver.hpp"
#include "BoundaryCondition.hpp"
#include "CSVWriter.hpp"
#include "FluxLimiter.hpp"
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <memory>
#include <string>

int main(int argc, char* argv[]) {
    // Default parameters
    std::string solverType = "projection";
    std::string limiterStr = "none";
    int nx = 128;
    int ny = 32;
    double U_max = 0.015;
    double endTime = 5.0;
    std::string outputDir = "output/channel_convergence";

    // Parse command-line arguments
    if (argc > 1) solverType = argv[1];
    if (argc > 2) limiterStr = argv[2];
    if (argc > 3) nx = std::atoi(argv[3]);
    if (argc > 4) ny = std::atoi(argv[4]);
    if (argc > 5) U_max = std::atof(argv[5]);
    if (argc > 6) endTime = std::atof(argv[6]);
    if (argc > 7) outputDir = argv[7];

    // Physical parameters (SI units)
    double L = 0.003;          // Channel height [m] (3 mm)
    double Lx = 0.03;          // Channel length [m] (30 mm)
    double rho = 1000.0;       // Density [kg/m³]
    double mu = 1.0e-3;        // Dynamic viscosity [Pa·s]
    double nu = mu / rho;      // Kinematic viscosity [m²/s]
    double U_mean = (2.0/3.0) * U_max;
    double Re = U_mean * L / nu;

    // Parse limiter type
    fluid::LimiterType limiterType = fluid::stringToLimiterType(limiterStr);

    std::cout << "=== Channel Flow Convergence Study ===" << std::endl;
    std::cout << "Solver: " << solverType << std::endl;
    std::cout << "Limiter: " << fluid::limiterTypeToString(limiterType) << std::endl;
    std::cout << "Grid: " << nx << " x " << ny << std::endl;
    std::cout << "Re: " << Re << std::endl;
    std::cout << "Output: " << outputDir << std::endl;
    std::cout << "=======================================" << std::endl;

    // Create grid
    fluid::Grid grid(nx, ny, Lx, L);

    // Create solver based on type
    std::unique_ptr<fluid::SolverBase> solver;
    if (solverType == "simple") {
        solver = std::make_unique<fluid::SimpleSolver>(rho, nu);
    } else if (solverType == "piso") {
        solver = std::make_unique<fluid::PisoSolver>(rho, nu);
    } else {
        solver = std::make_unique<fluid::Solver>(rho, nu);
    }

    solver->autoTimeStep = true;
    solver->cfl = 0.5;
    solver->setLimiter(limiterType);

    // Boundary conditions
    fluid::BoundaryCondition bc = fluid::BoundaryCondition::channelFlowParabolic(U_max);

    // CSV writer
    fluid::CSVWriter writer(outputDir);
    writer.p_ref = bc.p_ref;
    writer.createOutputDirectory();
    writer.writeMetadata(grid);

    // Output callback
    int outputCount = 0;
    auto outputCallback = [&](int /*step*/, double time, const fluid::Grid& g) {
        writer.writeAll(g, outputCount, time);
        outputCount++;
    };

    // Run simulation
    auto startTime = std::chrono::high_resolution_clock::now();
    solver->run(grid, bc, endTime, 100, outputCallback);
    auto endTimePoint = std::chrono::high_resolution_clock::now();
    double wallTime = std::chrono::duration<double>(endTimePoint - startTime).count();

    std::cout << "\nSimulation completed!" << std::endl;
    std::cout << "Total steps: " << solver->getStepCount() << std::endl;
    std::cout << "Wall time: " << wallTime << " s" << std::endl;

    // Write simulation log
    writer.writeSimulationLog(solverType + "+" + limiterStr, grid, rho, nu, endTime,
                              solver->getStepCount(), wallTime, Re);

    return 0;
}
