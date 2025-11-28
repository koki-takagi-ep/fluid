cmake_minimum_required(VERSION 3.16)
project(IncompressibleFluidSolver VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# コンパイルオプション
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -O3)
endif()

# OpenMPの検出（オプション）
find_package(OpenMP)

# ヘッダファイル
include_directories(${CMAKE_SOURCE_DIR}/include)

# ライブラリ
add_library(fluid_solver
    src/Grid.cpp
    src/SolverBase.cpp
    src/Solver.cpp
    src/SimpleSolver.cpp
    src/PressureSolver.cpp
    src/BoundaryCondition.cpp
    src/CSVWriter.cpp
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(fluid_solver PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(fluid_solver PUBLIC USE_OPENMP)
endif()

# サンプルプログラム
add_executable(cavity_flow examples/cavity_flow.cpp)
target_link_libraries(cavity_flow fluid_solver)

add_executable(channel_flow examples/channel_flow.cpp)
target_link_libraries(channel_flow fluid_solver)

add_executable(channel_flow_simple examples/channel_flow_simple.cpp)
target_link_libraries(channel_flow_simple fluid_solver)

add_executable(cavity_flow_simple examples/cavity_flow_simple.cpp)
target_link_libraries(cavity_flow_simple fluid_solver)

# インストール
install(TARGETS fluid_solver LIBRARY DESTINATION lib)
install(DIRECTORY include/ DESTINATION include)
