cmake_minimum_required(VERSION 2.8.12)
project(cripser CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3 -Wall")

add_subdirectory(pybind11)
pybind11_add_module(cripser cubicalripser_pybind.cpp compute_pairs.cpp coboundary_enumerator.cpp cube.cpp dense_cubical_grids.cpp joint_pairs.cpp)
