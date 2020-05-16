cmake_minimum_required(VERSION 2.8.12)
project(cubicalripser CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3 -Wall")

add_executable(cubicalripser cubicalripser.cpp compute_pairs.cpp coboundary_enumerator.cpp cube.cpp dense_cubical_grids.cpp joint_pairs.cpp)
