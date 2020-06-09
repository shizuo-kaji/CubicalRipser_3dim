cmake_minimum_required(VERSION 2.8.12)
project(cripser CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3 -march=native")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories("src/")

ADD_LIBRARY(mylib STATIC src/compute_pairs.cpp src/coboundary_enumerator.cpp src/cube.cpp src/dense_cubical_grids.cpp src/joint_pairs.cpp)

add_subdirectory(pybind11)

pybind11_add_module(cripser src/cubicalripser_pybind.cpp)
target_link_libraries(cripser PRIVATE mylib)

add_executable(cubicalripser src/cubicalripser.cpp)
target_link_libraries(cubicalripser PRIVATE mylib)
