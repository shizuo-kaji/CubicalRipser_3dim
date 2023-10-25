cmake_minimum_required(VERSION 3.5)
cmake_policy(SET CMP0028 NEW)

project(cripser CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include_directories("src/")


ADD_LIBRARY(mylib STATIC src/compute_pairs.cpp src/cube.cpp src/joint_pairs.cpp)

## V-construction
ADD_LIBRARY(vmylib STATIC src/coboundary_enumerator.cpp src/dense_cubical_grids.cpp)

## T-construction
ADD_LIBRARY(tmylib STATIC src/coboundary_enumerator_T.cpp src/dense_cubical_grids_T.cpp)

## Pytno modules
add_subdirectory(pybind11)
pybind11_add_module(cripser src/cubicalripser_pybind.cpp)
target_link_libraries(cripser PRIVATE mylib)
target_link_libraries(cripser PRIVATE vmylib)
pybind11_add_module(tcripser src/cubicalripser_pybind_T.cpp)
target_link_libraries(tcripser PRIVATE mylib)
target_link_libraries(tcripser PRIVATE tmylib)

## command-line executables
add_executable(cubicalripser src/cubicalripser.cpp)
add_executable(tcubicalripser src/cubicalripser.cpp)
target_link_libraries(cubicalripser PRIVATE mylib)
target_link_libraries(cubicalripser PRIVATE vmylib)
target_link_libraries(tcubicalripser PRIVATE mylib)
target_link_libraries(tcubicalripser PRIVATE tmylib)
