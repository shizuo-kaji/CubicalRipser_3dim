cmake_minimum_required(VERSION 3.5)
cmake_policy(SET CMP0028 NEW)

project(cripser CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# pybind11 パッケージを検索して必要に応じてインクルードする
find_package(pybind11 REQUIRED)

include_directories("src/")

# Add pybind11 subdirectory
add_subdirectory(C:/Users/shing/CubicalRipser_3dim/pybind11)

# Add libraries
add_library(mylib STATIC src/compute_pairs.cpp src/cube.cpp src/joint_pairs.cpp)
add_library(vmylib STATIC src/coboundary_enumerator.cpp src/dense_cubical_grids.cpp)
add_library(tmylib STATIC src/coboundary_enumerator_T.cpp src/dense_cubical_grids_T.cpp)

# Add pybind11 modules
pybind11_add_module(cripser src/cubicalripser_pybind.cpp)
target_link_libraries(cripser PRIVATE mylib vmylib)

pybind11_add_module(tcripser src/cubicalripser_pybind_T.cpp)
target_link_libraries(tcripser PRIVATE mylib tmylib)

# Add executables
add_executable(cubicalripser src/cubicalripser.cpp)
target_link_libraries(cubicalripser PRIVATE mylib vmylib)

add_executable(tcubicalripser src/cubicalripser.cpp)
target_link_libraries(tcubicalripser PRIVATE mylib tmylib)

# pybind11のインクルードディレクトリを指定
include_directories("C:/Users/shing/CubicalRipser_3dim/pybind11/include")
