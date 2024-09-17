cmake_minimum_required(VERSION 3.5)
cmake_policy(SET CMP0028 NEW)

project(cripser CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_PREFIX_PATH "C:/Users/shing/Documents/CubicalRipser_3dim")
set(pybind11_DIR "C:/Users/shing/AppData/Local/Packages/PythonSoftwareFoundation.Python.3.11_qbz5n2kfra8p0/LocalCache/local-packages/Python311/site-packages/pybind11/share/cmake")

# Find pybind11 package
find_package(pybind11 CONFIG REQUIRED)

# Add libraries
add_library(mylib STATIC src/compute_pairs.cpp src/cube.cpp src/joint_pairs.cpp)
add_library(vmylib STATIC src/coboundary_enumerator.cpp src/dense_cubical_grids.cpp)
add_library(tmylib STATIC src/coboundary_enumerator_T.cpp src/dense_cubical_grids_T.cpp)

# Add pybind11 modules
pybind11_add_module(cripser src/cubicalripser_pybind.cpp)
target_link_libraries(cripser PRIVATE mylib vmylib pybind11::module)

pybind11_add_module(tcripser src/cubicalripser_pybind_T.cpp)
target_link_libraries(tcripser PRIVATE mylib tmylib pybind11::module)

# Add executables
add_executable(cubicalripser src/cubicalripser.cpp)
target_link_libraries(cubicalripser PRIVATE mylib vmylib)

add_executable(tcubicalripser src/cubicalripser.cpp)
target_link_libraries(tcubicalripser PRIVATE mylib tmylib)

# Include pybind11 directory
include_directories(${pybind11_INCLUDE_DIRS})
