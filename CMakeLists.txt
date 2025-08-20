cmake_minimum_required(VERSION 3.15)  # Updated to a more recent version for better feature support
cmake_policy(SET CMP0028 NEW)

# Project settings
project(cripser LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_OSX_ARCHITECTURES arm64;x86_64)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
#find_package(pybind11 CONFIG REQUIRED)


# Include directories
include_directories("src/")

# Create static libraries
add_library(mylib STATIC
    src/compute_pairs.cpp
    src/coboundary_enumerator.cpp
    src/joint_pairs.cpp
)

# V-construction library
add_library(vmylib STATIC
    src/dense_cubical_grids.cpp
)

# T-construction library
add_library(tmylib STATIC
    src/dense_cubical_grids_T.cpp
)

# Python modules using pybind11
add_subdirectory(pybind11)

pybind11_add_module(_cripser src/cubicalripser_pybind.cpp)
target_link_libraries(_cripser PRIVATE mylib vmylib)

pybind11_add_module(tcripser src/cubicalripser_pybind_T.cpp)
target_link_libraries(tcripser PRIVATE mylib tmylib)

# Inject package version from setup.py if provided
if(DEFINED CRIPSER_VERSION)
  message(STATUS "Injecting version ${CRIPSER_VERSION} into pybind11 modules")
  target_compile_definitions(_cripser PRIVATE VERSION_INFO=\"${CRIPSER_VERSION}\")
  target_compile_definitions(tcripser PRIVATE VERSION_INFO=\"${CRIPSER_VERSION}\")
endif()

# Command-line executables
add_executable(cubicalripser src/cubicalripser.cpp)
target_link_libraries(cubicalripser PRIVATE mylib vmylib)

add_executable(tcubicalripser src/cubicalripser.cpp)
target_link_libraries(tcubicalripser PRIVATE mylib tmylib)

# Optional: Group libraries and executables in folders in IDEs
set_target_properties(mylib vmylib tmylib _cripser tcripser cubicalripser tcubicalripser PROPERTIES FOLDER "cripser")
