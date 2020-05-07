cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

if (MSVC)
    add_compile_options(/std=c++11)
else()
    add_compile_options(-std=c++11)
endif()

project(cal VERSION 1.0.0 LANGUAGES C CXX)

# Auxiliary files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Use GNUInstallDirs to install the library in the correct locations
# It's is always true for all the platforms
include(GNUInstallDirs)
#include(BuildType)
# The code should be build with PIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# External packages
find_package(OpenMP REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
    find_package(LAPACK REQUIRED)
    if(LAPACK_FOUND)
        find_package(LAPACKnames REQUIRED)
    else(LAPACK_FOUND)
        message(SEND_ERROR "Could not find a working LAPACK installation")
    endif(LAPACK_FOUND)
endif(BLAS_FOUND)

find_package(FFTW REQUIRED)
find_package(AATM REQUIRED)
find_package(SuiteSparse REQUIRED)
#find_package(RANDOM123 REQUIRED)

enable_testing()
add_subdirectory(src)
