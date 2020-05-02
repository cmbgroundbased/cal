cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

if (MSVC)
    add_compile_options(/W3 /WX)
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
find_package(OpenMP)
find_package(PkgConfig)
find_package(BLAS)
if(BLAS_FOUND)
    find_package(LAPACK)
    if(LAPACK_FOUND)
        find_package(LAPACKnames)
    else(LAPACK_FOUND)
        message(SEND_ERROR "Could not find a working LAPACK installation")
    endif(LAPACK_FOUND)
endif(BLAS_FOUND)

find_package(FFTW)
find_package(AATM)

enable_testing()
add_subdirectory(src)
