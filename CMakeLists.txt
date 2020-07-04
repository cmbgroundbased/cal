cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

if (MSVC)
    add_compile_options(/std=c++11)
else()
    add_compile_options(-std=c++11 -Wextra -pedantic)
endif()

project(cal VERSION 1.0.0 LANGUAGES C CXX)

# Auxiliary files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Use GNUInstallDirs to install the library in the correct locations
# It's is always true for all the platforms
include(GNUInstallDirs)

# Now, the defauld is: Debug. Convert it in Release when the repo
# will be ready to be released.
include(BuildType)

# The code should be build with PIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# PkgConfig - is used to retrieve information about installed libraries in the system
# For more details: `man pkg-config`
find_package(PkgConfig REQUIRED)

# Find OpenMP and MPI
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

# Find BLAS library with LAPACK extension
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
    find_package(LAPACK REQUIRED)
    if(LAPACK_FOUND)
        find_package(LAPACKnames REQUIRED)
    else(LAPACK_FOUND)
        message(SEND_ERROR "Could not find a working LAPACK installation")
    endif(LAPACK_FOUND)
endif(BLAS_FOUND)

# FFTW to perform Fourier Transorm
find_package(FFTW REQUIRED)

# Now, the cal baseline use ALMA - Atacama Large Millimiter Array
# by Instituto de Estructura de la Materia, 2009
find_package(AATM REQUIRED)

# The SuiteSparse library is used by the cal build and apply the
# atmospheric covariance matrix
find_package(SuiteSparse REQUIRED)

# Tests - work in progress
enable_testing()

# The directory that gather the whole source code.
add_subdirectory(src)
