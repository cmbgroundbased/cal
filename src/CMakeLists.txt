
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library directory - main and MPI version
# if MPI_FOUND is True
add_subdirectory(libcal)

if(MPI_FOUND)
    add_subdirectory(libcal_mpi)
endif()

# Tests subdirs
# add_subdirectory(test)

# add_subdirectory(test_mpi)

# python wrappers for CAL
add_subdirectory(pycal)

# Julia wappers for CAL (decomment when the julia-wraps will be written)
#if(JULIA)
#    message("-- Julia Wappers are enabled")
#    add_subdirectory(julia_wrappers)
#endif()
