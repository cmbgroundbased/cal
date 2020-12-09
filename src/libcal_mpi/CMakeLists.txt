string(STRIP "${MPI_CXX_COMPILE_FLAGS}" mpi_comp_flags)
string(STRIP "${MPI_CXX_LINK_FLAGS}" mpi_link_flags)

# Name of loadable module / shared library
set(MPI_CAL_MOD cal_mpi_core)

# Library sources
set(MPI_CAL_SOURCES
    src/CAL_MPI_AtmSim.cpp
    src/compress_volume.cpp
    src/coord_transform.cpp
    src/covariance_apply.cpp
    src/covariance_build.cpp
    src/covariance_eval.cpp
    src/covariance_square.cpp
    src/draw.cpp
    src/get_slice.cpp
    src/get_volume.cpp
    src/in_cone.cpp
    src/kolmogorov_init.cpp
    src/kolmovorov_autocov.cpp
    src/load_save_realization.cpp
    src/mpi_init.cpp
    src/observe.cpp
    src/print.cpp
    src/simulation.cpp
    src/smoothing_kernel.cpp
    src/smooth_interpolation.cpp
)

add_library(${MPI_CAL_MOD} ${MPI_CAL_SOURCES})

target_include_directories(${MPI_CAL_MOD} BEFORE PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../libcal"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libcal/include"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

# Dependencies

if(OpenMP_CXX_FOUND)
    target_compile_options(${MPI_CAL_MOD} PRIVATE "${OpenMP_CXX_FLAGS}")
    set_target_properties(${MPI_CAL_MOD} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    target_link_libraries(${MPI_CAL_MOD} "${OpenMP_CXX_LIBRARIES}")
endif(OpenMP_CXX_FOUND)

if(AATM_FOUND)
    target_compile_definitions(${MPI_CAL_MOD} PRIVATE HAVE_AATM=1)
    target_include_directories(${MPI_CAL_MOD} PUBLIC "${AATM_INCLUDE_DIRS}")
    target_link_libraries(${MPI_CAL_MOD} "${AATM_LIBRARIES}")
endif(AATM_FOUND)

if(CHOLMOD_FOUND)
    target_compile_definitions(${MPI_CAL_MOD} PRIVATE HAVE_CHOLMOD=1)
    target_include_directories(${MPI_CAL_MOD} PUBLIC "${CHOLMOD_INCLUDE_DIR}")
    target_link_libraries(${MPI_CAL_MOD} "${CHOLMOD_LIBRARY}")
endif(CHOLMOD_FOUND)

target_include_directories(${MPI_CAL_MOD} PRIVATE "${MPI_CXX_INCLUDE_PATH}")
target_compile_options(${MPI_CAL_MOD} PRIVATE "${mpi_comp_flags}")
set_target_properties(${MPI_CAL_MOD} PROPERTIES LINK_FLAGS "${mpi_link_flags}")
target_link_libraries(${MPI_CAL_MOD} "${MPI_CXX_LIBRARIES}")

# Build the total library

set(CAL_MPI_MOD_TOT cal_mpi)
add_library(${CAL_MPI_MOD_TOT} cal_mpi.cpp)


target_include_directories(${CAL_MPI_MOD_TOT} BEFORE PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_include_directories(${CAL_MPI_MOD_TOT} BEFORE PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../libcal/include"
)

target_include_directories(${CAL_MPI_MOD_TOT} PRIVATE "${MPI_CXX_INCLUDE_PATH}")
target_compile_options(${CAL_MPI_MOD_TOT} PRIVATE "${mpi_comp_flags}")

if(OpenMP_CXX_FOUND)
    target_compile_options(${CAL_MPI_MOD_TOT} PRIVATE "${OpenMP_CXX_FLAGS}")
    set_target_properties(${CAL_MPI_MOD_TOT} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
endif(OpenMP_CXX_FOUND)

target_link_libraries(${CAL_MPI_MOD_TOT} ${MPI_CAL_MOD} ${CAL_MOD_TOT} "${MPI_CXX_LIBRARIES}")

# Public headers
install(DIRECTORY include/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(TARGETS ${CAL_MPI_MOD_TOT} DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

# MPI unit test executable

# etc ...