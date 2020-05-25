
# Name of the internal static library
set(MPI_CAL_LIB cal_mpi_)

# Name of loadable module / shared library
set(MPI_CAL_MOD cal_mpi)

# Library sources
set(MPI_CAL_SOURCES
    src/AATM_fun.cpp
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
    src/math_rng.cpp
    src/math_sf.cpp
    src/observe.cpp
    src/print.cpp
    src/simulation.cpp
    src/smoothing_kernel.cpp
    src/smooth_interpolation.cpp
    src/sys_env.cpp
    src/sys_utils.cpp
)

# Add the internal library target

add_library(${MPI_CAL_LIB} STATIC ${MPI_CAL_SOURCES})

target_include_directories(${MPI_CAL_LIB} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src
)

if(MPI_FOUND)
    target_compile_definitions(${MPI_CAL_LIB} PRIVATE HAVE_MPI=1)
endif(MPI_FOUND)

# Dependencies

if(OpenMP_CXX_FOUND)
    target_compile_options(${MPI_CAL_LIB} PRIVATE "${OpenMP_CXX_FLAGS}")
    set_target_properties(${MPI_CAL_LIB} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    target_link_libraries(${MPI_CAL_LIB} "${OpenMP_CXX_LIBRARIES}")
endif(OpenMP_CXX_FOUND)

if(AATM_FOUND)
    target_compile_definitions(${MPI_CAL_LIB} PRIVATE HAVE_AATM=1)
    target_include_directories(${MPI_CAL_LIB} PUBLIC "${AATM_INCLUDE_DIRS}")
    target_link_libraries(${MPI_CAL_LIB} "${AATM_LIBRARIES}")
endif(AATM_FOUND)

if(CHOLMOD_FOUND)
    target_compile_definitions(${MPI_CAL_LIB} PRIVATE HAVE_CHOLMOD=1)
    target_include_directories(${MPI_CAL_LIB} PUBLIC "${CHOLMOD_INCLUDE_DIR}")
    target_link_libraries(${MPI_CAL_LIB} "${CHOLMOD_LIBRARY}")
endif(CHOLMOD_FOUND)

if(FFTW_FOUND)
    target_compile_definitions(${MPI_CAL_LIB} PRIVATE HAVE_FFTW=1)
    target_include_directories(${MPI_CAL_LIB} PUBLIC "${FFTW_INCLUDE_DIRS}")
    target_link_libraries(${MPI_CAL_LIB} "${FFTW_LIBRARIES}")
    if(FFTW_DOUBLE_THREADS_LIB_FOUND)
        target_compile_definitions(${MPI_CAL_LIB} PRIVATE HAVE_FFTW_THREADS=1)
    endif(FFTW_DOUBLE_THREADS_LIB_FOUND)
endif(FFTW_FOUND)

if(LAPACK_FOUND)
    target_compile_definitions(${MPI_CAL_LIB} PRIVATE HAVE_LAPACK=1)
    target_compile_definitions(${MPI_CAL_LIB} PRIVATE "LAPACK_NAMES_${LAPACK_NAMES}")
    set_target_properties(${MPI_CAL_LIB} PROPERTIES LINK_FLAGS
        "${LAPACK_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}"
    )
    target_link_libraries(${MPI_CAL_LIB} "${LAPACK_LIBRARIES}")
    target_link_libraries(${MPI_CAL_LIB} "${BLAS_LIBRARIES}")
endif(LAPACK_FOUND)



# Now define rules for the shared library

add_library(${MPI_CAL_MOD} SHARED ${MPI_CAL_SOURCES})

target_include_directories(${MPI_CAL_MOD} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS ${MPI_CAL_MOD} DESTINATION ${CMAKE_INSTALL_LIBDIR})
