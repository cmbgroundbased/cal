
# Name of loadable module / shared library
set(CAL_MOD cal)

# Library sources
set(CAL_SOURCES
    src/AATM_fun.cpp
    src/CALAtmSim.cpp
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

add_library(${CAL_MOD} SHARED ${CAL_SOURCES})

target_include_directories(${CAL_MOD} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src
)

if(MPI_FOUND)
    target_compile_definitions(${CAL_MOD} PRIVATE HAVE_MPI=1)
endif(MPI_FOUND)

# Dependencies

if(OpenMP_CXX_FOUND)
    target_compile_options(${CAL_MOD} PRIVATE "${OpenMP_CXX_FLAGS}")
    set_target_properties(${CAL_MOD} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    target_link_libraries(${CAL_MOD} "${OpenMP_CXX_LIBRARIES}")
endif(OpenMP_CXX_FOUND)

if(AATM_FOUND)
    target_compile_definitions(${CAL_MOD} PRIVATE HAVE_AATM=1)
    target_include_directories(${CAL_MOD} PUBLIC "${AATM_INCLUDE_DIRS}")
    target_link_libraries(${CAL_MOD} "${AATM_LIBRARIES}")
endif(AATM_FOUND)

if(CHOLMOD_FOUND)
    target_compile_definitions(${CAL_MOD} PRIVATE HAVE_CHOLMOD=1)
    target_include_directories(${CAL_MOD} PUBLIC "${CHOLMOD_INCLUDE_DIR}")
    target_link_libraries(${CAL_MOD} "${CHOLMOD_LIBRARY}")
endif(CHOLMOD_FOUND)

if(FFTW_FOUND)
    target_compile_definitions(${CAL_MOD} PRIVATE HAVE_FFTW=1)
    target_include_directories(${CAL_MOD} PUBLIC "${FFTW_INCLUDE_DIRS}")
    target_link_libraries(${CAL_MOD} "${FFTW_LIBRARIES}")
    if(FFTW_DOUBLE_THREADS_LIB_FOUND)
        target_compile_definitions(${CAL_MOD} PRIVATE HAVE_FFTW_THREADS=1)
    endif(FFTW_DOUBLE_THREADS_LIB_FOUND)
endif(FFTW_FOUND)

if(LAPACK_FOUND)
    target_compile_definitions(${CAL_MOD} PRIVATE HAVE_LAPACK=1)
    target_compile_definitions(${CAL_MOD} PRIVATE "LAPACK_NAMES_${LAPACK_NAMES}")
    set_target_properties(${CAL_MOD} PROPERTIES LINK_FLAGS
        "${LAPACK_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}"
    )
    target_link_libraries(${CAL_MOD} "${LAPACK_LIBRARIES}")
    target_link_libraries(${CAL_MOD} "${BLAS_LIBRARIES}")
endif(LAPACK_FOUND)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS ${CAL_MOD} DESTINATION ${CMAKE_INSTALL_LIBDIR})
