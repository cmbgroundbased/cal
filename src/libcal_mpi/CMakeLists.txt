string(STRIP "${MPI_CXX_COMPILE_FLAGS}" mpi_comp_flags)
string(STRIP "${MPI_CXX_LINK_FLAGS}" mpi_link_flags)

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
    src/mpi_init.cpp
)

# Add the internal library target
add_library(${MPI_CAL_LIB} STATIC ${MPI_CAL_SOURCES})
add_library(${MPI_CAL_MOD} SHARED ${MPI_CAL_SOURCES})

target_include_directories(${MPI_CAL_LIB} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src
)
target_include_directories(${MPI_CAL_MOD} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src
)

# Dependencies

if(OpenMP_CXX_FOUND)
    target_compile_options(${MPI_CAL_LIB} PRIVATE "${OpenMP_CXX_FLAGS}")
    target_compile_options(${MPI_CAL_MOD} PRIVATE "${OpenMP_CXX_FLAGS}")
    set_target_properties(${MPI_CAL_LIB} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    set_target_properties(${MPI_CAL_MOD} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    target_link_libraries(${MPI_CAL_LIB} "${OpenMP_CXX_LIBRARIES}")
    target_link_libraries(${MPI_CAL_MOD} "${OpenMP_CXX_LIBRARIES}")
endif(OpenMP_CXX_FOUND)

if(AATM_FOUND)
    target_compile_definitions(${MPI_CAL_LIB} PRIVATE HAVE_AATM=1)
    target_compile_definitions(${MPI_CAL_MOD} PRIVATE HAVE_AATM=1)
    target_include_directories(${MPI_CAL_LIB} PUBLIC "${AATM_INCLUDE_DIRS}")
    target_include_directories(${MPI_CAL_MOD} PUBLIC "${AATM_INCLUDE_DIRS}")
    target_link_libraries(${MPI_CAL_LIB} "${AATM_LIBRARIES}")
    target_link_libraries(${MPI_CAL_MOD} "${AATM_LIBRARIES}")
endif(AATM_FOUND)

if(CHOLMOD_FOUND)
    target_compile_definitions(${MPI_CAL_LIB} PRIVATE HAVE_CHOLMOD=1)
    target_compile_definitions(${MPI_CAL_MOD} PRIVATE HAVE_CHOLMOD=1)
    target_include_directories(${MPI_CAL_LIB} PUBLIC "${CHOLMOD_INCLUDE_DIR}")
    target_include_directories(${MPI_CAL_MOD} PUBLIC "${CHOLMOD_INCLUDE_DIR}")
    target_link_libraries(${MPI_CAL_LIB} "${CHOLMOD_LIBRARY}")
    target_link_libraries(${MPI_CAL_MOD} "${CHOLMOD_LIBRARY}")
endif(CHOLMOD_FOUND)

if(FFTW_FOUND)
    target_compile_definitions(${MPI_CAL_LIB} PRIVATE HAVE_FFTW=1)
    target_compile_definitions(${MPI_CAL_MOD} PRIVATE HAVE_FFTW=1)
    target_include_directories(${MPI_CAL_LIB} PUBLIC "${FFTW_INCLUDE_DIRS}")
    target_include_directories(${MPI_CAL_MOD} PUBLIC "${FFTW_INCLUDE_DIRS}")
    target_link_libraries(${MPI_CAL_LIB} "${FFTW_LIBRARIES}")
    target_link_libraries(${MPI_CAL_MOD} "${FFTW_LIBRARIES}")
    if(FFTW_DOUBLE_THREADS_LIB_FOUND)
        target_compile_definitions(${MPI_CAL_LIB} PRIVATE HAVE_FFTW_THREADS=1)
        target_compile_definitions(${MPI_CAL_MOD} PRIVATE HAVE_FFTW_THREADS=1)
    endif(FFTW_DOUBLE_THREADS_LIB_FOUND)
endif(FFTW_FOUND)

if(LAPACK_FOUND)
    target_compile_definitions(${MPI_CAL_LIB} PRIVATE HAVE_LAPACK=1)
    target_compile_definitions(${MPI_CAL_MOD} PRIVATE HAVE_LAPACK=1)
    target_compile_definitions(${MPI_CAL_LIB} PRIVATE "LAPACK_NAMES_${LAPACK_NAMES}")
    target_compile_definitions(${MPI_CAL_MOD} PRIVATE "LAPACK_NAMES_${LAPACK_NAMES}")
    set_target_properties(${MPI_CAL_LIB} PROPERTIES LINK_FLAGS
        "${LAPACK_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}"
    )
    set_target_properties(${MPI_CAL_MOD} PROPERTIES LINK_FLAGS
        "${LAPACK_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}"
    )
    target_link_libraries(${MPI_CAL_LIB} "${LAPACK_LIBRARIES}")
    target_link_libraries(${MPI_CAL_MOD} "${LAPACK_LIBRARIES}")

    target_link_libraries(${MPI_CAL_LIB} "${BLAS_LIBRARIES}")
    target_link_libraries(${MPI_CAL_MOD} "${BLAS_LIBRARIES}")
endif(LAPACK_FOUND)

# In principle, don't have to check if MPI_FOUND is true
# but we can...(evaluate to get rid of this statement)
if(MPI_FOUND)
    target_compile_definitions(${MPI_CAL_LIB} PRIVATE HAVE_MPI=1)
    target_compile_definitions(${MPI_CAL_MOD} PRIVATE HAVE_MPI=1)
    target_include_directories(${MPI_CAL_LIB} PRIVATE "${MPI_CXX_INCLUDE_PATH}")
    target_include_directories(${MPI_CAL_MOD} PRIVATE "${MPI_CXX_INCLUDE_PATH}")
    target_link_libraries(${MPI_CAL_LIB} PRIVATE ${MPI_C_LIBRARIES})
    target_link_libraries(${MPI_CAL_MOD} PRIVATE ${MPI_C_LIBRARIES})
    #set_target_properties(${MPI_CAL_LIB} PRIVATE "${mpi_comp_flags}")
    # set_target_properties(${MPI_CAL_MOD} PRIVATE "${mpi_comp_flags}")
    # set_target_properties(${MPI_CAL_LIB} PROPERTIES LINK_FLAGS "${mpi_link_flags}")
    # set_target_properties(${MPI_CAL_MOD} PROPERTIES LINK_FLAGS "${mpi_link_flags}")
endif()

# install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# install(TARGETS ${MPI_CAL_MOD} DESTINATION ${CMAKE_INSTALL_LIBDIR})
