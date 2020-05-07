
# Name of the internal static library
set(CAL_LIB cal_)

# Name of loadable module / shared library
set(CAL_MOD cal)

# Library sources
set(CAL_SOURCES
    src/CALAtmSim.cpp
    src/AATM_fun.cpp
    src/sys_env.cpp
    src/sys_utils.cpp
    src/math_rng.cpp
    src/math_sf.cpp
)

# Add the internal library target

add_library(${CAL_LIB} STATIC ${CAL_SOURCES})

target_include_directories(${CAL_LIB} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src
)

if(MPI_FOUND)
    target_compile_definitions(${CAL_LIB} PRIVATE HAVE_MPI=1)
endif(MPI_FOUND)

# Dependencies

if(OpenMP_CXX_FOUND)
    target_compile_options(${CAL_LIB} PRIVATE "${OpenMP_CXX_FLAGS}")
    set_target_properties(${CAL_LIB} PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    target_link_libraries(${CAL_LIB} "${OpenMP_CXX_LIBRARIES}")
endif(OpenMP_CXX_FOUND)

if(AATM_FOUND)
    target_compile_definitions(${CAL_LIB} PRIVATE HAVE_AATM=1)
    target_include_directories(${CAL_LIB} PUBLIC "${AATM_INCLUDE_DIRS}")
    target_link_libraries(${CAL_LIB} "${AATM_LIBRARIES}")
endif(AATM_FOUND)

if(CHOLMOD_FOUND)
    target_compile_definitions(${CAL_LIB} PRIVATE HAVE_CHOLMOD=1)
    target_include_directories(${CAL_LIB} PUBLIC "${CHOLMOD_INCLUDE_DIR}")
    target_link_libraries(${CAL_LIB} "${CHOLMOD_LIBRARY}")
endif(CHOLMOD_FOUND)

if(FFTW_FOUND)
    target_compile_definitions(${CAL_LIB} PRIVATE HAVE_FFTW=1)
    target_include_directories(${CAL_LIB} PUBLIC "${FFTW_INCLUDE_DIRS}")
    target_link_libraries(${CAL_LIB} "${FFTW_LIBRARIES}")
    if(FFTW_DOUBLE_THREADS_LIB_FOUND)
        target_compile_definitions(${CAL_LIB} PRIVATE HAVE_FFTW_THREADS=1)
    endif(FFTW_DOUBLE_THREADS_LIB_FOUND)
endif(FFTW_FOUND)

if(LAPACK_FOUND)
    target_compile_definitions(${CAL_LIB} PRIVATE HAVE_LAPACK=1)
    target_compile_definitions(${CAL_LIB} PRIVATE "LAPACK_NAMES_${LAPACK_NAMES}")
    set_target_properties(${CAL_LIB} PROPERTIES LINK_FLAGS
        "${LAPACK_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}"
    )
    target_link_libraries(${CAL_LIB} "${LAPACK_LIBRARIES}")
    target_link_libraries(${CAL_LIB} "${BLAS_LIBRARIES}")
endif(LAPACK_FOUND)



# Now define rules for the shared library

add_library(${CAL_MOD} SHARED ${CAL_SOURCES})

target_include_directories(${CAL_MOD} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE .
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS ${CAL_MOD} DESTINATION ${CMAKE_INSTALL_LIBDIR})
