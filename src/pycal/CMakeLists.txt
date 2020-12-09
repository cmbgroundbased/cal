# pybind11
add_subdirectory(pybind11)

# Copy the RELEASE file into the build directory
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/RELEASE"
    "${CMAKE_CURRENT_BINARY_DIR}/RELEASE"
    COPYONLY
)

set(MPI_CAL_LIB cal_mpi)
set(CAL_LIB cal)

# Determine our python install location and set up a rule to install the
# modules

execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
        import sys
        print('{}.{}'.format(sys.version_info.major, sys.version_info.minor))"
    OUTPUT_VARIABLE PYTHON_MAJORMINOR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "import site; print(site.getsitepackages()[0])"
    OUTPUT_VARIABLE PYTHON_SITE
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

message("-- Python exe:  " ${PYTHON_EXECUTABLE})
message("-- Python Path: " ${PYTHON_SITE})

pybind11_add_module(_libcal MODULE
    _libcal_common.cpp
    _libcal.cpp
    _libcal_sys.cpp
    _libcal_math_sf.cpp
    _libcal_math_rng.cpp
    _libcal_math_qarray.cpp
    _libcal_atm.cpp
)

if(OpenMP_CXX_FOUND)
    target_compile_options(_libcal PRIVATE "${OpenMP_CXX_FLAGS}")
    set_target_properties(_libcal PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    set(OMP_LIBS "${OpenMP_CXX_LIBRARIES}")
else()
    set(OMP_LIBS "")
endif(OpenMP_CXX_FOUND)

if(AATM_FOUND)
    target_compile_definitions(_libcal PRIVATE HAVE_AATM=1)
    target_include_directories(_libcal PUBLIC "${AATM_INCLUDE_DIRS}")
endif(AATM_FOUND)

if(CHOLMOD_FOUND)
    target_compile_definitions(_libcal PRIVATE HAVE_CHOLMOD=1)
    target_include_directories(_libcal PUBLIC "${CHOLMOD_INCLUDE_DIR}")
endif(CHOLMOD_FOUND)

target_link_libraries(_libcal PRIVATE ${CAL_LIB} "${OMP_LIBS}")

# Include path to the toast headers
target_include_directories(_libcal BEFORE PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libcal/include"
)

install(TARGETS _libcal DESTINATION ${PYTHON_SITE}/pycal)

# Install all the other python files

install(FILES
    __init__.py
    utils.py
    mpi.py
    timing.py
    cache.py
    dist.py
    op.py
    vis.py
    rng.py
    weather.py
    qarray.py
    RELEASE
    DESTINATION ${PYTHON_SITE}/pycal
)

# Optionally create a module for the MPI library
if(MPI_FOUND AND MPI4PY_FOUND)
    pybind11_add_module(_libcal_mpi MODULE
        _libcal_common.cpp
        _libcal_atm_mpi.cpp    
    )

    string(STRIP "${MPI_CXX_COMPILE_FLAGS}" mpi_comp_flags)
    target_compile_options(_libcal_mpi PRIVATE "${mpi_comp_flags}")

    string(STRIP "${MPI_CXX_LINK_FLAGS}" mpi_link_flags)
    set_target_properties(_libcal_mpi
        PROPERTIES LINK_FLAGS "${mpi_link_flags}"
    )

    target_include_directories(_libcal_mpi BEFORE PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${CMAKE_CURRENT_SOURCE_DIR}/../libcal/include"
        "${CMAKE_CURRENT_SOURCE_DIR}/../libcal_mpi/include"
        ${MPI_CXX_INCLUDE_PATH}
        ${MPI4PY_INCLUDE_DIR}
    )

    if(OpenMP_CXX_FOUND)
        target_compile_options(_libcal_mpi PRIVATE "${OpenMP_CXX_FLAGS}")
        set_target_properties(_libcal_mpi PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
    endif(OpenMP_CXX_FOUND)

    if(AATM_FOUND)
        target_compile_definitions(_libcal_mpi PRIVATE HAVE_AATM=1)
        target_include_directories(_libcal_mpi PUBLIC "${AATM_INCLUDE_DIRS}")
    endif(AATM_FOUND)

    if(CHOLMOD_FOUND)
        target_compile_definitions(_libcal_mpi PRIVATE HAVE_CHOLMOD=1)
        target_include_directories(_libcal_mpi PUBLIC "${CHOLMOD_INCLUDE_DIR}")
    endif(CHOLMOD_FOUND)

    target_link_libraries(_libcal_mpi
        PRIVATE libcal_mpi "${MPI_CXX_LIBRARIES}" "${OMP_LIBS}"
    )

    install(TARGETS _libcal_mpi DESTINATION ${PYTHON_SITE}/pycal)

endif(MPI_FOUND AND MPI4PY_FOUND)

# Process the sub directories
add_subdirectory(todmap)
