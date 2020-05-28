# Name of the internal static library
set(MPI_CAL_LIB cal_mpi_)

#======================================================

add_executable(mpi_atm_ctor
    Atmosphere_construction_MPI.cpp
)

add_executable(mpi_atm_obs
    Atmosphere_observations_MPI.cpp
)

target_include_directories(mpi_atm_ctor PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libcal_mpi/src"
)
target_include_directories(mpi_atm_obs PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libcal_mpi/src"
)

target_link_libraries(mpi_atm_ctor ${MPI_CAL_LIB} mpi)

target_link_libraries(mpi_atm_obs ${MPI_CAL_LIB})

# install(TARGETS aatm_test_basic DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_basic COMMAND aatm_test_basic)
