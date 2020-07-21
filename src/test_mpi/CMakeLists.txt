# Name of the internal static library
set(MPI_CAL_LIB cal_mpi)

#======================================================

add_executable(atm_dummy_simulation
    Atmosphere_construction_MPI.cpp
)

target_include_directories(atm_dummy_simulation PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libcal_mpi"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libcal_mpi/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libcal_mpi/include"
)
target_link_libraries(atm_dummy_simulation ${MPI_CAL_LIB} mpi)

# install(TARGETS aatm_test_basic DESTINATION ${CMAKE_INSTALL_BINDIR})

#add_test(NAME test_basic COMMAND aatm_test_basic)
