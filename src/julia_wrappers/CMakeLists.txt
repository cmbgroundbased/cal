set(JL_MODULE AtmSimJulia)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# find_package(JlCxx)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")

message(STATUS "Found JlCxx at ${JlCxx_location}")

set(JL_WR_SRC
    src/AtmSim_julia.cpp
)

add_library(${JL_MODULE} SHARED ${JL_WR_SRC})

target_include_directories(${JL_MODULE} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../libcal/include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(${JL_MODULE} cal JlCxx::cxxwrap_julia)

install(
TARGETS ${JL_MODULE}
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION lib
)
