
# Name of the internal static library
set(CAL_LIB cal_)

# Name of loadable module / shared library
set(CAL_MOD cal)

# Library sources
set(CAL_SOURCES
    src/CALAtmSim.cpp
    src/sys_env.cpp
    src/sys_utils.cpp
)

# Add the internal library target

add_library(${CAL_LIB} STATIC ${CAL_SOURCES})

target_include_directories(${CAL_LIB} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE src
)

# Now define rules for the shared library

add_library(${CAL_MOD} SHARED ${CAL_SOURCES})

target_include_directories(${CAL_MOD} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE .
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS ${CAL_MOD} DESTINATION ${CMAKE_INSTALL_LIBDIR})
