# Name of the internal static library
set(CAL_LIB cal_)

#======================================================

add_executable(cal_test_basic
    CalBasicTest.cpp
)

if(WIN32)
    target_compile_definitions(aatm_test_basic PRIVATE HAVE_WINDOWS=1)
endif(WIN32)

target_include_directories(cal_test_basic PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libcal/src"
)

target_link_libraries(cal_test_basic ${CAL_LIB})

# install(TARGETS aatm_test_basic DESTINATION ${CMAKE_INSTALL_BINDIR})

add_test(NAME test_basic COMMAND aatm_test_basic)
