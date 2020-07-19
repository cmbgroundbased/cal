# Name of the internal static library
set(CAL_LIB cal_)

#======================================================

add_executable(observe_test
    observe_test.cpp
)
add_executable(vector_test
    vector_test.cpp
)

if(WIN32)
    target_compile_definitions(aatm_test_basic PRIVATE HAVE_WINDOWS=1)
endif(WIN32)

target_include_directories(observe_test PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libcal/src"
)
target_include_directories(vector_test PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/../libcal/src"
)


target_link_libraries(observe_test ${CAL_LIB})
target_link_libraries(vector_test ${CAL_LIB})

# install(TARGETS aatm_test_basic DESTINATION ${CMAKE_INSTALL_BINDIR})

# add_test(NAME test_basic COMMAND aatm_test_basic)
