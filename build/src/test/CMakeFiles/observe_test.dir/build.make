# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/algebrato/Progetti/CMB4G/libcal

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/algebrato/Progetti/CMB4G/libcal/build

# Include any dependencies generated for this target.
include src/test/CMakeFiles/observe_test.dir/depend.make

# Include the progress variables for this target.
include src/test/CMakeFiles/observe_test.dir/progress.make

# Include the compile flags for this target's objects.
include src/test/CMakeFiles/observe_test.dir/flags.make

src/test/CMakeFiles/observe_test.dir/observe_test.cpp.o: src/test/CMakeFiles/observe_test.dir/flags.make
src/test/CMakeFiles/observe_test.dir/observe_test.cpp.o: ../src/test/observe_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/algebrato/Progetti/CMB4G/libcal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/test/CMakeFiles/observe_test.dir/observe_test.cpp.o"
	cd /home/algebrato/Progetti/CMB4G/libcal/build/src/test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/observe_test.dir/observe_test.cpp.o -c /home/algebrato/Progetti/CMB4G/libcal/src/test/observe_test.cpp

src/test/CMakeFiles/observe_test.dir/observe_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/observe_test.dir/observe_test.cpp.i"
	cd /home/algebrato/Progetti/CMB4G/libcal/build/src/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/algebrato/Progetti/CMB4G/libcal/src/test/observe_test.cpp > CMakeFiles/observe_test.dir/observe_test.cpp.i

src/test/CMakeFiles/observe_test.dir/observe_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/observe_test.dir/observe_test.cpp.s"
	cd /home/algebrato/Progetti/CMB4G/libcal/build/src/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/algebrato/Progetti/CMB4G/libcal/src/test/observe_test.cpp -o CMakeFiles/observe_test.dir/observe_test.cpp.s

# Object files for target observe_test
observe_test_OBJECTS = \
"CMakeFiles/observe_test.dir/observe_test.cpp.o"

# External object files for target observe_test
observe_test_EXTERNAL_OBJECTS =

src/test/observe_test: src/test/CMakeFiles/observe_test.dir/observe_test.cpp.o
src/test/observe_test: src/test/CMakeFiles/observe_test.dir/build.make
src/test/observe_test: src/libcal/libcal_.a
src/test/observe_test: /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.so
src/test/observe_test: /usr/lib/x86_64-linux-gnu/libpthread.so
src/test/observe_test: /usr/local/lib/libaatm.so
src/test/observe_test: /usr/local/lib/libcholmod.so
src/test/observe_test: /usr/lib/x86_64-linux-gnu/libfftw3.so
src/test/observe_test: /usr/lib/x86_64-linux-gnu/libfftw3f.so
src/test/observe_test: /usr/lib/x86_64-linux-gnu/libfftw3l.so
src/test/observe_test: /usr/lib/x86_64-linux-gnu/libfftw3_threads.so
src/test/observe_test: /usr/lib/x86_64-linux-gnu/libfftw3f_threads.so
src/test/observe_test: /usr/lib/x86_64-linux-gnu/libfftw3l_threads.so
src/test/observe_test: /usr/lib/x86_64-linux-gnu/libfftw3_omp.so
src/test/observe_test: /usr/lib/x86_64-linux-gnu/libfftw3f_omp.so
src/test/observe_test: /usr/lib/x86_64-linux-gnu/libfftw3l_omp.so
src/test/observe_test: /usr/lib/x86_64-linux-gnu/libopenblas.so
src/test/observe_test: src/test/CMakeFiles/observe_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/algebrato/Progetti/CMB4G/libcal/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable observe_test"
	cd /home/algebrato/Progetti/CMB4G/libcal/build/src/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/observe_test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/test/CMakeFiles/observe_test.dir/build: src/test/observe_test

.PHONY : src/test/CMakeFiles/observe_test.dir/build

src/test/CMakeFiles/observe_test.dir/clean:
	cd /home/algebrato/Progetti/CMB4G/libcal/build/src/test && $(CMAKE_COMMAND) -P CMakeFiles/observe_test.dir/cmake_clean.cmake
.PHONY : src/test/CMakeFiles/observe_test.dir/clean

src/test/CMakeFiles/observe_test.dir/depend:
	cd /home/algebrato/Progetti/CMB4G/libcal/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/algebrato/Progetti/CMB4G/libcal /home/algebrato/Progetti/CMB4G/libcal/src/test /home/algebrato/Progetti/CMB4G/libcal/build /home/algebrato/Progetti/CMB4G/libcal/build/src/test /home/algebrato/Progetti/CMB4G/libcal/build/src/test/CMakeFiles/observe_test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/test/CMakeFiles/observe_test.dir/depend
