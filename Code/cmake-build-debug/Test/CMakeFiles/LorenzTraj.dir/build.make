# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

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
CMAKE_COMMAND = /u/anmc/garegnan/Downloads/clion-2016.2.3/bin/cmake/bin/cmake

# The command to remove a file.
RM = /u/anmc/garegnan/Downloads/clion-2016.2.3/bin/cmake/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /u/anmc/garegnan/Desktop/Project/Code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug

# Include any dependencies generated for this target.
include Test/CMakeFiles/LorenzTraj.dir/depend.make

# Include the progress variables for this target.
include Test/CMakeFiles/LorenzTraj.dir/progress.make

# Include the compile flags for this target's objects.
include Test/CMakeFiles/LorenzTraj.dir/flags.make

Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o: Test/CMakeFiles/LorenzTraj.dir/flags.make
Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o: ../Test/LorenzPlotTrajectories.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o"
	cd /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/Test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o -c /u/anmc/garegnan/Desktop/Project/Code/Test/LorenzPlotTrajectories.cpp

Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.i"
	cd /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/Test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/anmc/garegnan/Desktop/Project/Code/Test/LorenzPlotTrajectories.cpp > CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.i

Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.s"
	cd /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/Test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/anmc/garegnan/Desktop/Project/Code/Test/LorenzPlotTrajectories.cpp -o CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.s

Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o.requires:

.PHONY : Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o.requires

Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o.provides: Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o.requires
	$(MAKE) -f Test/CMakeFiles/LorenzTraj.dir/build.make Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o.provides.build
.PHONY : Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o.provides

Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o.provides.build: Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o


Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o: Test/CMakeFiles/LorenzTraj.dir/flags.make
Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o: ../Test/problems.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o"
	cd /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/Test && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/LorenzTraj.dir/problems.cpp.o -c /u/anmc/garegnan/Desktop/Project/Code/Test/problems.cpp

Test/CMakeFiles/LorenzTraj.dir/problems.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LorenzTraj.dir/problems.cpp.i"
	cd /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/Test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /u/anmc/garegnan/Desktop/Project/Code/Test/problems.cpp > CMakeFiles/LorenzTraj.dir/problems.cpp.i

Test/CMakeFiles/LorenzTraj.dir/problems.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LorenzTraj.dir/problems.cpp.s"
	cd /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/Test && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /u/anmc/garegnan/Desktop/Project/Code/Test/problems.cpp -o CMakeFiles/LorenzTraj.dir/problems.cpp.s

Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o.requires:

.PHONY : Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o.requires

Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o.provides: Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o.requires
	$(MAKE) -f Test/CMakeFiles/LorenzTraj.dir/build.make Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o.provides.build
.PHONY : Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o.provides

Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o.provides.build: Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o


# Object files for target LorenzTraj
LorenzTraj_OBJECTS = \
"CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o" \
"CMakeFiles/LorenzTraj.dir/problems.cpp.o"

# External object files for target LorenzTraj
LorenzTraj_EXTERNAL_OBJECTS =

Test/LorenzTraj: Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o
Test/LorenzTraj: Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o
Test/LorenzTraj: Test/CMakeFiles/LorenzTraj.dir/build.make
Test/LorenzTraj: Solver/libSolver.a
Test/LorenzTraj: Test/CMakeFiles/LorenzTraj.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable LorenzTraj"
	cd /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/Test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LorenzTraj.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Test/CMakeFiles/LorenzTraj.dir/build: Test/LorenzTraj

.PHONY : Test/CMakeFiles/LorenzTraj.dir/build

Test/CMakeFiles/LorenzTraj.dir/requires: Test/CMakeFiles/LorenzTraj.dir/LorenzPlotTrajectories.cpp.o.requires
Test/CMakeFiles/LorenzTraj.dir/requires: Test/CMakeFiles/LorenzTraj.dir/problems.cpp.o.requires

.PHONY : Test/CMakeFiles/LorenzTraj.dir/requires

Test/CMakeFiles/LorenzTraj.dir/clean:
	cd /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/Test && $(CMAKE_COMMAND) -P CMakeFiles/LorenzTraj.dir/cmake_clean.cmake
.PHONY : Test/CMakeFiles/LorenzTraj.dir/clean

Test/CMakeFiles/LorenzTraj.dir/depend:
	cd /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /u/anmc/garegnan/Desktop/Project/Code /u/anmc/garegnan/Desktop/Project/Code/Test /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/Test /u/anmc/garegnan/Desktop/Project/Code/cmake-build-debug/Test/CMakeFiles/LorenzTraj.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Test/CMakeFiles/LorenzTraj.dir/depend

