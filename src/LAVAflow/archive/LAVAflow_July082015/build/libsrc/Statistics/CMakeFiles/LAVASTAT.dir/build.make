# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_COMMAND = /opt/software/tools/cmake/2.8.12/bin/cmake

# The command to remove a file.
RM = /opt/software/tools/cmake/2.8.12/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /opt/software/tools/cmake/2.8.12/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/tah09e/code/workspace/lavaflow

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/tah09e/code/workspace/lavaflow/build

# Include any dependencies generated for this target.
include libsrc/Statistics/CMakeFiles/LAVASTAT.dir/depend.make

# Include the progress variables for this target.
include libsrc/Statistics/CMakeFiles/LAVASTAT.dir/progress.make

# Include the compile flags for this target's objects.
include libsrc/Statistics/CMakeFiles/LAVASTAT.dir/flags.make

libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o: libsrc/Statistics/CMakeFiles/LAVASTAT.dir/flags.make
libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o: ../libsrc/Statistics/Histogram/Histogram.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /home/tah09e/code/workspace/lavaflow/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o"
	cd /home/tah09e/code/workspace/lavaflow/build/libsrc/Statistics && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o -c /home/tah09e/code/workspace/lavaflow/libsrc/Statistics/Histogram/Histogram.cxx

libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.i"
	cd /home/tah09e/code/workspace/lavaflow/build/libsrc/Statistics && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/tah09e/code/workspace/lavaflow/libsrc/Statistics/Histogram/Histogram.cxx > CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.i

libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.s"
	cd /home/tah09e/code/workspace/lavaflow/build/libsrc/Statistics && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/tah09e/code/workspace/lavaflow/libsrc/Statistics/Histogram/Histogram.cxx -o CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.s

libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o.requires:
.PHONY : libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o.requires

libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o.provides: libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o.requires
	$(MAKE) -f libsrc/Statistics/CMakeFiles/LAVASTAT.dir/build.make libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o.provides.build
.PHONY : libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o.provides

libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o.provides.build: libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o

# Object files for target LAVASTAT
LAVASTAT_OBJECTS = \
"CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o"

# External object files for target LAVASTAT
LAVASTAT_EXTERNAL_OBJECTS =

libsrc/Statistics/libLAVASTAT.a: libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o
libsrc/Statistics/libLAVASTAT.a: libsrc/Statistics/CMakeFiles/LAVASTAT.dir/build.make
libsrc/Statistics/libLAVASTAT.a: libsrc/Statistics/CMakeFiles/LAVASTAT.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libLAVASTAT.a"
	cd /home/tah09e/code/workspace/lavaflow/build/libsrc/Statistics && $(CMAKE_COMMAND) -P CMakeFiles/LAVASTAT.dir/cmake_clean_target.cmake
	cd /home/tah09e/code/workspace/lavaflow/build/libsrc/Statistics && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/LAVASTAT.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libsrc/Statistics/CMakeFiles/LAVASTAT.dir/build: libsrc/Statistics/libLAVASTAT.a
.PHONY : libsrc/Statistics/CMakeFiles/LAVASTAT.dir/build

libsrc/Statistics/CMakeFiles/LAVASTAT.dir/requires: libsrc/Statistics/CMakeFiles/LAVASTAT.dir/Histogram/Histogram.cxx.o.requires
.PHONY : libsrc/Statistics/CMakeFiles/LAVASTAT.dir/requires

libsrc/Statistics/CMakeFiles/LAVASTAT.dir/clean:
	cd /home/tah09e/code/workspace/lavaflow/build/libsrc/Statistics && $(CMAKE_COMMAND) -P CMakeFiles/LAVASTAT.dir/cmake_clean.cmake
.PHONY : libsrc/Statistics/CMakeFiles/LAVASTAT.dir/clean

libsrc/Statistics/CMakeFiles/LAVASTAT.dir/depend:
	cd /home/tah09e/code/workspace/lavaflow/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/tah09e/code/workspace/lavaflow /home/tah09e/code/workspace/lavaflow/libsrc/Statistics /home/tah09e/code/workspace/lavaflow/build /home/tah09e/code/workspace/lavaflow/build/libsrc/Statistics /home/tah09e/code/workspace/lavaflow/build/libsrc/Statistics/CMakeFiles/LAVASTAT.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libsrc/Statistics/CMakeFiles/LAVASTAT.dir/depend
