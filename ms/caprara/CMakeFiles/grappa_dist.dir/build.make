# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/kuko/SVN/vlaciky/ms/caprara/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/kuko/SVN/vlaciky/ms/caprara

# Include any dependencies generated for this target.
include CMakeFiles/grappa_dist.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/grappa_dist.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/grappa_dist.dir/flags.make

CMakeFiles/grappa_dist.dir/invdist.o: CMakeFiles/grappa_dist.dir/flags.make
CMakeFiles/grappa_dist.dir/invdist.o: src/invdist.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/kuko/SVN/vlaciky/ms/caprara/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/grappa_dist.dir/invdist.o"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/grappa_dist.dir/invdist.o   -c /home/kuko/SVN/vlaciky/ms/caprara/src/invdist.c

CMakeFiles/grappa_dist.dir/invdist.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/grappa_dist.dir/invdist.i"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/kuko/SVN/vlaciky/ms/caprara/src/invdist.c > CMakeFiles/grappa_dist.dir/invdist.i

CMakeFiles/grappa_dist.dir/invdist.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/grappa_dist.dir/invdist.s"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/kuko/SVN/vlaciky/ms/caprara/src/invdist.c -o CMakeFiles/grappa_dist.dir/invdist.s

CMakeFiles/grappa_dist.dir/invdist.o.requires:
.PHONY : CMakeFiles/grappa_dist.dir/invdist.o.requires

CMakeFiles/grappa_dist.dir/invdist.o.provides: CMakeFiles/grappa_dist.dir/invdist.o.requires
	$(MAKE) -f CMakeFiles/grappa_dist.dir/build.make CMakeFiles/grappa_dist.dir/invdist.o.provides.build
.PHONY : CMakeFiles/grappa_dist.dir/invdist.o.provides

CMakeFiles/grappa_dist.dir/invdist.o.provides.build: CMakeFiles/grappa_dist.dir/invdist.o
.PHONY : CMakeFiles/grappa_dist.dir/invdist.o.provides.build

CMakeFiles/grappa_dist.dir/med_util.o: CMakeFiles/grappa_dist.dir/flags.make
CMakeFiles/grappa_dist.dir/med_util.o: src/med_util.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/kuko/SVN/vlaciky/ms/caprara/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/grappa_dist.dir/med_util.o"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/grappa_dist.dir/med_util.o   -c /home/kuko/SVN/vlaciky/ms/caprara/src/med_util.c

CMakeFiles/grappa_dist.dir/med_util.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/grappa_dist.dir/med_util.i"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/kuko/SVN/vlaciky/ms/caprara/src/med_util.c > CMakeFiles/grappa_dist.dir/med_util.i

CMakeFiles/grappa_dist.dir/med_util.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/grappa_dist.dir/med_util.s"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/kuko/SVN/vlaciky/ms/caprara/src/med_util.c -o CMakeFiles/grappa_dist.dir/med_util.s

CMakeFiles/grappa_dist.dir/med_util.o.requires:
.PHONY : CMakeFiles/grappa_dist.dir/med_util.o.requires

CMakeFiles/grappa_dist.dir/med_util.o.provides: CMakeFiles/grappa_dist.dir/med_util.o.requires
	$(MAKE) -f CMakeFiles/grappa_dist.dir/build.make CMakeFiles/grappa_dist.dir/med_util.o.provides.build
.PHONY : CMakeFiles/grappa_dist.dir/med_util.o.provides

CMakeFiles/grappa_dist.dir/med_util.o.provides.build: CMakeFiles/grappa_dist.dir/med_util.o
.PHONY : CMakeFiles/grappa_dist.dir/med_util.o.provides.build

CMakeFiles/grappa_dist.dir/uf.o: CMakeFiles/grappa_dist.dir/flags.make
CMakeFiles/grappa_dist.dir/uf.o: src/uf.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/kuko/SVN/vlaciky/ms/caprara/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/grappa_dist.dir/uf.o"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/grappa_dist.dir/uf.o   -c /home/kuko/SVN/vlaciky/ms/caprara/src/uf.c

CMakeFiles/grappa_dist.dir/uf.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/grappa_dist.dir/uf.i"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/kuko/SVN/vlaciky/ms/caprara/src/uf.c > CMakeFiles/grappa_dist.dir/uf.i

CMakeFiles/grappa_dist.dir/uf.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/grappa_dist.dir/uf.s"
	/usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/kuko/SVN/vlaciky/ms/caprara/src/uf.c -o CMakeFiles/grappa_dist.dir/uf.s

CMakeFiles/grappa_dist.dir/uf.o.requires:
.PHONY : CMakeFiles/grappa_dist.dir/uf.o.requires

CMakeFiles/grappa_dist.dir/uf.o.provides: CMakeFiles/grappa_dist.dir/uf.o.requires
	$(MAKE) -f CMakeFiles/grappa_dist.dir/build.make CMakeFiles/grappa_dist.dir/uf.o.provides.build
.PHONY : CMakeFiles/grappa_dist.dir/uf.o.provides

CMakeFiles/grappa_dist.dir/uf.o.provides.build: CMakeFiles/grappa_dist.dir/uf.o
.PHONY : CMakeFiles/grappa_dist.dir/uf.o.provides.build

# Object files for target grappa_dist
grappa_dist_OBJECTS = \
"CMakeFiles/grappa_dist.dir/invdist.o" \
"CMakeFiles/grappa_dist.dir/med_util.o" \
"CMakeFiles/grappa_dist.dir/uf.o"

# External object files for target grappa_dist
grappa_dist_EXTERNAL_OBJECTS =

libgrappa_dist.so: CMakeFiles/grappa_dist.dir/invdist.o
libgrappa_dist.so: CMakeFiles/grappa_dist.dir/med_util.o
libgrappa_dist.so: CMakeFiles/grappa_dist.dir/uf.o
libgrappa_dist.so: CMakeFiles/grappa_dist.dir/build.make
libgrappa_dist.so: CMakeFiles/grappa_dist.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C shared library libgrappa_dist.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/grappa_dist.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/grappa_dist.dir/build: libgrappa_dist.so
.PHONY : CMakeFiles/grappa_dist.dir/build

CMakeFiles/grappa_dist.dir/requires: CMakeFiles/grappa_dist.dir/invdist.o.requires
CMakeFiles/grappa_dist.dir/requires: CMakeFiles/grappa_dist.dir/med_util.o.requires
CMakeFiles/grappa_dist.dir/requires: CMakeFiles/grappa_dist.dir/uf.o.requires
.PHONY : CMakeFiles/grappa_dist.dir/requires

CMakeFiles/grappa_dist.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/grappa_dist.dir/cmake_clean.cmake
.PHONY : CMakeFiles/grappa_dist.dir/clean

CMakeFiles/grappa_dist.dir/depend:
	cd /home/kuko/SVN/vlaciky/ms/caprara && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kuko/SVN/vlaciky/ms/caprara/src /home/kuko/SVN/vlaciky/ms/caprara/src /home/kuko/SVN/vlaciky/ms/caprara /home/kuko/SVN/vlaciky/ms/caprara /home/kuko/SVN/vlaciky/ms/caprara/CMakeFiles/grappa_dist.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/grappa_dist.dir/depend

