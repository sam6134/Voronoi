# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.19

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.19.1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.19.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/samarthsingh/Voronoi

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/samarthsingh/Voronoi

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/local/Cellar/cmake/3.19.1/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/local/Cellar/cmake/3.19.1/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/samarthsingh/Voronoi/CMakeFiles /Users/samarthsingh/Voronoi//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/samarthsingh/Voronoi/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named cgal_check_build_system

# Build rule for target.
cgal_check_build_system: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 cgal_check_build_system
.PHONY : cgal_check_build_system

# fast build rule for target.
cgal_check_build_system/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/cgal_check_build_system.dir/build.make CMakeFiles/cgal_check_build_system.dir/build
.PHONY : cgal_check_build_system/fast

#=============================================================================
# Target rules for targets named ALL_CGAL_TARGETS

# Build rule for target.
ALL_CGAL_TARGETS: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 ALL_CGAL_TARGETS
.PHONY : ALL_CGAL_TARGETS

# fast build rule for target.
ALL_CGAL_TARGETS/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/ALL_CGAL_TARGETS.dir/build.make CMakeFiles/ALL_CGAL_TARGETS.dir/build
.PHONY : ALL_CGAL_TARGETS/fast

#=============================================================================
# Target rules for targets named voronoi

# Build rule for target.
voronoi: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 voronoi
.PHONY : voronoi

# fast build rule for target.
voronoi/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/voronoi.dir/build.make CMakeFiles/voronoi.dir/build
.PHONY : voronoi/fast

voronoi.o: voronoi.cpp.o

.PHONY : voronoi.o

# target to build an object file
voronoi.cpp.o:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/voronoi.dir/build.make CMakeFiles/voronoi.dir/voronoi.cpp.o
.PHONY : voronoi.cpp.o

voronoi.i: voronoi.cpp.i

.PHONY : voronoi.i

# target to preprocess a source file
voronoi.cpp.i:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/voronoi.dir/build.make CMakeFiles/voronoi.dir/voronoi.cpp.i
.PHONY : voronoi.cpp.i

voronoi.s: voronoi.cpp.s

.PHONY : voronoi.s

# target to generate assembly for a file
voronoi.cpp.s:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/voronoi.dir/build.make CMakeFiles/voronoi.dir/voronoi.cpp.s
.PHONY : voronoi.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... rebuild_cache"
	@echo "... ALL_CGAL_TARGETS"
	@echo "... cgal_check_build_system"
	@echo "... voronoi"
	@echo "... voronoi.o"
	@echo "... voronoi.i"
	@echo "... voronoi.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

