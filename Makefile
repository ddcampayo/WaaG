# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_SOURCE_DIR = /home/daniel/WaaG

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/daniel/WaaG

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# Special rule for the target test
test:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running tests..."
	/usr/bin/ctest --force-new-ctest-process $(ARGS)
.PHONY : test

# Special rule for the target test
test/fast: test

.PHONY : test/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/daniel/WaaG/CMakeFiles /home/daniel/WaaG/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/daniel/WaaG/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named pPart

# Build rule for target.
pPart: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 pPart
.PHONY : pPart

# fast build rule for target.
pPart/fast:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/build
.PHONY : pPart/fast

create.o: create.cpp.o

.PHONY : create.o

# target to build an object file
create.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/create.cpp.o
.PHONY : create.cpp.o

create.i: create.cpp.i

.PHONY : create.i

# target to preprocess a source file
create.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/create.cpp.i
.PHONY : create.cpp.i

create.s: create.cpp.s

.PHONY : create.s

# target to generate assembly for a file
create.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/create.cpp.s
.PHONY : create.cpp.s

draw.o: draw.cpp.o

.PHONY : draw.o

# target to build an object file
draw.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/draw.cpp.o
.PHONY : draw.cpp.o

draw.i: draw.cpp.i

.PHONY : draw.i

# target to preprocess a source file
draw.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/draw.cpp.i
.PHONY : draw.cpp.i

draw.s: draw.cpp.s

.PHONY : draw.s

# target to generate assembly for a file
draw.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/draw.cpp.s
.PHONY : draw.cpp.s

fields.o: fields.cpp.o

.PHONY : fields.o

# target to build an object file
fields.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/fields.cpp.o
.PHONY : fields.cpp.o

fields.i: fields.cpp.i

.PHONY : fields.i

# target to preprocess a source file
fields.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/fields.cpp.i
.PHONY : fields.cpp.i

fields.s: fields.cpp.s

.PHONY : fields.s

# target to generate assembly for a file
fields.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/fields.cpp.s
.PHONY : fields.cpp.s

linear_DD_scalar_prod.o: linear_DD_scalar_prod.cpp.o

.PHONY : linear_DD_scalar_prod.o

# target to build an object file
linear_DD_scalar_prod.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_DD_scalar_prod.cpp.o
.PHONY : linear_DD_scalar_prod.cpp.o

linear_DD_scalar_prod.i: linear_DD_scalar_prod.cpp.i

.PHONY : linear_DD_scalar_prod.i

# target to preprocess a source file
linear_DD_scalar_prod.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_DD_scalar_prod.cpp.i
.PHONY : linear_DD_scalar_prod.cpp.i

linear_DD_scalar_prod.s: linear_DD_scalar_prod.cpp.s

.PHONY : linear_DD_scalar_prod.s

# target to generate assembly for a file
linear_DD_scalar_prod.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_DD_scalar_prod.cpp.s
.PHONY : linear_DD_scalar_prod.cpp.s

linear_aux.o: linear_aux.cpp.o

.PHONY : linear_aux.o

# target to build an object file
linear_aux.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_aux.cpp.o
.PHONY : linear_aux.cpp.o

linear_aux.i: linear_aux.cpp.i

.PHONY : linear_aux.i

# target to preprocess a source file
linear_aux.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_aux.cpp.i
.PHONY : linear_aux.cpp.i

linear_aux.s: linear_aux.cpp.s

.PHONY : linear_aux.s

# target to generate assembly for a file
linear_aux.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_aux.cpp.s
.PHONY : linear_aux.cpp.s

linear_fill_Delta_DD.o: linear_fill_Delta_DD.cpp.o

.PHONY : linear_fill_Delta_DD.o

# target to build an object file
linear_fill_Delta_DD.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_fill_Delta_DD.cpp.o
.PHONY : linear_fill_Delta_DD.cpp.o

linear_fill_Delta_DD.i: linear_fill_Delta_DD.cpp.i

.PHONY : linear_fill_Delta_DD.i

# target to preprocess a source file
linear_fill_Delta_DD.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_fill_Delta_DD.cpp.i
.PHONY : linear_fill_Delta_DD.cpp.i

linear_fill_Delta_DD.s: linear_fill_Delta_DD.cpp.s

.PHONY : linear_fill_Delta_DD.s

# target to generate assembly for a file
linear_fill_Delta_DD.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_fill_Delta_DD.cpp.s
.PHONY : linear_fill_Delta_DD.cpp.s

linear_p_equation.o: linear_p_equation.cpp.o

.PHONY : linear_p_equation.o

# target to build an object file
linear_p_equation.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_p_equation.cpp.o
.PHONY : linear_p_equation.cpp.o

linear_p_equation.i: linear_p_equation.cpp.i

.PHONY : linear_p_equation.i

# target to preprocess a source file
linear_p_equation.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_p_equation.cpp.i
.PHONY : linear_p_equation.cpp.i

linear_p_equation.s: linear_p_equation.cpp.s

.PHONY : linear_p_equation.s

# target to generate assembly for a file
linear_p_equation.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_p_equation.cpp.s
.PHONY : linear_p_equation.cpp.s

linear_s_equation.o: linear_s_equation.cpp.o

.PHONY : linear_s_equation.o

# target to build an object file
linear_s_equation.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_s_equation.cpp.o
.PHONY : linear_s_equation.cpp.o

linear_s_equation.i: linear_s_equation.cpp.i

.PHONY : linear_s_equation.i

# target to preprocess a source file
linear_s_equation.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_s_equation.cpp.i
.PHONY : linear_s_equation.cpp.i

linear_s_equation.s: linear_s_equation.cpp.s

.PHONY : linear_s_equation.s

# target to generate assembly for a file
linear_s_equation.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_s_equation.cpp.s
.PHONY : linear_s_equation.cpp.s

linear_solve_for_moments.o: linear_solve_for_moments.cpp.o

.PHONY : linear_solve_for_moments.o

# target to build an object file
linear_solve_for_moments.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_solve_for_moments.cpp.o
.PHONY : linear_solve_for_moments.cpp.o

linear_solve_for_moments.i: linear_solve_for_moments.cpp.i

.PHONY : linear_solve_for_moments.i

# target to preprocess a source file
linear_solve_for_moments.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_solve_for_moments.cpp.i
.PHONY : linear_solve_for_moments.cpp.i

linear_solve_for_moments.s: linear_solve_for_moments.cpp.s

.PHONY : linear_solve_for_moments.s

# target to generate assembly for a file
linear_solve_for_moments.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_solve_for_moments.cpp.s
.PHONY : linear_solve_for_moments.cpp.s

linear_solve_for_weights.o: linear_solve_for_weights.cpp.o

.PHONY : linear_solve_for_weights.o

# target to build an object file
linear_solve_for_weights.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_solve_for_weights.cpp.o
.PHONY : linear_solve_for_weights.cpp.o

linear_solve_for_weights.i: linear_solve_for_weights.cpp.i

.PHONY : linear_solve_for_weights.i

# target to preprocess a source file
linear_solve_for_weights.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_solve_for_weights.cpp.i
.PHONY : linear_solve_for_weights.cpp.i

linear_solve_for_weights.s: linear_solve_for_weights.cpp.s

.PHONY : linear_solve_for_weights.s

# target to generate assembly for a file
linear_solve_for_weights.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_solve_for_weights.cpp.s
.PHONY : linear_solve_for_weights.cpp.s

linear_vect_to_field.o: linear_vect_to_field.cpp.o

.PHONY : linear_vect_to_field.o

# target to build an object file
linear_vect_to_field.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_vect_to_field.cpp.o
.PHONY : linear_vect_to_field.cpp.o

linear_vect_to_field.i: linear_vect_to_field.cpp.i

.PHONY : linear_vect_to_field.i

# target to preprocess a source file
linear_vect_to_field.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_vect_to_field.cpp.i
.PHONY : linear_vect_to_field.cpp.i

linear_vect_to_field.s: linear_vect_to_field.cpp.s

.PHONY : linear_vect_to_field.s

# target to generate assembly for a file
linear_vect_to_field.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_vect_to_field.cpp.s
.PHONY : linear_vect_to_field.cpp.s

linear_w_equation.o: linear_w_equation.cpp.o

.PHONY : linear_w_equation.o

# target to build an object file
linear_w_equation.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_w_equation.cpp.o
.PHONY : linear_w_equation.cpp.o

linear_w_equation.i: linear_w_equation.cpp.i

.PHONY : linear_w_equation.i

# target to preprocess a source file
linear_w_equation.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_w_equation.cpp.i
.PHONY : linear_w_equation.cpp.i

linear_w_equation.s: linear_w_equation.cpp.s

.PHONY : linear_w_equation.s

# target to generate assembly for a file
linear_w_equation.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_w_equation.cpp.s
.PHONY : linear_w_equation.cpp.s

linear_w_equation2.o: linear_w_equation2.cpp.o

.PHONY : linear_w_equation2.o

# target to build an object file
linear_w_equation2.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_w_equation2.cpp.o
.PHONY : linear_w_equation2.cpp.o

linear_w_equation2.i: linear_w_equation2.cpp.i

.PHONY : linear_w_equation2.i

# target to preprocess a source file
linear_w_equation2.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_w_equation2.cpp.i
.PHONY : linear_w_equation2.cpp.i

linear_w_equation2.s: linear_w_equation2.cpp.s

.PHONY : linear_w_equation2.s

# target to generate assembly for a file
linear_w_equation2.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_w_equation2.cpp.s
.PHONY : linear_w_equation2.cpp.s

linear_w_equation3.o: linear_w_equation3.cpp.o

.PHONY : linear_w_equation3.o

# target to build an object file
linear_w_equation3.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_w_equation3.cpp.o
.PHONY : linear_w_equation3.cpp.o

linear_w_equation3.i: linear_w_equation3.cpp.i

.PHONY : linear_w_equation3.i

# target to preprocess a source file
linear_w_equation3.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_w_equation3.cpp.i
.PHONY : linear_w_equation3.cpp.i

linear_w_equation3.s: linear_w_equation3.cpp.s

.PHONY : linear_w_equation3.s

# target to generate assembly for a file
linear_w_equation3.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/linear_w_equation3.cpp.s
.PHONY : linear_w_equation3.cpp.s

lloyds.o: lloyds.cpp.o

.PHONY : lloyds.o

# target to build an object file
lloyds.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/lloyds.cpp.o
.PHONY : lloyds.cpp.o

lloyds.i: lloyds.cpp.i

.PHONY : lloyds.i

# target to preprocess a source file
lloyds.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/lloyds.cpp.i
.PHONY : lloyds.cpp.i

lloyds.s: lloyds.cpp.s

.PHONY : lloyds.s

# target to generate assembly for a file
lloyds.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/lloyds.cpp.s
.PHONY : lloyds.cpp.s

move.o: move.cpp.o

.PHONY : move.o

# target to build an object file
move.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/move.cpp.o
.PHONY : move.cpp.o

move.i: move.cpp.i

.PHONY : move.i

# target to preprocess a source file
move.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/move.cpp.i
.PHONY : move.cpp.i

move.s: move.cpp.s

.PHONY : move.s

# target to generate assembly for a file
move.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/move.cpp.s
.PHONY : move.cpp.s

number.o: number.cpp.o

.PHONY : number.o

# target to build an object file
number.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/number.cpp.o
.PHONY : number.cpp.o

number.i: number.cpp.i

.PHONY : number.i

# target to preprocess a source file
number.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/number.cpp.i
.PHONY : number.cpp.i

number.s: number.cpp.s

.PHONY : number.s

# target to generate assembly for a file
number.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/number.cpp.s
.PHONY : number.cpp.s

pParticles_dG_rev.o: pParticles_dG_rev.cpp.o

.PHONY : pParticles_dG_rev.o

# target to build an object file
pParticles_dG_rev.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/pParticles_dG_rev.cpp.o
.PHONY : pParticles_dG_rev.cpp.o

pParticles_dG_rev.i: pParticles_dG_rev.cpp.i

.PHONY : pParticles_dG_rev.i

# target to preprocess a source file
pParticles_dG_rev.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/pParticles_dG_rev.cpp.i
.PHONY : pParticles_dG_rev.cpp.i

pParticles_dG_rev.s: pParticles_dG_rev.cpp.s

.PHONY : pParticles_dG_rev.s

# target to generate assembly for a file
pParticles_dG_rev.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/pParticles_dG_rev.cpp.s
.PHONY : pParticles_dG_rev.cpp.s

volumes.o: volumes.cpp.o

.PHONY : volumes.o

# target to build an object file
volumes.cpp.o:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/volumes.cpp.o
.PHONY : volumes.cpp.o

volumes.i: volumes.cpp.i

.PHONY : volumes.i

# target to preprocess a source file
volumes.cpp.i:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/volumes.cpp.i
.PHONY : volumes.cpp.i

volumes.s: volumes.cpp.s

.PHONY : volumes.s

# target to generate assembly for a file
volumes.cpp.s:
	$(MAKE) -f CMakeFiles/pPart.dir/build.make CMakeFiles/pPart.dir/volumes.cpp.s
.PHONY : volumes.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... pPart"
	@echo "... test"
	@echo "... create.o"
	@echo "... create.i"
	@echo "... create.s"
	@echo "... draw.o"
	@echo "... draw.i"
	@echo "... draw.s"
	@echo "... fields.o"
	@echo "... fields.i"
	@echo "... fields.s"
	@echo "... linear_DD_scalar_prod.o"
	@echo "... linear_DD_scalar_prod.i"
	@echo "... linear_DD_scalar_prod.s"
	@echo "... linear_aux.o"
	@echo "... linear_aux.i"
	@echo "... linear_aux.s"
	@echo "... linear_fill_Delta_DD.o"
	@echo "... linear_fill_Delta_DD.i"
	@echo "... linear_fill_Delta_DD.s"
	@echo "... linear_p_equation.o"
	@echo "... linear_p_equation.i"
	@echo "... linear_p_equation.s"
	@echo "... linear_s_equation.o"
	@echo "... linear_s_equation.i"
	@echo "... linear_s_equation.s"
	@echo "... linear_solve_for_moments.o"
	@echo "... linear_solve_for_moments.i"
	@echo "... linear_solve_for_moments.s"
	@echo "... linear_solve_for_weights.o"
	@echo "... linear_solve_for_weights.i"
	@echo "... linear_solve_for_weights.s"
	@echo "... linear_vect_to_field.o"
	@echo "... linear_vect_to_field.i"
	@echo "... linear_vect_to_field.s"
	@echo "... linear_w_equation.o"
	@echo "... linear_w_equation.i"
	@echo "... linear_w_equation.s"
	@echo "... linear_w_equation2.o"
	@echo "... linear_w_equation2.i"
	@echo "... linear_w_equation2.s"
	@echo "... linear_w_equation3.o"
	@echo "... linear_w_equation3.i"
	@echo "... linear_w_equation3.s"
	@echo "... lloyds.o"
	@echo "... lloyds.i"
	@echo "... lloyds.s"
	@echo "... move.o"
	@echo "... move.i"
	@echo "... move.s"
	@echo "... number.o"
	@echo "... number.i"
	@echo "... number.s"
	@echo "... pParticles_dG_rev.o"
	@echo "... pParticles_dG_rev.i"
	@echo "... pParticles_dG_rev.s"
	@echo "... volumes.o"
	@echo "... volumes.i"
	@echo "... volumes.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

