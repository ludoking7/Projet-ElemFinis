# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2023.1.3\bin\cmake\win\x64\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2023.1.3\bin\cmake\win\x64\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\cmake-build-debug"

# Include any dependencies generated for this target.
include CMakeFiles/myFem.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/myFem.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/myFem.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/myFem.dir/flags.make

CMakeFiles/myFem.dir/src/fem.c.obj: CMakeFiles/myFem.dir/flags.make
CMakeFiles/myFem.dir/src/fem.c.obj: C:/Users/ludov/OneDrive/Documents/UCLouvain/Q6/ElemFinis/Project\ ElemFinis/Project/src/fem.c
CMakeFiles/myFem.dir/src/fem.c.obj: CMakeFiles/myFem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/myFem.dir/src/fem.c.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/myFem.dir/src/fem.c.obj -MF CMakeFiles\myFem.dir\src\fem.c.obj.d -o CMakeFiles\myFem.dir\src\fem.c.obj -c "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\src\fem.c"

CMakeFiles/myFem.dir/src/fem.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myFem.dir/src/fem.c.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\src\fem.c" > CMakeFiles\myFem.dir\src\fem.c.i

CMakeFiles/myFem.dir/src/fem.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myFem.dir/src/fem.c.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\src\fem.c" -o CMakeFiles\myFem.dir\src\fem.c.s

CMakeFiles/myFem.dir/src/homework.c.obj: CMakeFiles/myFem.dir/flags.make
CMakeFiles/myFem.dir/src/homework.c.obj: C:/Users/ludov/OneDrive/Documents/UCLouvain/Q6/ElemFinis/Project\ ElemFinis/Project/src/homework.c
CMakeFiles/myFem.dir/src/homework.c.obj: CMakeFiles/myFem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/myFem.dir/src/homework.c.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/myFem.dir/src/homework.c.obj -MF CMakeFiles\myFem.dir\src\homework.c.obj.d -o CMakeFiles\myFem.dir\src\homework.c.obj -c "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\src\homework.c"

CMakeFiles/myFem.dir/src/homework.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myFem.dir/src/homework.c.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\src\homework.c" > CMakeFiles\myFem.dir\src\homework.c.i

CMakeFiles/myFem.dir/src/homework.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myFem.dir/src/homework.c.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\src\homework.c" -o CMakeFiles\myFem.dir\src\homework.c.s

CMakeFiles/myFem.dir/src/main.c.obj: CMakeFiles/myFem.dir/flags.make
CMakeFiles/myFem.dir/src/main.c.obj: C:/Users/ludov/OneDrive/Documents/UCLouvain/Q6/ElemFinis/Project\ ElemFinis/Project/src/main.c
CMakeFiles/myFem.dir/src/main.c.obj: CMakeFiles/myFem.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building C object CMakeFiles/myFem.dir/src/main.c.obj"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/myFem.dir/src/main.c.obj -MF CMakeFiles\myFem.dir\src\main.c.obj.d -o CMakeFiles\myFem.dir\src\main.c.obj -c "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\src\main.c"

CMakeFiles/myFem.dir/src/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/myFem.dir/src/main.c.i"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\src\main.c" > CMakeFiles\myFem.dir\src\main.c.i

CMakeFiles/myFem.dir/src/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/myFem.dir/src/main.c.s"
	C:\PROGRA~1\JETBRA~1\CLION2~1.3\bin\mingw\bin\gcc.exe $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\src\main.c" -o CMakeFiles\myFem.dir\src\main.c.s

# Object files for target myFem
myFem_OBJECTS = \
"CMakeFiles/myFem.dir/src/fem.c.obj" \
"CMakeFiles/myFem.dir/src/homework.c.obj" \
"CMakeFiles/myFem.dir/src/main.c.obj"

# External object files for target myFem
myFem_EXTERNAL_OBJECTS =

myFem.exe: CMakeFiles/myFem.dir/src/fem.c.obj
myFem.exe: CMakeFiles/myFem.dir/src/homework.c.obj
myFem.exe: CMakeFiles/myFem.dir/src/main.c.obj
myFem.exe: CMakeFiles/myFem.dir/build.make
myFem.exe: CMakeFiles/myFem.dir/linkLibs.rsp
myFem.exe: CMakeFiles/myFem.dir/objects1
myFem.exe: CMakeFiles/myFem.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Linking C executable myFem.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\myFem.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/myFem.dir/build: myFem.exe
.PHONY : CMakeFiles/myFem.dir/build

CMakeFiles/myFem.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\myFem.dir\cmake_clean.cmake
.PHONY : CMakeFiles/myFem.dir/clean

CMakeFiles/myFem.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project" "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project" "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\cmake-build-debug" "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\cmake-build-debug" "C:\Users\ludov\OneDrive\Documents\UCLouvain\Q6\ElemFinis\Project ElemFinis\Project\cmake-build-debug\CMakeFiles\myFem.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/myFem.dir/depend

