#
# This is an example CMakeLists.txt file to compile a program against mathplot
#
# Mathplot provides all its functionality as modules. What you do is
# compile these modules - the mathplot core - into a library. You then
# compile your own program, importing the mathplot features you need
# and linking to this library.
#
cmake_minimum_required(VERSION 3.28.5 FATAL_ERROR)

# Give your project a name. Specify language CXX (for maths, mathplot and your code) and also C (for glad)
project(mathplot_template LANGUAGES CXX C)

# mathplot uses c++-20 language features
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#
# Find the libraries which will be needed by mathplot
#

set(OpenGL_GL_PREFERENCE "GLVND") # Following `cmake --help-policy CMP0072`
find_package(OpenGL REQUIRED)
find_package(glfw3 3.3 REQUIRED)
find_package(Freetype REQUIRED)

#
# Define locations of our three submodules, located in extern/
#

set(MPLOT_SOURCE_DIR "${PROJECT_SOURCE_DIR}/mathplot" CACHE PATH "The path to sebsjames/mathplot")
set(MATHS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/maths" CACHE PATH "The path to sebsjames/maths")
set(JSON_SOURCE_DIR "${PROJECT_SOURCE_DIR}/json" CACHE PATH "The path to nlohmann/json")
set(MATHS_CONFIG_SOURCE_DIR "${PROJECT_SOURCE_DIR}/maths_config" CACHE PATH "The path to sebsjames/maths_config")

#
# Set compiler flags
#

# Assume a gcc compiler (or a gcc mimic like Clang).
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -O3 -Wno-unused-result -Wno-unknown-pragmas -Wfatal-errors")

# Tell the program where the mathplot fonts are, to compile them into the binary
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMPLOT_FONTS_DIR=\"\\\"${MPLOT_SOURCE_DIR}/fonts\\\"\"")

# Additional GL compiler flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL3_PROTOTYPES -DGL_GLEXT_PROTOTYPES")

#
# Include paths for module compilation
#

# mathplot needs these includes
include_directories(${OPENGL_INCLUDE_DIR} ${FREETYPE_INCLUDE_DIRS} ${MPLOT_SOURCE_DIR})

# This is required for compilation of the nlohmann/json module
include_directories ("${JSON_SOURCE_DIR}/include")

#
# Define lists of module files to build. We pull in some definitions
# of variables that list the modules in maths and mathplot
#

# Include maths module defs by running the sebsjames/maths macro 'setup_module_variables_for_maths'
include(${MATHS_SOURCE_DIR}/cmake/module_definitions.cmake)
setup_module_variables_for_maths (${MATHS_SOURCE_DIR})
# Include mathplot module defs by running two macros from sebsjames/mathplot
include(${MPLOT_SOURCE_DIR}/cmake/module_definitions.cmake)
setup_module_variables_for_mathplot (${MPLOT_SOURCE_DIR} ${MATHS_SOURCE_DIR} ${JSON_SOURCE_DIR})
# maths_config (for sm::config)
include(${MATHS_CONFIG_SOURCE_DIR}/cmake/module_definitions.cmake)
setup_module_variables_for_maths_config (${MATHS_CONFIG_SOURCE_DIR} ${JSON_SOURCE_DIR} ${MATHS_SOURCE_DIR})

#
# This program uses sm::hdfdata
#
find_package(HDF5 REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HDF5_DEFINITIONS}")
include_directories(${HDF5_INCLUDE_DIR})

#
# Compile two 'boilerplate' libraries - glad and the mathplot fontdata
#

# For mathplot, we must build the glad library here
add_library(glad STATIC ${MPLOT_SOURCE_DIR}/mplot/glad/gl.c)

# We have to build the mathplot fonts data library - this packs up
# some TTF fonts into an object file; they're unpacked at runtime.
add_library(mplot_fontdata STATIC ${MPLOT_SOURCE_DIR}/mplot/VisualFaceAsm.cpp)

#
# Compile a library with the mathplot and locally defined features we want
#

set(MY_RD_MODULES
  ${SM_CONFIG_MODULES}
  ${SM_HDFDATA_MODULES}
  ${SM_HEXGRID_HDF_MODULES}
  ${MPLOT_CORE_MODULES}
  ${MPLOT_MATHS_READCURVES_MODULES}
  ${MPLOT_READCURVES_MODULES}
  ${MPLOT_SOURCE_DIR}/mplot/HexGridVisual.cppm
  ${PROJECT_SOURCE_DIR}/rd_base.cppm
  ${PROJECT_SOURCE_DIR}/rd_schnakenberg.cppm
)
list(REMOVE_DUPLICATES MY_RD_MODULES)
add_library(rd_modules SHARED)
target_sources(rd_modules
  PUBLIC FILE_SET CXX_MODULES BASE_DIRS ${PROJECT_SOURCE_DIR}
  FILES ${MY_RD_MODULES}
)
target_link_libraries(rd_modules OpenGL::GL Freetype::Freetype glfw glad mplot_fontdata ${HDF5_LIBRARIES})

#
# Define the cmake target schnakenberg
#

# The sample is compiled from a c++ file, all the mathplot modules, and the reaction diffusion system module.
add_executable(schnakenberg schnakenberg.cpp)

# Link our internal library rd_modules
target_link_libraries(schnakenberg rd_modules)

#
# All done!
#
