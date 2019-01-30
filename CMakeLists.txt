
cmake_minimum_required (VERSION 3.9)

add_subdirectory(src/libmrc)

project(PSC)

enable_testing()

option(USE_CUDA "Build CUDA components" OFF)

if(USE_CUDA)
  enable_language(CUDA)
endif()

option(USE_VPIC "Interface with VPIC" OFF)
if(USE_VPIC)
  add_library(VPIC::VPIC INTERFACE IMPORTED)
  set_target_properties(VPIC::VPIC PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "/Users/kai/build/vpic/src"
    INTERFACE_LINK_LIBRARIES "/Users/kai/build/vpic/build/libvpic.a"
  )
endif()

find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)

find_package(GTest)

configure_file (
  "${PROJECT_SOURCE_DIR}/src/include/psc_config.h.in"
  "${PROJECT_BINARY_DIR}/src/include/psc_config.h"
  )
include_directories(${CMAKE_CURRENT_BINARY_DIR}/src/include)
# FIXME, this seems too ugly to find mrc_config.h
include_directories(${CMAKE_CURRENT_BINARY_DIR}/src/libmrc/include)

add_subdirectory(src)



