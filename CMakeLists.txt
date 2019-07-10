
cmake_minimum_required (VERSION 3.9)

add_subdirectory(src/libmrc)

project(PSC)

enable_testing()

function(psc_option name description default)
  set(PSC_USE_${name} ${default} CACHE STRING "${description}")
  set_property(CACHE PSC_USE_${name} PROPERTY
    STRINGS "ON;TRUE;AUTO;OFF;FALSE"
  )
endfunction()

option(USE_CUDA "Build CUDA components" OFF)
option(USE_VPIC "Interface with VPIC" OFF)
psc_option(ADIOS2 "Build with adios2 support" AUTO)

# CUDA
if(USE_CUDA)
  enable_language(CUDA)
endif()

find_package(MPI REQUIRED)
add_definitions(-DOMPI_SKIP_MPICXX)

# ADIOS2
if(PSC_USE_ADIOS2 STREQUAL AUTO)
  find_package(ADIOS2)
elseif(PSC_USE_ADIOS2)
  find_package(ADIOS2 REQUIRED)
endif()
if(ADIOS2_FOUND)
  set(PSC_HAVE_ADIOS2 1)
endif()

find_package(GTest)

function(GenerateHeaderConfig)
  set(PSC_CONFIG_DEFINES)
  foreach(OPT IN LISTS ARGN)
    string(APPEND PSC_CONFIG_DEFINES "
/* CMake Option: PSC_USE_${OPT}=${PSC_USE_${OPT}} */
#cmakedefine PSC_HAVE_${OPT}
")
  endforeach()

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/src/include/PscConfig.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/src/include/PscConfig.h.in
  )
  configure_file(
    ${CMAKE_CURRENT_BINARY_DIR}/src/include/PscConfig.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/src/include/PscConfig.h
  )
endfunction()

GenerateHeaderConfig(ADIOS2)

# FIXME, merge with other config
# FIXME, unify options / autodetect
# FIXME, mv helpers into separate file
configure_file (
  "${PROJECT_SOURCE_DIR}/src/include/psc_config.h.in"
  "${PROJECT_BINARY_DIR}/src/include/psc_config.h"
  )
include_directories(${CMAKE_CURRENT_BINARY_DIR}/src/include)
# FIXME, this seems too ugly to find mrc_config.h
include_directories(${CMAKE_CURRENT_BINARY_DIR}/src/libmrc/include)

add_subdirectory(src)



