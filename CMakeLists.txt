
cmake_minimum_required (VERSION 3.9)

foreach(policy
    CMP0074 # CMake 3.12
    )
  if(POLICY ${policy})
    cmake_policy(SET ${policy} NEW)
  endif()
endforeach()

add_subdirectory(src/libmrc)

project(PSC)

include(CTest)

function(psc_option name description default)
  set(PSC_USE_${name} ${default} CACHE STRING "${description}")
  set_property(CACHE PSC_USE_${name} PROPERTY
    STRINGS "ON;TRUE;AUTO;OFF;FALSE"
  )
endfunction()

option(USE_CUDA "Build CUDA components" OFF)
option(USE_VPIC "Interface with VPIC" OFF)
option(USE_GTEST_DISCOVER_TESTS "Run tests to discover contained googletest cases" OFF)
psc_option(ADIOS2 "Build with adios2 support" AUTO)

# CUDA
if(USE_CUDA)
  list(APPEND CMAKE_CUDA_FLAGS "-std c++14")

  enable_language(CUDA)

  list(REMOVE_ITEM CMAKE_CUDA_FLAGS "-std c++14")
endif()


set(MPI_CXX_SKIP_MPICXX ON)
find_package(MPI REQUIRED)
foreach(lang IN ITEMS CXX C)
  get_property(tmp TARGET MPI::MPI_${lang} PROPERTY INTERFACE_COMPILE_OPTIONS)
  list(TRANSFORM tmp REPLACE "-fexceptions" "-abc $<$<COMPILE_LANGUAGE:CUDA>:SHELL:-Xcompiler >-fexceptions")
  set_property(TARGET MPI::MPI_${lang} PROPERTY INTERFACE_COMPILE_OPTIONS "${tmp}")
  unset(tmp)
endforeach()

# ADIOS2
if(PSC_USE_ADIOS2 STREQUAL AUTO)
  find_package(ADIOS2)
elseif(PSC_USE_ADIOS2)
  find_package(ADIOS2 REQUIRED)
endif()
if(ADIOS2_FOUND)
  set(PSC_HAVE_ADIOS2 1)
endif()

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

# FIXME, unify USE_CUDA, USE_VPIC options / autodetect
# FIXME, mv helpers into separate file
GenerateHeaderConfig(ADIOS2)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/src/include)
# FIXME, this seems too ugly to find mrc_config.h
include_directories(${CMAKE_CURRENT_BINARY_DIR}/src/libmrc/include)

add_subdirectory(external)
if (BUILD_TESTING)
  include(GoogleTest)
endif()

add_subdirectory(src)



