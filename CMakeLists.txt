
cmake_minimum_required (VERSION 3.17)

foreach(policy
    CMP0074 # CMake 3.12
    CMP0104 # CMAKE_CUDA_ARCHITECTURES
    )
  if(POLICY ${policy})
    cmake_policy(SET ${policy} NEW)
  endif()
endforeach()

project(PSC)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(cmake/CPM.cmake)

include(CTest)

function(psc_option name description default)
  set(PSC_USE_${name} ${default} CACHE STRING "${description}")
  set_property(CACHE PSC_USE_${name} PROPERTY
    STRINGS "ON;TRUE;AUTO;OFF;FALSE"
  )
endfunction()

set(PSC_GPU "host" CACHE STRING "Build with GPU supoort (host/cuda/hip)")
set_property(CACHE PSC_GPU PROPERTY STRINGS "host;cuda;hip")
option(USE_VPIC "Interface with VPIC" OFF)
option(USE_GTEST_DISCOVER_TESTS "Run tests to discover contained googletest cases" OFF)
psc_option(ADIOS2 "Build with adios2 support" AUTO)
option(PSC_USE_NVTX "Build with NVTX support" OFF)
option(PSC_USE_RMM "Build with RMM memory manager support" OFF)
option(PSC_BOUNDS_CHECK "Turn on bounds-checking" OFF)
option(PSC_USE_PERFETTO "Turn on perfetto-based tracing" OFF)

option(USE_CUDA "Build CUDA components" OFF)

if (USE_CUDA)
  set(PSC_GPU "cuda" CACHE STRING  "Build with GPU supoort (none/cuda/hip)")
endif()

if(PSC_GPU STREQUAL "cuda")
  # This is needed on Summit when not using modules (in a spack build) to
  # make sure that nvcc uses the same host compiler that was otherwise
  # specified
  set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")
  enable_language(CUDA)

  find_package(Thrust 1.10.0 REQUIRED)
  thrust_create_target(Thrust)
  set(USE_CUDA ON)
elseif(PSC_GPU STREQUAL "hip")
  set(USE_CUDA ON)
else()
  set(USE_CUDA OFF)
endif()
set(GTENSOR_DEVICE "${PSC_GPU}" CACHE STRING "")

CPMAddPackage(
  NAME gtensor
  GITHUB_REPOSITORY wdmapp/gtensor
  GIT_TAG main
  VERSION 0.0.0
)

if (BUILD_TESTING)
  CPMAddPackage(
    NAME GTest
    GITHUB_REPOSITORY google/googletest
    VERSION 1.14.0
    OPTIONS
        "INSTALL_GTEST OFF"
        "gtest_force_shared_crt ON"
  )
  if (GTest_ADDED)
     add_library(GTest::gtest ALIAS gtest)
     add_library(GTest::gtest_main ALIAS gtest_main)
  endif()
  include(GoogleTest)
endif()

set(MPI_CXX_SKIP_MPICXX ON)
find_package(MPI REQUIRED)

# ADIOS2
if(PSC_USE_ADIOS2 STREQUAL AUTO)
  find_package(ADIOS2 CONFIG)
elseif(PSC_USE_ADIOS2)
  find_package(ADIOS2 CONFIG REQUIRED)
endif()
if(ADIOS2_FOUND)
  set(PSC_HAVE_ADIOS2 1)
endif()

# NVTX
if (PSC_USE_NVTX)
  find_package(CUDAToolkit REQUIRED)
  set(PSC_HAVE_NVTX 1)
endif()

# RMM
if (PSC_USE_RMM)
  find_package(rmm 0.18.0 CONFIG REQUIRED)
  set(PSC_HAVE_RMM 1)
endif()

if (PSC_USE_PERFETTO)
  CPMAddPackage("gh:google/perfetto@24.2")
  find_package(Threads REQUIRED)
  add_library(perfetto STATIC ${perfetto_SOURCE_DIR}/sdk/perfetto.cc)
  add_library(perfetto::perfetto ALIAS perfetto)
  target_include_directories(perfetto PUBLIC ${perfetto_SOURCE_DIR}/sdk)
  target_link_libraries(perfetto PRIVATE Threads::Threads)
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
GenerateHeaderConfig(ADIOS2 NVTX RMM)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/src/include)
# FIXME, this seems too ugly to find mrc_config.h
include_directories(${CMAKE_CURRENT_BINARY_DIR}/src/libmrc/include)

add_subdirectory(src/libmrc)
add_subdirectory(src)
