
cmake_minimum_required(VERSION 3.9)

foreach(policy
    CMP0074 # CMake 3.12
    )
  if(POLICY ${policy})
    cmake_policy(SET ${policy} NEW)
  endif()
endforeach()

project(libmrc)

find_package(MPI REQUIRED C)

find_package(HDF5 REQUIRED C HL)
if(HDF5_FOUND AND NOT TARGET HDF5::C)
  add_library(HDF5::C INTERFACE IMPORTED)
  set_target_properties(HDF5::C PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "${HDF5_C_DEFINTIONS}"
    INTERFACE_INCLUDE_DIRECTORIES "${HDF5_C_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${HDF5_C_LIBRARIES}"
    )
  set(HAVE_HDF5 1)
  if(HDF5_IS_PARALLEL)
    set(HAVE_HDF5_PARALLEL 1)
  endif(HDF5_IS_PARALLEL)
endif()
if(HDF5_FOUND AND NOT TARGET HDF5::HL)
  add_library(HDF5::HL INTERFACE IMPORTED)
  set_target_properties(HDF5::HL PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "${HDF5_C_DEFINTIONS}"
    INTERFACE_INCLUDE_DIRECTORIES "${HDF5_C_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${HDF5_C_HL_LIBRARIES};${HDF5_C_LIBRARIES}"
  )
endif()

# FIXME!!!
file(GLOB LIBMRC_SOURCES src/*)
file(GLOB LIBMRC_REMOVE src/mrc_ddc_mb.c src/mrc_fld_common.c src/mrc_mat_petsc.c
  src/mrc_ndarray_common.c src/mrc_ts_petsc.c src/mrc_vec_petsc.c
  src/mrc_perfetto.cxx)
list(REMOVE_ITEM LIBMRC_SOURCES ${LIBMRC_REMOVE})

add_library(mrc ${LIBMRC_SOURCES})
target_include_directories(mrc
  PUBLIC include .
)
target_link_libraries(mrc
  PUBLIC
    m
    MPI::MPI_C
    HDF5::HL
)
target_compile_features(mrc
  PRIVATE
    c_std_99
)

if (PSC_HAVE_NVTX)
  target_link_libraries(mrc PUBLIC CUDA::nvToolsExt)
  set(HAVE_NVTX 1)
endif()

if (PSC_USE_PERFETTO)
  target_sources(mrc PRIVATE src/mrc_perfetto.cxx)
  target_link_libraries(mrc PRIVATE perfetto::perfetto)
  target_compile_definitions(mrc PRIVATE HAVE_PERFETTO)
endif()

configure_file (
  "${PROJECT_SOURCE_DIR}/include/mrc_config.h.in"
  "${PROJECT_BINARY_DIR}/include/mrc_config.h"
  )
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)

add_subdirectory(tests)
