
cmake_minimum_required (VERSION 3.9)

add_subdirectory(src/libmrc)

project(PSC)

option (USE_CUDA "Build CUDA components" OFF) 

if(USE_CUDA)
  enable_language(CUDA)
endif()

find_package(MPI REQUIRED)
include_directories(${MPI_C_INCLUDE_PATH})
link_libraries(${MPI_C_LIBRARIES})
add_definitions(-DOMPI_SKIP_MPICXX)

find_package(HDF5 REQUIRED C HL) # FIXME, don't make required / how about parallel
include_directories(${HDF5_INCLUDE_DIRS})
link_libraries(${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

configure_file (
  "${PROJECT_SOURCE_DIR}/src/include/psc_config.h.in"
  "${PROJECT_BINARY_DIR}/src/include/psc_config.h"
  )
include_directories(${CMAKE_CURRENT_BINARY_DIR}/src/include)

add_subdirectory(src)



