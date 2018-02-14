
cmake_minimum_required (VERSION 3.9)

project(PSC)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
include_directories(${MPI_C_INCLUDE_PATH})
link_libraries(${MPI_C_LIBRARIES})
add_definitions(-DOMPI_SKIP_MPICXX)

find_package(HDF5 REQUIRED C HL)
include_directories(${HDF5_INCLUDE_DIRS})
link_libraries(${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

add_subdirectory(src/libmrc)
add_subdirectory(src)



