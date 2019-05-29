
add_library(kg INTERFACE)
target_compile_features(kg INTERFACE cxx_std_11)
# FIXME, hack for mrc_common.h
target_include_directories(kg INTERFACE include ${CMAKE_SOURCE_DIR}/src/libmrc/include)

if (PSC_HAVE_ADIOS2)
  target_link_libraries(kg INTERFACE adios2::adios2)
endif()

add_subdirectory(testing)
