
add_subdirectory(libpsc)

macro(add_psc_executable name)
  add_executable(${name} ${name}.cxx)
  target_link_libraries(${name} psc)
endmacro(add_psc_executable)

add_psc_executable(psc_bubble_yz)
add_psc_executable(psc_flatfoil_yz)
add_psc_executable(psc_test_io_2)
add_executable(psc_test_io psc_test_io.cxx psc_test_io_xdmf.c)
target_link_libraries(psc_test_io psc)

if(USE_VPIC)
  add_psc_executable(psc_harris_xz)
endif()
