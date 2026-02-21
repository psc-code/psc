
add_subdirectory(kg)
add_subdirectory(libpsc)

macro(add_psc_executable name)
  add_executable(${name})
  target_gtensor_sources(${name} PRIVATE ${name}.cxx)
  target_link_libraries(${name} psc)
endmacro(add_psc_executable)

if(NOT USE_CUDA)
  add_psc_executable(psc_bgk)
endif()

add_psc_executable(psc_bubble_yz)
add_psc_executable(psc_flatfoil_yz)
add_psc_executable(psc_flatfoil_yz_small)
add_psc_executable(psc_whistler)
add_psc_executable(psc_harris_yz)
add_psc_executable(psc_2d_shock)
add_psc_executable(psc_radiation)
add_psc_executable(psc_shock)

if(NOT USE_CUDA)
  install(
    TARGETS psc_bgk
    RUNTIME DESTINATION bin
  )
endif()

install(
  TARGETS psc_bubble_yz psc_flatfoil_yz psc_whistler psc_shock
  RUNTIME DESTINATION bin
)
