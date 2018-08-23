
add_library(pscvpic
  psc_fields_vpic.cxx
  psc_particles_vpic.cxx
  psc_vpic_bits.cxx
  vpic_base.cxx
  )

target_include_directories(pscvpic
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/include
    ${PROJECT_SOURCE_DIR}/src/libpsc/vpic
  )
target_link_libraries(pscvpic
  mrc
  VPIC::VPIC
  )
target_compile_features(pscvpic
  PUBLIC cxx_std_11
  )

