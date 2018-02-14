
add_library(pscvpic
  psc_vpic_bits.cxx
  simulation.cxx
  vpic_base.cxx
  vpic_mparticles.cxx
  )

target_include_directories(pscvpic
  PRIVATE ${PROJECT_SOURCE_DIR}/src/include
  )
target_link_libraries(pscvpic
  mrc
  VPIC::VPIC
  )
target_compile_features(pscvpic
  PUBLIC cxx_std_11
  )

