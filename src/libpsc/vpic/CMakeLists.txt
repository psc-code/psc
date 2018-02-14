
add_library(pscvpic
  psc_fields_vpic.cxx
  psc_particles_vpic.cxx
  psc_push_fields_vpic.cxx
  psc_push_particles_vpic.cxx
  psc_sort_vpic.cxx
  psc_collision_vpic.cxx
  psc_marder_vpic.cxx
  psc_bnd_vpic.cxx
  psc_bnd_fields_vpic.cxx
  psc_bnd_particles_vpic.cxx
  psc_output_fields_item_vpic.cxx
  psc_method_vpic.cxx
  psc_vpic_bits.cxx
  simulation.cxx
  vpic_base.cxx
  vpic_mparticles.cxx
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

