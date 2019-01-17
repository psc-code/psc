
add_library(psc
  integrate.cxx
  mfields.cxx
  mparticles.cxx
  psc.cxx
  psc_checkpoint.cxx
  psc_fields_c.cxx
  psc_fields_single.cxx
  psc_particles_impl.cxx
  psc_stats.cxx
  rngpool.cxx

  psc_collision/psc_collision_impl.cxx

  psc_balance/psc_balance.cxx

  psc_diag/psc_diag.cxx
  psc_diag/psc_diag_item.cxx
  psc_diag/psc_diag_item_field_energy.cxx
  psc_diag/psc_diag_item_particle_energy.cxx
  )
target_include_directories(psc PUBLIC ../include)
target_link_libraries(psc PUBLIC mrc)
target_compile_features(psc PUBLIC cxx_std_11)

if (USE_CUDA)
  add_subdirectory(cuda)
  target_link_libraries(psc PRIVATE psccuda)
  target_include_directories(psc PRIVATE vpic) # FIXME!
endif()

add_subdirectory(vpic)
target_link_libraries(psc PUBLIC pscvpic)

add_subdirectory(tests)
