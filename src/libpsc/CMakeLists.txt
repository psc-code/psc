
add_library(psc
  bk_mparticles.cxx
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

  psc_output_fields/psc_output_fields_item.cxx
  psc_output_fields/psc_output_fields_item_jeh.cxx
  psc_output_fields/psc_output_fields_item_moments_1st.cxx
  psc_output_fields/psc_output_fields_item_moments_1st_nc_double.cxx
  psc_output_fields/psc_output_fields_item_moments_1st_nc_single.cxx
  psc_output_fields/psc_output_fields_item_moments_2nd_nc_double.cxx

  psc_output_particles/psc_output_particles.cxx
  psc_output_particles/psc_output_particles_ascii.cxx
  psc_output_particles/psc_output_particles_hdf5_impl.cxx
  psc_output_particles/psc_output_particles_none.cxx
  )
target_include_directories(psc PUBLIC ../include)
target_link_libraries(psc PUBLIC mrc)
target_compile_features(psc PUBLIC cxx_std_11)

if (USE_CUDA)
  add_subdirectory(cuda)
  target_link_libraries(psc PRIVATE psccuda)
  target_include_directories(psc PRIVATE vpic) # FIXME!
endif()

if (USE_VPIC)
  add_subdirectory(vpic)
  target_link_libraries(psc PUBLIC pscvpic)
endif()

add_subdirectory(tests)
