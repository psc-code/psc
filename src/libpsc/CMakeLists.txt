
add_library(psc
  bk_mparticles.cxx
  integrate.cxx
  mfields.cxx
  mparticles.cxx
  psc.cxx
  psc_checkpoint.cxx
  psc_fields_c.cxx
  psc_fields_single.cxx
  psc_main.cxx
  psc_method.cxx
  psc_method_default.cxx
  psc_particles_double.cxx
  psc_particles_single.cxx
  psc_particles_single_by_kind.cxx
  psc_stats.cxx
  rngpool.cxx

  psc_push_particles/psc_push_particles.cxx
  psc_push_particles/psc_push_particles_1st.cxx
  psc_push_particles/psc_push_particles_2nd.cxx
  psc_push_particles/psc_push_particles_generic_c.cxx
  psc_push_particles/push_part_1st_xz.cxx
  psc_push_particles/push_part_1st_yz.cxx
  psc_push_particles/push_part_single2_1vb_yz.cxx
  psc_push_particles/push_part_xy.cxx
  psc_push_particles/push_part_xyz.cxx
  psc_push_particles/push_part_xz.cxx
  psc_push_particles/push_part_y.cxx
  psc_push_particles/push_part_yz.cxx
  psc_push_particles/push_part_z.cxx
  psc_push_particles/push_part_1.cxx
  psc_push_particles/1vb/psc_push_particles_1vb_double.cxx
  psc_push_particles/1vb/psc_push_particles_1vb_single.cxx
  psc_push_particles/1vb/psc_push_particles_1vbec_double.cxx
  psc_push_particles/1vb/psc_push_particles_1vbec_double_1.cxx
  psc_push_particles/1vb/psc_push_particles_1vbec_double_xyz.cxx
  psc_push_particles/1vb/psc_push_particles_1vbec_double_yz.cxx
  psc_push_particles/1vb/psc_push_particles_1vbec_single.cxx
  psc_push_particles/1vb/psc_push_particles_1vbec_single_1.cxx
  psc_push_particles/1vb/psc_push_particles_1vbec_single_xyz.cxx
  psc_push_particles/1vb/psc_push_particles_1vbec_single_xyz_xz.cxx
  psc_push_particles/1vb/psc_push_particles_1vbec_single_yz.cxx

  psc_push_fields/psc_marder.cxx
  psc_push_fields/psc_marder_c.cxx
  psc_push_fields/psc_marder_single.cxx
  psc_push_fields/psc_push_fields.cxx
  psc_push_fields/psc_push_fields_impl.cxx
  psc_push_fields/psc_push_fields_none.cxx

  psc_sort/psc_sort.cxx
  psc_sort/psc_sort_impl.cxx
  psc_sort/psc_sort_none.cxx

  psc_randomize/psc_randomize.cxx
  psc_randomize/psc_randomize_double.cxx
  psc_randomize/psc_randomize_none.cxx

  psc_collision/psc_collision.cxx
  psc_collision/psc_collision_double.cxx
  psc_collision/psc_collision_none.cxx
  psc_collision/psc_collision_single.cxx

  psc_bnd/psc_bnd.cxx
  psc_bnd/psc_bnd_impl.cxx

  psc_bnd_fields/psc_bnd_fields.cxx
  psc_bnd_fields/psc_bnd_fields_impl.cxx
  psc_bnd_fields/psc_bnd_fields_none.cxx

  psc_bnd_particles/psc_bnd_particles.cxx
  psc_bnd_particles/psc_bnd_particles_impl.cxx

  psc_balance/psc_balance.cxx
  psc_balance/psc_balance_double.cxx
  psc_balance/psc_balance_single.cxx

  psc_checks/psc_checks.cxx
  psc_checks/psc_checks_1st_double.cxx
  psc_checks/psc_checks_1st_single.cxx
  psc_checks/psc_checks_2nd_double.cxx
  psc_checks/psc_checks_2nd_single.cxx

  psc_event_generator/psc_event_generator.cxx
  psc_event_generator/psc_event_generator_none.cxx

  psc_diag/psc_diag.cxx
  psc_diag/psc_diag_item.cxx
  psc_diag/psc_diag_item_field_energy.cxx
  psc_diag/psc_diag_item_particle_energy.cxx

  psc_output_fields/psc_output_fields.cxx
  psc_output_fields/psc_output_fields_c.cxx
  psc_output_fields/psc_output_fields_item.cxx
  psc_output_fields/psc_output_fields_item_flds_c.cxx
  psc_output_fields/psc_output_fields_item_flds_single.cxx
  psc_output_fields/psc_output_fields_item_jeh.cxx
  psc_output_fields/psc_output_fields_item_moments_1st_double.cxx
  psc_output_fields/psc_output_fields_item_moments_1st_nc_double.cxx
  psc_output_fields/psc_output_fields_item_moments_1st_nc_single.cxx
  psc_output_fields/psc_output_fields_item_moments_1st_single.cxx
  psc_output_fields/psc_output_fields_item_moments_2nd_nc_double.cxx

  psc_output_fields_collection/psc_output_fields_collection.cxx

  psc_output_particles/psc_output_particles.cxx
  psc_output_particles/psc_output_particles_ascii.cxx
  psc_output_particles/psc_output_particles_hdf5_double.cxx
  psc_output_particles/psc_output_particles_hdf5_single.cxx
  psc_output_particles/psc_output_particles_none.cxx

  psc_heating/psc_heating.cxx
  psc_heating/psc_heating_single.cxx
  psc_heating/psc_heating_double.cxx
  psc_heating/psc_heating_spot.cxx

  psc_inject/psc_inject.cxx
  psc_inject/psc_inject_single.cxx
  psc_inject/psc_inject_double.cxx
  psc_inject/psc_target.cxx
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
