
set(MPIEXEC_COMMAND
  ${MPIEXEC_EXECUTABLE}

  # ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
)

macro(add_psc_test name)
  add_executable(${name})
  target_gtensor_sources(${name} PRIVATE ${name}.cxx)
  target_link_libraries(${name} psc GTest::gtest)

  if(USE_GTEST_DISCOVER_TESTS)
    gtest_discover_tests(${name} EXEC_WRAPPER ${MPIEXEC_COMMAND})
  else()
    gtest_add_tests(TARGET ${name} EXEC_WRAPPER ${MPIEXEC_COMMAND})
  endif()
endmacro()

add_psc_test(test_grid)
add_psc_test(test_grid_2)
add_psc_test(test_rng)
add_psc_test(test_mparticles_cuda)
add_psc_test(test_mparticles)
add_psc_test(test_output_particles)
add_psc_test(test_reflective_bcs)
add_psc_test(test_reflective_bcs_integration)
add_psc_test(test_mfields)
add_psc_test(test_mfields_cuda)
add_psc_test(test_bnd)
add_psc_test(test_bnd_prt_inflow)
add_psc_test(test_deposit)
add_psc_test(test_current_deposition)
add_psc_test(test_push_particles)
add_psc_test(test_push_particles_2)
add_psc_test(test_push_fields)
add_psc_test(test_moments)
add_psc_test(test_collision)

if(USE_CUDA AND NOT USE_VPIC)
  add_psc_test(test_collision_cuda)
endif()

add_psc_test(test_balance)
add_psc_test(TestUniqueIdGenerator)
add_psc_test(test_mfields_io)
