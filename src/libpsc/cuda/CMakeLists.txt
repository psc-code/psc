
add_library(psccuda
  cuda_base.cu
  
  psc_particles_cuda.cu
  cuda_mparticles_base.cu
  cuda_mparticles.cu
  cuda_mparticles_iface.cu

  psc_fields_cuda.cxx
  cuda_mfields.cu
  cuda_mfields_iface.cu

  cuda_push_mprts_yz.cu

  cuda_push_fields_yz.cu
  cuda_push_fields_xyz.cu

  bnd_particles_cuda_impl.cu
  cuda_bndp.cu
  cuda_bndp_scan.cu
  cuda_bndp_spine.cu
  
  cuda_bnd.cu
  bnd_cuda_3_impl.cu

  psc_bnd_fields_cuda.cxx
  cuda_mfields_bnd.cu

  psc_output_fields_item_flds_cuda.cxx
  psc_output_fields_item_moments_1st_nc_cuda.cxx
  cuda_moments_1st_nc.cu

  collision_cuda_impl.cu

  psc_heating_cuda.cxx
  cuda_heating.cu
  )
target_include_directories(psccuda
  PUBLIC
  ${PROJECT_SOURCE_DIR}/src/include
  ${PROJECT_SOURCE_DIR}/src/libpsc/cuda
  )
target_link_libraries(psccuda
  mrc
  )
target_compile_features(psccuda
  PUBLIC cxx_std_11
  )

macro(add_cuda_test name)
  if(GTEST_FOUND)
    add_executable(${name} tests/${name}.cu)
    target_link_libraries(${name} psccuda GTest::GTest GTest::Main)
    gtest_discover_tests(${name})
  endif()
endmacro()

add_cuda_test(range_unittest)
add_cuda_test(test_cuda_base)
add_cuda_test(test_cuda_accel)
add_cuda_test(test_cuda_bndp)
add_cuda_test(test_cuda_bndp_2)
add_cuda_test(test_cuda_mfields)
add_cuda_test(test_cuda_mparticles)

