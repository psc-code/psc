
add_library(psccuda
  cuda_base.cu
  
  psc_particles_cuda.cxx
  cuda_mparticles_base.cu
  cuda_mparticles.cu
  cuda_mparticles_iface.cu

  psc_fields_cuda.cxx
  cuda_mfields.cu
  cuda_mfields_iface.cu

  psc_push_particles_cuda.cxx
  cuda_push_mprts_yz.cu

  psc_push_fields_cuda.cxx
  cuda_push_fields_yz.cu

  psc_marder_cuda.cxx

  psc_bnd_particles_cuda.cu
  cuda_bndp.cu
  cuda_bndp_scan.cu
  cuda_bndp_spine.cu
  
  psc_bnd_cuda.cxx
  cuda_bnd.cu

  psc_bnd_fields_cuda.cxx
  cuda_mfields_bnd.cu

  psc_output_fields_item_flds_cuda.cxx
  psc_output_fields_item_moments_1st_nc_cuda.cxx
  cuda_moments_1st_nc.cu

  psc_heating_cuda.cxx
  cuda_heating.cu
  
  psc_inject_cuda.cxx
  )
target_include_directories(psccuda
  PRIVATE ${PROJECT_SOURCE_DIR}/src/include
  )
target_link_libraries(psccuda
  mrc
  )
target_compile_features(psccuda
  PUBLIC cxx_std_11
  )

# FIXME, only hackery here...
set(GTEST_ROOT /Users/kai/build/googletest/googletest)

macro(add_cuda_test name)
add_executable(${name} tests/${name}.cu)
target_include_directories(${name} PRIVATE ${PROJECT_SOURCE_DIR}/src/include ${PROJECT_SOURCE_DIR}/src/libpsc/cuda)
target_link_libraries(${name} psccuda)
target_include_directories(${name} PRIVATE ${GTEST_ROOT}/include)
target_link_libraries(${name} ${GTEST_ROOT}/src/gtest-all.o ${GTEST_ROOT}/src/gtest_main.o)
endmacro()

add_cuda_test(test_cuda_accel)
add_cuda_test(test_cuda_bndp)
add_cuda_test(test_cuda_mfields)
add_cuda_test(test_cuda_mparticles_1)

