
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

