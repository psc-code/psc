
add_library(psc)
target_gtensor_sources(psc PRIVATE
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

  vpic/psc_particles_vpic.cxx
  vpic/psc_fields_vpic.cxx
  vpic/psc_vpic_bits.cxx
  vpic/vpic_base.cxx
  )
target_include_directories(psc PUBLIC ../include)
target_link_libraries(psc PUBLIC kg mrc)
target_link_libraries(psc PUBLIC MPI::MPI_CXX)
target_link_libraries(psc PUBLIC gtensor::gtensor)
if(USE_CUDA)
  target_compile_definitions(psc PUBLIC GTENSOR_USE_THRUST)
endif()
target_compile_features(psc PUBLIC cxx_std_11)

if (USE_CUDA)
  target_gtensor_sources(psc PRIVATE
    cuda/cuda_base.cxx

    cuda/psc_particles_cuda.cxx
    cuda/cuda_mparticles_base.cxx
    cuda/cuda_mparticles.cxx

    cuda/psc_fields_cuda.cxx

    cuda/cuda_push_mprts_yz.cxx

    cuda/push_fields_cuda_impl.cxx
    cuda/cuda_push_fields_yz.cxx

    cuda/bnd_particles_cuda_impl.cxx
    cuda/cuda_bndp.cxx
    cuda/cuda_bndp_scan.cxx
    cuda/cuda_bndp_spine.cxx

    cuda/cuda_bnd.cxx
    cuda/bnd_cuda_3_impl.cxx

    cuda/psc_bnd_fields_cuda.cxx

    cuda/cuda_moments_1st_nc.cxx

    cuda/collision_cuda_impl.cxx

    cuda/cuda_heating.cxx

    cuda/mem.cxx
    )
  target_include_directories(psc PUBLIC cuda)
endif()

if (PSC_HAVE_RMM)
  target_link_libraries(psc PUBLIC rmm::rmm)
endif()

if (PSC_BOUNDS_CHECK)
  target_compile_definitions(psc PUBLIC BOUNDS_CHECK)
endif()

if (USE_VPIC)
  add_library(VPIC::VPIC INTERFACE IMPORTED)
  set_target_properties(VPIC::VPIC PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "/Users/kai/build/vpic/src"
    INTERFACE_LINK_LIBRARIES "/Users/kai/build/vpic/build/libvpic.a"
  )

  target_link_libraries(psc PUBLIC VPIC::VPIC)
endif()

if (BUILD_TESTING)
  add_subdirectory(tests)
endif()
