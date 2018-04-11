
#pragma once

// ----------------------------------------------------------------------
// cuda_push_mprts

template<typename BS>
struct cuda_mparticles;

template<typename Config>
struct CudaPushParticles_
{
  using BS = typename Config::Bs;
  
  static void push_mprts_yz(cuda_mparticles<BS>* cmprts, struct cuda_mfields *cmflds);
  static void push_mprts_xyz(cuda_mparticles<BS>* cmprts, struct cuda_mfields *cmflds);
};
