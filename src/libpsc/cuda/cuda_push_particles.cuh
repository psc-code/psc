
#pragma once

enum DEPOSIT {
  DEPOSIT_VB_2D,
  DEPOSIT_VB_3D,
};

// ----------------------------------------------------------------------
// cuda_push_mprts

template<typename BS>
struct cuda_mparticles;

template<typename Config>
struct CudaPushParticles_
{
  using BS = typename Config::Bs;
  using CudaMparticles = cuda_mparticles<BS>;
  
  static void push_mprts(CudaMparticles* cmprts, struct cuda_mfields *cmflds);

  template<bool REORDER>
  static void push_mprts_ab(CudaMparticles* cmprts, struct cuda_mfields *cmflds);
};
