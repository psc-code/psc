
#pragma once

enum DEPOSIT {
  DEPOSIT_VB_2D,
  DEPOSIT_VB_3D,
};

enum CURRMEM {
  CURRMEM_SHARED,
  CURRMEM_GLOBAL,
};

// ----------------------------------------------------------------------
// cuda_push_mprts

template<typename BS>
struct cuda_mparticles;

template<typename Config>
struct CudaPushParticles_
{
  using BS = typename Config::Bs;
  using Current = typename Config::Current;
  using CudaMparticles = cuda_mparticles<BS>;
  
  static void push_mprts_yz(CudaMparticles* cmprts, struct cuda_mfields *cmflds);
  static void push_mprts_xyz(CudaMparticles* cmprts, struct cuda_mfields *cmflds);

  template<typename IP, enum DEPOSIT DEPOSIT, enum CURRMEM CURRMEM>
  static void push_mprts_yz_reorder(CudaMparticles* cmprts, struct cuda_mfields *cmflds);

  template<bool REORDER, typename OPT_IP, enum DEPOSIT DEPOSIT, enum CURRMEM CURRMEM>
  static void push_mprts_ab(CudaMparticles* cmprts, struct cuda_mfields *cmflds);
};
