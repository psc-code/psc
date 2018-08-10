
#pragma once

#include "push_particles.hxx"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"
#include "cuda_push_particles.cuh"
#include "interpolate.hxx"

struct DepositVb3d : std::integral_constant<int, DEPOSIT_VB_3D> {};
struct DepositVb2d : std::integral_constant<int, DEPOSIT_VB_2D> {};

struct CurrmemGlobal;
struct CurrmemShared;

template<typename DIM, typename BS, typename IP, typename DEPOSIT, typename CURRMEM>
struct CudaPushpConfig
{
  using dim = DIM;
  using Bs = BS;
  using Ip = IP;
  using Deposit = DEPOSIT;
  using Currmem = CURRMEM;
};

template<typename dim, typename BS>
using CudaConfig1vbec3d = CudaPushpConfig<dim, BS, opt_ip_1st_ec, DepositVb3d, CurrmemShared>;

template<typename dim>
using CudaConfig1vb = CudaPushpConfig<dim, BS144, opt_ip_1st, DepositVb2d, CurrmemShared>;

template<typename dim, typename BS>
using CudaConfig1vbec3dGmem = CudaPushpConfig<dim, BS, opt_ip_1st_ec, DepositVb3d, CurrmemGlobal>;

// ======================================================================
// psc_push_particles: "1vb_4x4_cuda"

template<typename Config>
class PushParticlesCuda : PushParticlesBase
{
public:
  using BS = typename Config::Bs;
  using Mparticles = MparticlesCuda<BS>;
  using MfieldsState = MfieldsStateCuda;
  
  void push_mprts(Mparticles& mprts, MfieldsState& mflds)
  {
    CudaPushParticles_<Config>::push_mprts(mprts.cmprts(), mflds.cmflds());
  }
  
  void push_mprts_yz(MparticlesBase& mprts_base, MfieldsStateBase& mflds_base) override
  {
    auto& mflds = mflds_base.get_as<MfieldsState>(EX, EX + 6);
    auto& mprts = mprts_base.get_as<Mparticles>();
    push_mprts(mprts, mflds);
    mprts_base.put_as(mprts);
    mflds_base.put_as(mflds, JXI, JXI + 3);
  }
};

