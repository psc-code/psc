
#pragma once

#include "push_particles.hxx"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"
#include "cuda_push_particles.cuh"
#include "interpolate.hxx"

struct IpEc { using type = opt_ip_1st_ec; };
struct IpRegular { using type = opt_ip_1st; };

struct DepositVb3d : std::true_type {};
struct DepositVb2d : std::false_type {};

struct CurrmemGlobal;
struct CurrmemShared;

template<typename BS, typename IP, typename DEPOSIT, typename CURRMEM>
struct PushParticlesConfig
{
  using Bs = BS;
  using Ip = IP;
  using Deposit = DEPOSIT;
  using Currmem = CURRMEM;
};

using Config1vb = PushParticlesConfig<BS144, IpRegular, DepositVb2d, CurrmemShared>;
using Config1vbec3d = PushParticlesConfig<BS144, IpEc, DepositVb3d, CurrmemShared>;
using Config1vbec3dGmem = PushParticlesConfig<BS144, IpEc, DepositVb3d, CurrmemShared>;

// ======================================================================
// psc_push_particles: "1vb_4x4_cuda"

template<typename Config>
class PushParticlesCuda : PushParticlesBase
{
public:
  void push_mprts(MparticlesCuda& mprts, MfieldsCuda& mflds)
  {
    CudaPushParticles_<Config>::push_mprts_yz(mprts.cmprts(), mflds.cmflds);
  }
  
  void push_mprts_yz(PscMparticlesBase mprts_base, PscMfieldsBase mflds_base) override
  {
    auto& mflds = mflds_base->get_as<MfieldsCuda>(EX, EX + 6);
    auto& mprts = mprts_base->get_as<MparticlesCuda>();
    push_mprts(mprts, mflds);
    mprts_base->put_as(mprts);
    mflds_base->put_as(mflds, JXI, JXI + 3);
  }
};

