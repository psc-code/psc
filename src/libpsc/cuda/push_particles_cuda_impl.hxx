
#pragma once

#include "push_particles.hxx"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"
#include "cuda_iface.h"

struct IpEc : std::true_type {};
struct IpRegular : std::false_type {};

struct DepositVb3d : std::true_type {};
struct DepositVb2d : std::false_type {};

struct CurrentGlobal : std::true_type {};
struct CurrentShared : std::false_type {};

template<typename IP, typename DEPOSIT, typename CURRENT>
struct Config
{
  using Ip = IP;
  using Deposit = DEPOSIT;
  using Current = CURRENT;
};

// ======================================================================
// psc_push_particles: "1vb_4x4_cuda"

template<typename Config>
class PushParticlesCuda : PushParticlesBase
{
public:
  void push_mprts(MparticlesCuda& mprts, MfieldsCuda& mflds)
  {
    int bs[3] = { BS144::x::value, BS144::y::value, BS144::z::value };
    cuda_push_mprts_yz(mprts.cmprts(), mflds.cmflds, bs, Config::Ip::value, Config::Deposit::value,
		       Config::Current::value);
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

