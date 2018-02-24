
#include "cuda_iface.h"
#include "psc_push_particles_private.h"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"

#include "push_particles.hxx"

#include <string.h>

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
  void push_mprts_yz(struct psc_mparticles *mprts,
		     struct psc_mfields *mflds_base) override
  {
    /* it's difficult to convert mprts due to ordering constraints (?) */
    assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);
    
    mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(EX, EX + 6);
    struct cuda_mparticles *cmprts = mparticles_cuda_t(mprts)->cmprts();
    int bs[3] = { BS::x::value, BS::y::value, BS::z::value };
    cuda_push_mprts_yz(cmprts, mf->cmflds, bs, Config::Ip::value, Config::Deposit::value,
		       Config::Current::value);
    mf.put_as(mflds_base, JXI, JXI + 3);
  }
};

// ----------------------------------------------------------------------
// psc_push_particles: subclasses

#define MAKE_1VBEC3D_YZ(MP_BS, NAME, CONFIG)				\
  struct psc_push_particles_ops_ ## NAME :				\
    psc_push_particles_ops {						\
    psc_push_particles_ops_## NAME () {					\
      using PushParticles_t = PushParticlesCuda<CONFIG>;		\
      using PushParticlesWrapper_t = PushParticlesWrapper<PushParticles_t>; \
      name          = #NAME;						\
      size          = PushParticlesWrapper_t::size;			\
      setup         = PushParticlesWrapper_t::setup;			\
      destroy       = PushParticlesWrapper_t::destroy;			\
      mp_flags      = MP_BS;						\
    }									\
  } psc_push_particles_## NAME ##_ops;

using Config1vb = Config<IpRegular, DepositVb2d, CurrentShared>;
MAKE_1VBEC3D_YZ(MP_BLOCKSIZE_4X4X4, 1vb_cuda, Config1vb);
using Config1vbec3d = Config<IpEc, DepositVb3d, CurrentShared>;
MAKE_1VBEC3D_YZ(MP_BLOCKSIZE_4X4X4, 1vbec3d_cuda, Config1vbec3d);
using Config1vbec3dGmem = Config<IpEc, DepositVb3d, CurrentShared>;
MAKE_1VBEC3D_YZ(MP_BLOCKSIZE_4X4X4, 1vbec3d_gmem_cuda, Config1vbec3dGmem);

