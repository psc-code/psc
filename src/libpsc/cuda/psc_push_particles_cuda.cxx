
#include "cuda_iface.h"
#include "psc_push_particles_private.h"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"

#include <string.h>

struct CurrentGlobal : std::true_type {};
struct CurrentShared : std::false_type {};

template<typename CURRENT>
struct Config
{
  using Current = CURRENT;
};

// ======================================================================
// psc_push_particles: "1vb_4x4_cuda"

template<typename Config>
class PushParticlesCuda1vb
{
public:
  static void push_mprts_yz(struct psc_push_particles *push,
			    struct psc_mparticles *mprts,
			    struct psc_mfields *mflds_base)
  {
    // it's difficult to convert mprts because of the ordering constraints (?)
    assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);
    
    mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(EX, EX + 6);
    struct cuda_mparticles *cmprts = mparticles_cuda_t(mprts)->cmprts();
    int bs[3] = { BS::x::value, BS::y::value, BS::z::value };
    cuda_push_mprts_yz(cmprts, mf->cmflds, bs, false, false, Config::Current::value);
    mf.put_as(mflds_base, JXI, JXI + 3);
  }
};

template<typename Config>
class PushParticlesCuda1vbec3d
{
public:
  static void push_mprts_yz(struct psc_push_particles *push,
			    struct psc_mparticles *mprts,
			    struct psc_mfields *mflds_base)
  {
    /* it's difficult to convert mprts due to ordering constraints (?) */
    assert(strcmp(psc_mparticles_type(mprts), "cuda") == 0);
    
    mfields_cuda_t mf = mflds_base->get_as<mfields_cuda_t>(EX, EX + 6);
    struct cuda_mparticles *cmprts = mparticles_cuda_t(mprts)->cmprts();
    int bs[3] = { BS::x::value, BS::y::value, BS::z::value };
    cuda_push_mprts_yz(cmprts, mf->cmflds, bs, true, true, Config::Current::value);
    mf.put_as(mflds_base, JXI, JXI + 3);
  }
};

// ----------------------------------------------------------------------
// psc_push_particles: subclass "1vb_cuda"

struct psc_push_particles_ops_cuda : psc_push_particles_ops {
  psc_push_particles_ops_cuda() {
    name                  = "1vb_cuda";
    push_mprts_yz         = PushParticlesCuda1vb<Config<CurrentShared>>::push_mprts_yz;
    mp_flags              = MP_BLOCKSIZE_4X4X4;
  }
} psc_push_particles_1vb_cuda_ops;

#define MAKE_1VBEC3D_YZ(MP_BS, NAME, CONFIG)				\
  struct psc_push_particles_ops_ ## NAME :				\
    psc_push_particles_ops {						\
    psc_push_particles_ops_## NAME () {					\
      using PushParticlesCuda_t = PushParticlesCuda1vbec3d<CONFIG>;	\
      name                  = #NAME;					\
      push_mprts_yz         = PushParticlesCuda_t::push_mprts_yz;	\
      mp_flags              = MP_BS;					\
    }									\
  } psc_push_particles_## NAME ##_ops;


MAKE_1VBEC3D_YZ(MP_BLOCKSIZE_4X4X4, 1vbec3d_cuda     , Config<CurrentShared>);
MAKE_1VBEC3D_YZ(MP_BLOCKSIZE_4X4X4, 1vbec3d_gmem_cuda, Config<CurrentGlobal>);

