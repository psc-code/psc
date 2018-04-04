
#include "cuda_iface.h"
#include "psc_push_particles_private.h"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"

#include "push_particles_cuda_impl.hxx"

#include <string.h>

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

