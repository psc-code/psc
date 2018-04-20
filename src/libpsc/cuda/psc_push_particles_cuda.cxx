
#include "cuda_iface.h"
#include "psc_push_particles_private.h"
#include "psc_particles_cuda.h"
#include "psc_fields_cuda.h"

#include "push_particles_cuda_impl.hxx"

#include <string.h>

// ----------------------------------------------------------------------
// psc_push_particles: subclasses

#define MAKE_1VBEC3D_YZ(NAME, CONFIG)				\
  struct psc_push_particles_ops_ ## NAME :				\
    psc_push_particles_ops {						\
    psc_push_particles_ops_## NAME () {					\
      using PushParticles_t = PushParticlesCuda<CONFIG>;		\
      using PushParticlesWrapper_t = PushParticlesWrapper<PushParticles_t>; \
      name          = #NAME;						\
      size          = PushParticlesWrapper_t::size;			\
      setup         = PushParticlesWrapper_t::setup;			\
      destroy       = PushParticlesWrapper_t::destroy;			\
    }									\
  } psc_push_particles_## NAME ##_ops;

MAKE_1VBEC3D_YZ(1vb_cuda, Config1vb);
MAKE_1VBEC3D_YZ(1vbec3d_cuda, Config1vbec3d);
MAKE_1VBEC3D_YZ(1vbec3d_gmem_cuda, Config1vbec3dGmem);

