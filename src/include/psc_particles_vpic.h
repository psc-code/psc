
#ifndef PSC_PARTICLES_VPIC_H
#define PSC_PARTICLES_VPIC_H

#include "psc_particles_private.h"
#include "psc_particles_single.h"

#include "particles.hxx"
#include "particles_traits.hxx"

#include "../libpsc/vpic/vpic_iface.h" // FIXME path

struct psc_mparticles_vpic {
  Particles *vmprts;
  Simulation *sim;
};

using mparticles_vpic_t = mparticles_base<psc_mparticles_vpic>;

template<>
struct mparticles_traits<mparticles_vpic_t>
{
  static constexpr const char* name = "vpic";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#define psc_mparticles_vpic(mprts)({					\
      assert((struct psc_mparticles_ops *) mprts->obj.ops == &psc_mparticles_vpic_ops); \
      mrc_to_subobj(mprts, struct psc_mparticles_vpic);			\
})

#endif
