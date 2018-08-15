
#ifndef PSC_PARTICLES_VPIC_H
#define PSC_PARTICLES_VPIC_H

#include "psc_particles_private.h"
#include "psc_particles_single.h"

#include "particles.hxx"
#include "particles_traits.hxx"
#include "../libpsc/vpic/psc_vpic_bits.h"

#include "../libpsc/vpic/vpic_iface.h"

using MparticlesVpic = MparticlesVpic_<Particles>;

template<>
struct Mparticles_traits<MparticlesVpic>
{
  static constexpr const char* name = "vpic";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

// ======================================================================

#endif
