
#ifndef PSC_PARTICLE_DOUBLE_H
#define PSC_PARTICLE_DOUBLE_H

#include "psc_particles_private.h"
#include "psc.h"

#include "particles_simple.hxx"
#include "particles_traits.hxx"

using particle_double_t = psc_particle<double>;

using MparticlesDouble = Mparticles<particle_double_t>;

template<>
struct Mparticles_traits<MparticlesDouble>
{
  static constexpr const char* name = "double";
  static MPI_Datatype mpi_dtype() { return MPI_DOUBLE; }
};

#endif
