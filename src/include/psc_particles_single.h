
#ifndef PSC_PARTICLE_SINGLE_H
#define PSC_PARTICLE_SINGLE_H

#include "psc_particles_private.h"

#include "particles_simple.hxx"
#include "particles_traits.hxx"

using particle_single_t = psc_particle<float>;

using MparticlesSingle = Mparticles<particle_single_t>;

template<>
struct Mparticles_traits<MparticlesSingle>
{
  static constexpr const char* name = "single";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#define particle_single_qni_wni(prt) ((prt)->qni_wni_)

#endif
