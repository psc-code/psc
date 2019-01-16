
#ifndef PSC_PARTICLE_SINGLE_H
#define PSC_PARTICLE_SINGLE_H

#include "particles_simple.hxx"
#include "particles_traits.hxx"

using MparticlesSingle = MparticlesSimple<ParticleSimple<float>>;

template<>
struct Mparticles_traits<MparticlesSingle>
{
  static constexpr const char* name = "single";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

#endif
