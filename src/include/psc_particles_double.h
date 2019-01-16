
#ifndef PSC_PARTICLE_DOUBLE_H
#define PSC_PARTICLE_DOUBLE_H

#include "particles_simple.hxx"
#include "particles_traits.hxx"

using MparticlesDouble = MparticlesSimple<ParticleSimple<double>>;

template<>
struct Mparticles_traits<MparticlesDouble>
{
  static constexpr const char* name = "double";
  static MPI_Datatype mpi_dtype() { return MPI_DOUBLE; }
};

#endif
