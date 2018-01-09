
#ifndef PSC_PARTICLES_CUDA_H
#define PSC_PARTICLES_CUDA_H

#include "psc_particles_private.h"
#include "psc_particles_single.h"

#include "particles_traits.hxx"

#include "psc_particle_buf_cuda.h"

#define PTYPE PTYPE_CUDA
#include "psc_particles_common.h"
#undef PTYPE

template<>
struct mparticles_traits<particle_cuda_t>
{
  static constexpr const char* name = "cuda";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};



#endif
