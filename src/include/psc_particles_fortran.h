
#ifndef PSC_PARTICLE_FORTRAN_H
#define PSC_PARTICLE_FORTRAN_H

#include "psc_particles_private.h"

#include "particles_traits.hxx"

#define PTYPE PTYPE_FORTRAN
#include "psc_particle_buf_common.h"
#include "psc_particles_common.h"
#undef PTYPE

template<>
struct mparticles_traits<particle_fortran_t>
{
  static constexpr const char* name = "fortran";
  static MPI_Datatype mpi_dtype() { return MPI_DOUBLE; }
};

#endif
