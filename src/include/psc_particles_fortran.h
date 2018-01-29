
#ifndef PSC_PARTICLE_FORTRAN_H
#define PSC_PARTICLE_FORTRAN_H

#include "psc_particles_private.h"

#include "particles.hxx"
#include "particles_traits.hxx"

using particle_fortran_real_t = double;

struct particle_fortran_t
{
  using real_t = particle_fortran_real_t;

  real_t xi, yi, zi;
  real_t pxi, pyi, pzi;
  real_t qni;
  real_t mni;
  real_t cni;
  real_t lni;
  real_t wni;
};

#define psc_mparticles_fortran(mprts) mrc_to_subobj(mprts, struct psc_mparticles_fortran)

using psc_mparticles_fortran = psc_mparticles_<particle_fortran_t>;
using mparticles_fortran_t = mparticles<psc_mparticles_fortran>;

template<>
struct mparticles_traits<mparticles_fortran_t>
{
  static constexpr const char* name = "fortran";
  static MPI_Datatype mpi_dtype() { return MPI_DOUBLE; }
};

#endif
