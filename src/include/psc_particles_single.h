
#ifndef PSC_PARTICLE_SINGLE_H
#define PSC_PARTICLE_SINGLE_H

#include "psc_particles_private.h"

#include "particles.hxx"
#include "particles_traits.hxx"

using particle_single_real_t = float;

struct particle_single_t : psc_particle<particle_single_real_t> {};

using psc_particle_single_buf_t = psc_particle_buf<particle_single_t>;

#define psc_mparticles_single(mprts) mrc_to_subobj(mprts, struct psc_mparticles_single)

#define PTYPE PTYPE_SINGLE
#include "psc_particles_common.h"
#undef PTYPE

using mparticles_single_t = mparticles<psc_mparticles_single>;

template<>
struct mparticles_traits<mparticles_single_t>
{
  static constexpr const char* name = "single";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

// FIXME, needs to eventually
template<>
struct mparticles_traits<particle_single_t>
{
  static constexpr const char* name = "single";
  static MPI_Datatype mpi_dtype() { return MPI_FLOAT; }
};

// can't do this as inline function since struct psc isn't known yet
#define particle_single_qni_div_mni(p) ({			\
      particle_single_real_t rv;				\
      rv = ppsc->kinds[p->kind_].q / ppsc->kinds[p->kind_].m;	\
      rv;							\
    })

#define particle_single_qni(p) ({				\
      particle_single_real_t rv;				\
      rv = ppsc->kinds[p->kind_].q;				\
      rv;							\
    })

#define particle_single_mni(p) ({				\
      particle_single_real_t rv;				\
      rv = ppsc->kinds[p->kind_].m;				\
      rv;							\
    })

#define particle_single_wni(p) ({				\
      particle_single_real_t rv;				\
      rv = (p)->qni_wni / ppsc->kinds[(p)->kind_].q;		\
      rv;							\
    })

#define particle_single_qni_wni(prt) ((prt)->qni_wni)

#endif
