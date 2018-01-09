
#ifndef PSC_PARTICLE_DOUBLE_H
#define PSC_PARTICLE_DOUBLE_H

#include "psc_particles_private.h"
#include "psc.h"

#include "particles.hxx"
#include "particles_traits.hxx"

#define PTYPE PTYPE_DOUBLE
#include "psc_particle_buf_common.h"
#include "psc_particles_common.h"
#undef PTYPE

using mparticles_double_t = mparticles<psc_particle_double_buf_t, psc_mparticles_double>;

template<>
struct mparticles_traits<mparticles_double_t>
{
  static constexpr const char* name = "double";
  static MPI_Datatype mpi_dtype() { return MPI_DOUBLE; }
};

template<>
struct mparticles_traits<particle_double_t>
{
  static constexpr const char* name = "double";
  static MPI_Datatype mpi_dtype() { return MPI_DOUBLE; }
};

// can't do this as inline function since struct psc isn't known yet
#define particle_double_qni_div_mni(p) ({			\
      particle_double_real_t rv;				\
      rv = ppsc->kinds[p->kind_].q / ppsc->kinds[p->kind_].m;	\
      rv;							\
    })

#define particle_double_qni(p) ({				\
      particle_double_real_t rv;				\
      rv = ppsc->kinds[p->kind_].q;				\
      rv;							\
    })

#define particle_double_mni(p) ({				\
      particle_double_real_t rv;				\
      rv = ppsc->kinds[p->kind_].m;				\
      rv;							\
    })

#define particle_double_wni(p) ({				\
      particle_double_real_t rv;				\
      rv = p->qni_wni / ppsc->kinds[p->kind_].q;			\
      rv;							\
    })

static inline particle_double_real_t
particle_double_qni_wni(particle_double_t *p)
{
  return p->qni_wni;
}

static inline int
particle_double_kind(particle_double_t *prt)
{
  return prt->kind_;
}

#endif
