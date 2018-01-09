
#ifndef PSC_PARTICLE_SINGLE_H
#define PSC_PARTICLE_SINGLE_H

#include "psc_particles_private.h"
#include "psc.h"

#include "particles_traits.hxx"

#define PTYPE PTYPE_SINGLE
#include "psc_particle_buf_common.h"
#include "psc_particles_common.h"
#undef PTYPE

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
      rv = p->qni_wni / ppsc->kinds[p->kind_].q;			\
      rv;							\
    })

#define particle_single_qni_wni(prt) ((prt)->qni_wni)

#endif
