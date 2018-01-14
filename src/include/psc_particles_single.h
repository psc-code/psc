
#ifndef PSC_PARTICLE_SINGLE_H
#define PSC_PARTICLE_SINGLE_H

#include "psc_particles_private.h"

#include "particles.hxx"
#include "particles_traits.hxx"

#include <vector>

using particle_single_real_t = float;

struct particle_single_t : psc_particle<particle_single_real_t> {};

#define psc_mparticles_single(mprts) mrc_to_subobj(mprts, struct psc_mparticles_single)

template<>
struct psc_mparticles_patch<particle_single_t> : particles_base<particle_single_t>
{
  using Base = particles_base<particle_single_t>;
  
  ~psc_mparticles_patch()
  {
    free(prt_array_alt);
    free(b_idx);
    free(b_ids);
    free(b_cnt);
  }

  void reserve(unsigned int new_capacity)
  {
    unsigned int old_capacity = Base::buf.capacity();
    Base::reserve(new_capacity);

    new_capacity = Base::buf.capacity();
    if (new_capacity != old_capacity) {
      free(prt_array_alt);
      prt_array_alt = (particle_single_t *) malloc(new_capacity * sizeof(*prt_array_alt));
      b_idx = (unsigned int *) realloc(b_idx, new_capacity * sizeof(*b_idx));
      b_ids = (unsigned int *) realloc(b_ids, new_capacity * sizeof(*b_ids));
    }
  }

  particle_single_t *prt_array_alt;
  int nr_blocks;
  unsigned int *b_idx;
  unsigned int *b_ids;
  unsigned int *b_cnt;
  bool need_reorder;
};

using psc_mparticles_single = psc_mparticles_<particle_single_t>;
using mparticles_single_t = mparticles<psc_mparticles_single>;

template<>
struct mparticles_traits<mparticles_single_t>
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
