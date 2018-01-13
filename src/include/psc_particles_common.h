
#include <cstdlib>

#define PTYPE_SINGLE          1
#define PTYPE_DOUBLE          2
#define PTYPE_FORTRAN         5
#define PTYPE_CUDA            6

#if PTYPE == PTYPE_SINGLE

#define particle_PTYPE_real_t particle_single_real_t
#define particle_PTYPE_t particle_single_t

#define psc_mparticles_PTYPE_patch psc_mparticles_single_patch
#define psc_mparticles_PTYPE psc_mparticles_single
#define psc_mparticles_PTYPE_ops psc_mparticles_single_ops
#define psc_particle_PTYPE_iter_t psc_particle_single_iter_t
#define psc_particle_PTYPE_range_t psc_particle_single_range_t

#elif PTYPE == PTYPE_DOUBLE

#define particle_PTYPE_real_t particle_double_real_t
#define particle_PTYPE_t particle_double_t

#define psc_mparticles_PTYPE_patch psc_mparticles_double_patch
#define psc_mparticles_PTYPE psc_mparticles_double
#define psc_mparticles_PTYPE_ops psc_mparticles_double_ops
#define psc_particle_PTYPE_iter_t psc_particle_double_iter_t
#define psc_particle_PTYPE_range_t psc_particle_double_range_t 

#elif PTYPE == PTYPE_FORTRAN

#define particle_PTYPE_real_t particle_fortran_real_t
#define particle_PTYPE_t particle_fortran_t

#define psc_mparticles_PTYPE_patch psc_mparticles_fortran_patch
#define psc_mparticles_PTYPE psc_mparticles_fortran
#define psc_mparticles_PTYPE_ops psc_mparticles_fortran_ops
#define psc_particle_PTYPE_iter_t psc_particle_fortran_iter_t
#define psc_particle_PTYPE_range_t psc_particle_fortran_range_t 

#endif

// ======================================================================

template<typename P>
struct psc_mparticles_PTYPE_patch;

using psc_particle_PTYPE_range_t = psc_mparticles_PTYPE_patch<particle_PTYPE_t>&;

template<typename P>
struct psc_mparticles_PTYPE_patch : particles_base<P>
{
  using Base = particles_base<P>;
  
  ~psc_mparticles_PTYPE_patch()
  {
#if PTYPE == PTYPE_SINGLE
    free(prt_array_alt);
    free(b_idx);
    free(b_ids);
    free(b_cnt);
#endif
  }

  psc_particle_PTYPE_range_t range()
  {
    return psc_particle_PTYPE_range_t(*this);
  }

  void reserve(unsigned int new_capacity)
  {
    unsigned int old_capacity = Base::buf.capacity();
    Base::reserve(new_capacity);

    new_capacity = Base::buf.capacity();
    if (new_capacity != old_capacity) {
#if PTYPE == PTYPE_SINGLE
      free(prt_array_alt);
      prt_array_alt = (particle_PTYPE_t *) malloc(new_capacity * sizeof(*prt_array_alt));
      b_idx = (unsigned int *) realloc(b_idx, new_capacity * sizeof(*b_idx));
      b_ids = (unsigned int *) realloc(b_ids, new_capacity * sizeof(*b_ids));
#endif
    }
  }

#if PTYPE == PTYPE_SINGLE
  particle_PTYPE_t *prt_array_alt;
  int nr_blocks;
  unsigned int *b_idx;
  unsigned int *b_ids;
  unsigned int *b_cnt;
  bool need_reorder;
#endif

};

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE

struct psc_mparticles_PTYPE
{
  using particles_t = psc_mparticles_PTYPE_patch<particle_PTYPE_t>;
  
  particles_t *patch;
};

#include <math.h>

#undef particle_PTYPE_real_t
#undef particle_PTYPE_t

#undef psc_mparticles_PTYPE_patch
#undef psc_mparticles_PTYPE
#undef psc_mparticles_PTYPE_ops
#undef psc_mparticles_PTYPE_patch_capacity
#undef psc_particle_PTYPE_iter_t
#undef psc_particle_PTYPE_range_t 

