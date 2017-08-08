
#define PTYPE_SINGLE          1
#define PTYPE_DOUBLE          2
#define PTYPE_SINGLE_BY_BLOCK 3
#define PTYPE_C               4
#define PTYPE_FORTRAN         5

#if PTYPE == PTYPE_SINGLE

#define particle_PTYPE_real_t particle_single_real_t
#define particle_PTYPE_t particle_single_t
#define psc_particle_PTYPE psc_particle_PTYPE_single
#define psc_mparticles_PTYPE_patch psc_mparticles_single_patch
#define psc_mparticles_PTYPE psc_mparticles_single
#define psc_mparticles_PTYPE_ops psc_mparticles_single_ops
#define psc_mparticles_PTYPE_get_one psc_mparticles_single_get_one

#elif PTYPE == PTYPE_DOUBLE

#define particle_PTYPE_real_t particle_double_real_t
#define particle_PTYPE_t particle_double_t
#define psc_particle_PTYPE psc_particle_double
#define psc_mparticles_PTYPE_patch psc_mparticles_double_patch
#define psc_mparticles_PTYPE psc_mparticles_double
#define psc_mparticles_PTYPE_ops psc_mparticles_double_ops
#define psc_mparticles_PTYPE_get_one psc_mparticles_double_get_one

#elif PTYPE == PTYPE_SINGLE_BY_BLOCK

#define particle_PTYPE_real_t particle_single_by_block_real_t
#define particle_PTYPE_t particle_single_by_block_t
#define psc_particle_PTYPE psc_particle_PTYPE_single_by_block
#define psc_mparticles_PTYPE_patch psc_mparticles_single_by_block_patch
#define psc_mparticles_PTYPE psc_mparticles_single_by_block
#define psc_mparticles_PTYPE_ops psc_mparticles_single_by_block_ops
#define psc_mparticles_PTYPE_get_one psc_mparticles_single_by_block_get_one

#elif PTYPE == PTYPE_C

#define particle_PTYPE_real_t particle_c_real_t
#define particle_PTYPE_t particle_c_t
#define psc_particle_PTYPE psc_particle_PTYPE_c
#define psc_mparticles_PTYPE_patch psc_mparticles_c_patch
#define psc_mparticles_PTYPE psc_mparticles_c
#define psc_mparticles_PTYPE_ops psc_mparticles_c_ops
#define psc_mparticles_PTYPE_get_one psc_mparticles_c_get_one

#elif PTYPE == PTYPE_FORTRAN

#define particle_PTYPE_real_t particle_fortran_real_t
#define particle_PTYPE_t particle_fortran_t
#define psc_particle_PTYPE psc_particle_PTYPE_fortran
#define psc_mparticles_PTYPE_patch psc_mparticles_fortran_patch
#define psc_mparticles_PTYPE psc_mparticles_fortran
#define psc_mparticles_PTYPE_ops psc_mparticles_fortran_ops
#define psc_mparticles_PTYPE_get_one psc_mparticles_fortran_get_one

#endif

// ----------------------------------------------------------------------
// particle_PYTPE_real_t

#if PTYPE == PTYPE_SINGLE || PTYPE == PTYPE_SINGLE_BY_BLOCK

typedef float particle_PTYPE_real_t;

#elif PTYPE == PTYPE_DOUBLE || PTYPE == PTYPE_C || PTYPE == PTYPE_FORTRAN

typedef double particle_PTYPE_real_t;

#endif

// ----------------------------------------------------------------------
// MPI_PARTICLES_PTYPE_REAL
// annoying, but need to use a macro, which means we can't consolidate float/double

#if PTYPE == PTYPE_SINGLE

#define MPI_PARTICLES_SINGLE_REAL MPI_FLOAT
#define psc_mparticles_single(mprts) mrc_to_subobj(mprts, struct psc_mparticles_single)

#elif PTYPE == PTYPE_DOUBLE

#define MPI_PARTICLES_DOUBLE_REAL MPI_DOUBLE
#define psc_mparticles_double(prts) mrc_to_subobj(prts, struct psc_mparticles_double)

#elif PTYPE == PTYPE_SINGLE_BY_BLOCK

#define MPI_PARTICLES_SINGLE_REAL MPI_FLOAT
#define psc_mparticles_single_by_block(prts) mrc_to_subobj(prts, struct psc_mparticles_single_by_block)

#elif PTYPE == PTYPE_C

#define MPI_PARTICLES_C_REAL MPI_DOUBLE
#define psc_mparticles_c(prts) mrc_to_subobj(prts, struct psc_mparticles_c)

#elif PTYPE == PTYPE_FORTRAN

#define MPI_PARTICLES_FORTRAN_REAL MPI_DOUBLE
#define psc_mparticles_fortran(prts) mrc_to_subobj(prts, struct psc_mparticles_fortran)

#endif

// ----------------------------------------------------------------------
// particle_PTYPE_t

#if PTYPE == PTYPE_SINGLE || PTYPE == PTYPE_SINGLE_BY_BLOCK || PTYPE == PTYPE_DOUBLE

typedef struct psc_particle_PTYPE {
  particle_PTYPE_real_t xi, yi, zi;
  particle_PTYPE_real_t qni_wni;
  particle_PTYPE_real_t pxi, pyi, pzi;
  int kind;
} particle_PTYPE_t;

#elif PTYPE == PTYPE_C

typedef struct psc_particle_PTYPE {
  particle_PTYPE_real_t xi, yi, zi;
  particle_PTYPE_real_t pxi, pyi, pzi;
  particle_PTYPE_real_t qni;
  particle_PTYPE_real_t mni;
  particle_PTYPE_real_t wni;
  long long kind; // 64 bits to match the other members, for bnd exchange
} particle_PTYPE_t;

#elif PTYPE == PTYPE_FORTRAN

typedef struct psc_particle_PTYPE {
  particle_PTYPE_real_t xi, yi, zi;
  particle_PTYPE_real_t pxi, pyi, pzi;
  particle_PTYPE_real_t qni;
  particle_PTYPE_real_t mni;
  particle_PTYPE_real_t cni;
  particle_PTYPE_real_t lni;
  particle_PTYPE_real_t wni;
} particle_PTYPE_t;

#endif

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch

struct psc_mparticles_PTYPE_patch {
  particle_PTYPE_t *prt_array;
  int n_prts;
  int n_alloced;

#if PTYPE == PTYPE_SINGLE
  particle_PTYPE_t *prt_array_alt;
  int b_mx[3];
  int nr_blocks;
  particle_PTYPE_real_t b_dxi[3];
  unsigned int *b_idx;
  unsigned int *b_ids;
  unsigned int *b_cnt;
  unsigned int n_send;
  unsigned int n_part_save;
  bool need_reorder;
#endif
  
#if PTYPE == PTYPE_SINGLE_BY_BLOCK
  particle_PTYPE_t *prt_array_alt;
  int b_mx[3];
  int nr_blocks;
  particle_PTYPE_real_t b_dxi[3];
  unsigned int *b_idx;
  unsigned int *b_ids;
  unsigned int *b_cnt;
  unsigned int *b_off;
  bool need_reorder;
#endif
};

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE

struct psc_mparticles_PTYPE {
  struct psc_mparticles_PTYPE_patch *patch;
};

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_get_one

static inline particle_PTYPE_t *
psc_mparticles_PTYPE_get_one(struct psc_mparticles *mprts, int p, int n)
{
  assert(psc_mparticles_ops(mprts) == &psc_mparticles_PTYPE_ops);
  return &psc_mparticles_PTYPE(mprts)->patch[p].prt_array[n];
}

// ----------------------------------------------------------------------
//

#if PTYPE == PTYPE_SINGLE

#elif PTYPE == PTYPE_DOUBLE

#elif PTYPE == PTYPE_SINGLE_BY_BLOCK

#elif PTYPE == PTYPE_C

#elif PTYPE == PTYPE_FORTRAN

#endif



#include <math.h>

#undef particle_PTYPE_real_t
#undef particle_PTYPE_t
#undef psc_particle_PTYPE
#undef psc_mparticles_PTYPE_patch
#undef psc_mparticles_PTYPE
#undef psc_mparticles_PTYPE_ops
#undef psc_mparticles_PTYPE_get_one
