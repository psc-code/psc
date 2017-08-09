
#include <mrc_bits.h>

#include <stdlib.h>

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
#define psc_mparticles_PTYPE_get_n_prts psc_mparticles_single_get_n_prts
#define psc_mparticles_PTYPE_patch_reserve psc_mparticles_single_patch_reserve
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_single_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_single_patch_resize
#define psc_particle_PTYPE_iter_t psc_particle_single_iter_t
#define psc_particle_PTYPE_iter_equal psc_particle_single_iter_equal
#define psc_particle_PTYPE_iter_next psc_particle_single_iter_next 
#define psc_particle_PTYPE_iter_deref psc_particle_single_iter_deref 
#define psc_particle_PTYPE_iter_at psc_particle_single_iter_at 
#define psc_particle_PTYPE_range_t psc_particle_single_range_t
#define psc_particle_PTYPE_range_mprts psc_particle_single_range_mprts
#define psc_particle_PTYPE_range_size psc_particle_single_range_size

#elif PTYPE == PTYPE_DOUBLE

#define particle_PTYPE_real_t particle_double_real_t
#define particle_PTYPE_t particle_double_t
#define psc_particle_PTYPE psc_particle_double
#define psc_mparticles_PTYPE_patch psc_mparticles_double_patch
#define psc_mparticles_PTYPE psc_mparticles_double
#define psc_mparticles_PTYPE_ops psc_mparticles_double_ops
#define psc_mparticles_PTYPE_get_one psc_mparticles_double_get_one
#define psc_mparticles_PTYPE_get_n_prts psc_mparticles_double_get_n_prts
#define psc_mparticles_PTYPE_patch_reserve psc_mparticles_double_patch_reserve
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_double_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_double_patch_resize
#define psc_particle_PTYPE_iter_t psc_particle_double_iter_t
#define psc_particle_PTYPE_iter_equal psc_particle_double_iter_equal
#define psc_particle_PTYPE_iter_next psc_particle_double_iter_next 
#define psc_particle_PTYPE_iter_deref psc_particle_double_iter_deref 
#define psc_particle_PTYPE_iter_at psc_particle_double_iter_at 
#define psc_particle_PTYPE_range_t psc_particle_double_range_t 
#define psc_particle_PTYPE_range_mprts psc_particle_double_range_mprts 
#define psc_particle_PTYPE_range_size psc_particle_double_range_size 

#elif PTYPE == PTYPE_SINGLE_BY_BLOCK

#define particle_PTYPE_real_t particle_single_by_block_real_t
#define particle_PTYPE_t particle_single_by_block_t
#define psc_particle_PTYPE psc_particle_PTYPE_single_by_block
#define psc_mparticles_PTYPE_patch psc_mparticles_single_by_block_patch
#define psc_mparticles_PTYPE psc_mparticles_single_by_block
#define psc_mparticles_PTYPE_ops psc_mparticles_single_by_block_ops
#define psc_mparticles_PTYPE_get_one psc_mparticles_single_by_block_get_one
#define psc_mparticles_PTYPE_get_n_prts psc_mparticles_single_by_block_get_n_prts
#define psc_mparticles_PTYPE_patch_reserve psc_mparticles_single_by_block_patch_reserve
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_single_by_block_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_single_by_block_patch_resize
#define psc_particle_PTYPE_iter_t psc_particle_single_by_block_iter_t
#define psc_particle_PTYPE_iter_equal psc_particle_single_by_block_iter_equal
#define psc_particle_PTYPE_iter_next psc_particle_single_by_block_iter_next 
#define psc_particle_PTYPE_iter_deref psc_particle_single_by_block_iter_deref 
#define psc_particle_PTYPE_iter_at psc_particle_single_by_block_iter_at 
#define psc_particle_PTYPE_range_t psc_particle_single_by_block_range_t 
#define psc_particle_PTYPE_range_mprts psc_particle_single_by_block_range_mprts 
#define psc_particle_PTYPE_range_size psc_particle_single_by_block_range_size 

#elif PTYPE == PTYPE_C

#define particle_PTYPE_real_t particle_c_real_t
#define particle_PTYPE_t particle_c_t
#define psc_particle_PTYPE psc_particle_PTYPE_c
#define psc_mparticles_PTYPE_patch psc_mparticles_c_patch
#define psc_mparticles_PTYPE psc_mparticles_c
#define psc_mparticles_PTYPE_ops psc_mparticles_c_ops
#define psc_mparticles_PTYPE_get_one psc_mparticles_c_get_one
#define psc_mparticles_PTYPE_get_n_prts psc_mparticles_c_get_n_prts
#define psc_mparticles_PTYPE_patch_reserve psc_mparticles_c_patch_reserve
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_c_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_c_patch_resize
#define psc_particle_PTYPE_iter_t psc_particle_c_iter_t
#define psc_particle_PTYPE_iter_equal psc_particle_c_iter_equal
#define psc_particle_PTYPE_iter_next psc_particle_c_iter_next 
#define psc_particle_PTYPE_iter_deref psc_particle_c_iter_deref 
#define psc_particle_PTYPE_iter_at psc_particle_c_iter_at 
#define psc_particle_PTYPE_range_t psc_particle_c_range_t 
#define psc_particle_PTYPE_range_mprts psc_particle_c_range_mprts 
#define psc_particle_PTYPE_range_size psc_particle_c_range_size 

#elif PTYPE == PTYPE_FORTRAN

#define particle_PTYPE_real_t particle_fortran_real_t
#define particle_PTYPE_t particle_fortran_t
#define psc_particle_PTYPE psc_particle_PTYPE_fortran
#define psc_mparticles_PTYPE_patch psc_mparticles_fortran_patch
#define psc_mparticles_PTYPE psc_mparticles_fortran
#define psc_mparticles_PTYPE_ops psc_mparticles_fortran_ops
#define psc_mparticles_PTYPE_get_one psc_mparticles_fortran_get_one
#define psc_mparticles_PTYPE_get_n_prts psc_mparticles_fortran_get_n_prts
#define psc_mparticles_PTYPE_patch_reserve psc_mparticles_fortran_patch_reserve
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_fortran_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_fortran_patch_resize
#define psc_particle_PTYPE_iter_t psc_particle_fortran_iter_t
#define psc_particle_PTYPE_iter_equal psc_particle_fortran_iter_equal
#define psc_particle_PTYPE_iter_next psc_particle_fortran_iter_next 
#define psc_particle_PTYPE_iter_deref psc_particle_fortran_iter_deref 
#define psc_particle_PTYPE_iter_at psc_particle_fortran_iter_at 
#define psc_particle_PTYPE_range_t psc_particle_fortran_range_t 
#define psc_particle_PTYPE_range_mprts psc_particle_fortran_range_mprts 
#define psc_particle_PTYPE_range_size psc_particle_fortran_range_size 

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
// psc_mparticles_PTYPE_get_n_prts

static inline int
psc_mparticles_PTYPE_get_n_prts(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);

  return sub->patch[p].n_prts;
}

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch_reserve

static inline void
psc_mparticles_PTYPE_patch_reserve(struct psc_mparticles *mprts, int p, int new_capacity)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);
  struct psc_mparticles_PTYPE_patch *patch = &sub->patch[p];

  if (new_capacity <= patch->n_alloced)
    return;

  int n_alloced = MAX(new_capacity, patch->n_alloced * 2);
  patch->n_alloced = n_alloced;

  patch->prt_array = realloc(patch->prt_array, n_alloced * sizeof(*patch->prt_array));

#if PTYPE == PTYPE_SINGLE
  free(patch->prt_array_alt);
  patch->prt_array_alt = malloc(n_alloced * sizeof(*patch->prt_array_alt));
  patch->b_idx = realloc(patch->b_idx, n_alloced * sizeof(*patch->b_idx));
  patch->b_ids = realloc(patch->b_ids, n_alloced * sizeof(*patch->b_ids));
#endif

#if PTYPE == PTYPE_SINGLE_BY_BLOCK
  free(patch->prt_array_alt);
  patch->prt_array_alt = malloc(n_alloced * sizeof(*patch->prt_array_alt));
  patch->b_idx = realloc(patch->b_idx, n_alloced * sizeof(*patch->b_idx));
  patch->b_ids = realloc(patch->b_ids, n_alloced * sizeof(*patch->b_ids));
#endif
}

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch_resize

static inline void
psc_mparticles_PTYPE_patch_resize(struct psc_mparticles *mprts, int p, int n_prts)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);
  struct psc_mparticles_PTYPE_patch *patch = &sub->patch[p];

  assert(n_prts <= patch->n_alloced);
  patch->n_prts = n_prts;
}

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch_push_back

static inline void
psc_mparticles_PTYPE_patch_push_back(struct psc_mparticles *mprts, int p,
			       particle_PTYPE_t prt)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);
  struct psc_mparticles_PTYPE_patch *patch = &sub->patch[p];
  
  int n = patch->n_prts;
  if (n == patch->n_alloced) {
    psc_mparticles_PTYPE_patch_reserve(mprts, p, n + 1);
  }
  patch->prt_array[n++] = prt;
  patch->n_prts = n;
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_iter_t

typedef struct {
  int n;
  int p;
  const struct psc_mparticles *mprts;
} psc_particle_PTYPE_iter_t;

// ----------------------------------------------------------------------
// psc_particle_PTYPE_iter_equal

static inline bool
psc_particle_PTYPE_iter_equal(psc_particle_PTYPE_iter_t iter, psc_particle_PTYPE_iter_t iter2)
{
  assert(iter.mprts == iter2.mprts && iter.p == iter2.p);
  return iter.n == iter2.n;
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_iter_next

static inline psc_particle_PTYPE_iter_t
psc_particle_PTYPE_iter_next(psc_particle_PTYPE_iter_t iter)
{
  return (psc_particle_PTYPE_iter_t) {
    .n    = iter.n + 1,
    .p    = iter.p,
    .mprts = iter.mprts,
  };
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_iter_deref

static inline particle_PTYPE_t *
psc_particle_PTYPE_iter_deref(psc_particle_PTYPE_iter_t iter)
{
  // FIXME, shouldn't have to cast away const
  return psc_mparticles_PTYPE_get_one((struct psc_mparticles *) iter.mprts, iter.p, iter.n);
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_iter_at

static inline particle_PTYPE_t *
psc_particle_PTYPE_iter_at(psc_particle_PTYPE_iter_t iter, int m)
{
  // FIXME, shouldn't have to cast away const
  return psc_mparticles_PTYPE_get_one((struct psc_mparticles *) iter.mprts, iter.p, iter.n + m);
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_range_t

typedef struct {
  psc_particle_PTYPE_iter_t begin;
  psc_particle_PTYPE_iter_t end;
} psc_particle_PTYPE_range_t;

// ----------------------------------------------------------------------
// psc_particle_PTYPE_range_mprts

static inline psc_particle_PTYPE_range_t
psc_particle_PTYPE_range_mprts(struct psc_mparticles *mprts, int p)
{
  return (psc_particle_PTYPE_range_t) {
    .begin = { .n = 0                                        , .p = p, .mprts = mprts },
    .end   = { .n = psc_mparticles_PTYPE_get_n_prts(mprts, p), .p = p, .mprts = mprts },
  };
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_range_size

static inline unsigned int
psc_particle_PTYPE_range_size(psc_particle_PTYPE_range_t prts)
{
  return prts.end.n - prts.begin.n;
}

#include <math.h>

#undef particle_PTYPE_real_t
#undef particle_PTYPE_t
#undef psc_particle_PTYPE
#undef psc_mparticles_PTYPE_patch
#undef psc_mparticles_PTYPE
#undef psc_mparticles_PTYPE_ops
#undef psc_mparticles_PTYPE_get_one
#undef psc_mparticles_PTYPE_get_n_prts
#undef psc_mparticles_PTYPE_patch_reserve
#undef psc_mparticles_PTYPE_patch_push_back
#undef psc_mparticles_PTYPE_patch_resize
#undef psc_particle_PTYPE_iter_t
#undef psc_particle_PTYPE_iter_equal
#undef psc_particle_PTYPE_iter_next 
#undef psc_particle_PTYPE_iter_deref 
#undef psc_particle_PTYPE_iter_at 
#undef psc_particle_PTYPE_range_t 
#undef psc_particle_PTYPE_range_mprts 
#undef psc_particle_PTYPE_range_size 
