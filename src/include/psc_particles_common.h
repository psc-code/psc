
#include <mrc_bits.h>

#include <stdlib.h>

#include "psc_particle_buf_common.h"

#if PTYPE == PTYPE_SINGLE

#define particle_PTYPE_real_t particle_single_real_t
#define particle_PTYPE_t particle_single_t

#define psc_particle_PTYPE_buf_t psc_particle_single_buf_t
#define psc_particle_PTYPE_buf_ctor psc_particle_single_buf_ctor
#define psc_particle_PTYPE_buf_dtor psc_particle_single_buf_dtor
#define psc_particle_PTYPE_buf_size psc_particle_single_buf_size
#define psc_particle_PTYPE_buf_resize psc_particle_single_buf_resize
#define psc_particle_PTYPE_buf_reserve psc_particle_single_buf_reserve
#define psc_particle_PTYPE_buf_capacity psc_particle_single_buf_capacity
#define psc_particle_PTYPE_buf_push_back psc_particle_single_buf_push_back
#define psc_particle_PTYPE_buf_at_ptr psc_particle_single_buf_at_ptr

#define psc_mparticles_PTYPE_patch psc_mparticles_single_patch
#define psc_mparticles_PTYPE psc_mparticles_single
#define psc_mparticles_PTYPE_ops psc_mparticles_single_ops
#define psc_mparticles_PTYPE_patch_get_buf psc_mparticles_single_patch_get_buf
#define psc_mparticles_PTYPE_get_one psc_mparticles_single_get_one
#define psc_mparticles_PTYPE_get_n_prts psc_mparticles_single_get_n_prts
#define psc_mparticles_PTYPE_patch_reserve psc_mparticles_single_patch_reserve
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_single_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_single_patch_resize
#define psc_mparticles_PTYPE_patch_capacity psc_mparticles_single_patch_capacity
#define psc_mparticles_PTYPE_patch_get_b_dxi psc_mparticles_single_patch_get_b_dxi 
#define psc_mparticles_PTYPE_patch_get_b_mx psc_mparticles_single_patch_get_b_mx
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

#define psc_particle_PTYPE_buf_t psc_particle_double_buf_t
#define psc_particle_PTYPE_buf_ctor psc_particle_double_buf_ctor
#define psc_particle_PTYPE_buf_dtor psc_particle_double_buf_dtor
#define psc_particle_PTYPE_buf_size psc_particle_double_buf_size
#define psc_particle_PTYPE_buf_resize psc_particle_double_buf_resize
#define psc_particle_PTYPE_buf_reserve psc_particle_double_buf_reserve
#define psc_particle_PTYPE_buf_capacity psc_particle_double_buf_capacity
#define psc_particle_PTYPE_buf_push_back psc_particle_double_buf_push_back
#define psc_particle_PTYPE_buf_at_ptr psc_particle_double_buf_at_ptr

#define psc_mparticles_PTYPE_patch psc_mparticles_double_patch
#define psc_mparticles_PTYPE psc_mparticles_double
#define psc_mparticles_PTYPE_ops psc_mparticles_double_ops
#define psc_mparticles_PTYPE_patch_get_buf psc_mparticles_double_patch_get_buf
#define psc_mparticles_PTYPE_get_one psc_mparticles_double_get_one
#define psc_mparticles_PTYPE_get_n_prts psc_mparticles_double_get_n_prts
#define psc_mparticles_PTYPE_patch_reserve psc_mparticles_double_patch_reserve
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_double_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_double_patch_resize
#define psc_mparticles_PTYPE_patch_capacity psc_mparticles_double_patch_capacity
#define psc_mparticles_PTYPE_patch_get_b_dxi psc_mparticles_double_patch_get_b_dxi 
#define psc_mparticles_PTYPE_patch_get_b_mx psc_mparticles_double_patch_get_b_mx
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

#define psc_particle_PTYPE_buf_t psc_particle_single_by_block_buf_t
#define psc_particle_PTYPE_buf_ctor psc_particle_single_by_block_buf_ctor
#define psc_particle_PTYPE_buf_dtor psc_particle_single_by_block_buf_dtor
#define psc_particle_PTYPE_buf_size psc_particle_single_by_block_buf_size
#define psc_particle_PTYPE_buf_resize psc_particle_single_by_block_buf_resize
#define psc_particle_PTYPE_buf_reserve psc_particle_single_by_block_buf_reserve
#define psc_particle_PTYPE_buf_capacity psc_particle_single_by_block_buf_capacity
#define psc_particle_PTYPE_buf_push_back psc_particle_single_by_block_buf_push_back
#define psc_particle_PTYPE_buf_at_ptr psc_particle_single_by_block_buf_at_ptr

#define psc_mparticles_PTYPE_patch psc_mparticles_single_by_block_patch
#define psc_mparticles_PTYPE psc_mparticles_single_by_block
#define psc_mparticles_PTYPE_ops psc_mparticles_single_by_block_ops
#define psc_mparticles_PTYPE_patch_get_buf psc_mparticles_single_by_block_patch_get_buf
#define psc_mparticles_PTYPE_get_one psc_mparticles_single_by_block_get_one
#define psc_mparticles_PTYPE_get_n_prts psc_mparticles_single_by_block_get_n_prts
#define psc_mparticles_PTYPE_patch_reserve psc_mparticles_single_by_block_patch_reserve
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_single_by_block_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_single_by_block_patch_resize
#define psc_mparticles_PTYPE_patch_capacity psc_mparticles_single_by_block_patch_capacity
#define psc_mparticles_PTYPE_patch_get_b_dxi psc_mparticles_single_by_block_patch_get_b_dxi 
#define psc_mparticles_PTYPE_patch_get_b_mx psc_mparticles_single_by_block_patch_get_b_mx
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

#define psc_particle_PTYPE_buf_t psc_particle_c_buf_t
#define psc_particle_PTYPE_buf_ctor psc_particle_c_buf_ctor
#define psc_particle_PTYPE_buf_dtor psc_particle_c_buf_dtor
#define psc_particle_PTYPE_buf_size psc_particle_c_buf_size
#define psc_particle_PTYPE_buf_resize psc_particle_c_buf_resize
#define psc_particle_PTYPE_buf_reserve psc_particle_c_buf_reserve
#define psc_particle_PTYPE_buf_capacity psc_particle_c_buf_capacity
#define psc_particle_PTYPE_buf_push_back psc_particle_c_buf_push_back
#define psc_particle_PTYPE_buf_at_ptr psc_particle_c_buf_at_ptr

#define psc_mparticles_PTYPE_patch psc_mparticles_c_patch
#define psc_mparticles_PTYPE psc_mparticles_c
#define psc_mparticles_PTYPE_ops psc_mparticles_c_ops
#define psc_mparticles_PTYPE_patch_get_buf psc_mparticles_c_patch_get_buf
#define psc_mparticles_PTYPE_get_one psc_mparticles_c_get_one
#define psc_mparticles_PTYPE_get_n_prts psc_mparticles_c_get_n_prts
#define psc_mparticles_PTYPE_patch_reserve psc_mparticles_c_patch_reserve
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_c_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_c_patch_resize
#define psc_mparticles_PTYPE_patch_capacity psc_mparticles_c_patch_capacity
#define psc_mparticles_PTYPE_patch_get_b_dxi psc_mparticles_c_patch_get_b_dxi 
#define psc_mparticles_PTYPE_patch_get_b_mx psc_mparticles_c_patch_get_b_mx
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

#define psc_particle_PTYPE_buf_t psc_particle_fortran_buf_t
#define psc_particle_PTYPE_buf_ctor psc_particle_fortran_buf_ctor
#define psc_particle_PTYPE_buf_dtor psc_particle_fortran_buf_dtor
#define psc_particle_PTYPE_buf_size psc_particle_fortran_buf_size
#define psc_particle_PTYPE_buf_resize psc_particle_fortran_buf_resize
#define psc_particle_PTYPE_buf_reserve psc_particle_fortran_buf_reserve
#define psc_particle_PTYPE_buf_capacity psc_particle_fortran_buf_capacity
#define psc_particle_PTYPE_buf_push_back psc_particle_fortran_buf_push_back
#define psc_particle_PTYPE_buf_at_ptr psc_particle_fortran_buf_at_ptr

#define psc_mparticles_PTYPE_patch psc_mparticles_fortran_patch
#define psc_mparticles_PTYPE psc_mparticles_fortran
#define psc_mparticles_PTYPE_ops psc_mparticles_fortran_ops
#define psc_mparticles_PTYPE_patch_get_buf psc_mparticles_fortran_patch_get_buf
#define psc_mparticles_PTYPE_get_one psc_mparticles_fortran_get_one
#define psc_mparticles_PTYPE_get_n_prts psc_mparticles_fortran_get_n_prts
#define psc_mparticles_PTYPE_patch_reserve psc_mparticles_fortran_patch_reserve
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_fortran_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_fortran_patch_resize
#define psc_mparticles_PTYPE_patch_capacity psc_mparticles_fortran_patch_capacity
#define psc_mparticles_PTYPE_patch_get_b_dxi psc_mparticles_fortran_patch_get_b_dxi 
#define psc_mparticles_PTYPE_patch_get_b_mx psc_mparticles_fortran_patch_get_b_mx
#define psc_particle_PTYPE_iter_t psc_particle_fortran_iter_t
#define psc_particle_PTYPE_iter_equal psc_particle_fortran_iter_equal
#define psc_particle_PTYPE_iter_next psc_particle_fortran_iter_next 
#define psc_particle_PTYPE_iter_deref psc_particle_fortran_iter_deref 
#define psc_particle_PTYPE_iter_at psc_particle_fortran_iter_at 
#define psc_particle_PTYPE_range_t psc_particle_fortran_range_t 
#define psc_particle_PTYPE_range_mprts psc_particle_fortran_range_mprts 
#define psc_particle_PTYPE_range_size psc_particle_fortran_range_size

#elif PTYPE == PTYPE_CUDA

#define particle_PTYPE_real_t particle_cuda_real_t
#define particle_PTYPE_t particle_cuda_t

#define psc_particle_PTYPE_buf_t psc_particle_cuda_buf_t
#define psc_particle_PTYPE_buf_ctor psc_particle_cuda_buf_ctor
#define psc_particle_PTYPE_buf_dtor psc_particle_cuda_buf_dtor
#define psc_particle_PTYPE_buf_size psc_particle_cuda_buf_size
#define psc_particle_PTYPE_buf_resize psc_particle_cuda_buf_resize
#define psc_particle_PTYPE_buf_reserve psc_particle_cuda_buf_reserve
#define psc_particle_PTYPE_buf_capacity psc_particle_cuda_buf_capacity
#define psc_particle_PTYPE_buf_push_back psc_particle_cuda_buf_push_back
#define psc_particle_PTYPE_buf_at_ptr psc_particle_cuda_buf_at_ptr

#define psc_mparticles_PTYPE psc_mparticles_cuda
#define psc_mparticles_PTYPE_patch_get_b_dxi psc_mparticles_cuda_patch_get_b_dxi 
#define psc_mparticles_PTYPE_patch_get_b_mx psc_mparticles_cuda_patch_get_b_mx

#endif

// ======================================================================

#if PTYPE == PTYPE_CUDA

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE

struct psc_mparticles_PTYPE {
  struct cuda_mparticles *cmprts;
};

const particle_PTYPE_real_t *psc_mparticles_PTYPE_patch_get_b_dxi(struct psc_mparticles *mprts, int p);
const int *psc_mparticles_PTYPE_patch_get_b_mx(struct psc_mparticles *mprts, int p);

#else // PTYPE != PTYPE_CUDA

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch

struct psc_mparticles_PTYPE_patch {
  psc_particle_PTYPE_buf_t buf;

  int b_mx[3];
  particle_PTYPE_real_t b_dxi[3];
#if PTYPE == PTYPE_SINGLE
  particle_PTYPE_t *prt_array_alt;
  int nr_blocks;
  unsigned int *b_idx;
  unsigned int *b_ids;
  unsigned int *b_cnt;
  bool need_reorder;
#endif
  
#if PTYPE == PTYPE_SINGLE_BY_BLOCK
  particle_PTYPE_t *prt_array_alt;
  int nr_blocks;
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
// psc_mparticles_PTYPE_patch_get_buf

static inline psc_particle_PTYPE_buf_t *
psc_mparticles_PTYPE_patch_get_buf(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);
  struct psc_mparticles_PTYPE_patch *patch = &sub->patch[p];

  return &patch->buf;
}

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_get_one

static inline particle_PTYPE_t *
psc_mparticles_PTYPE_get_one(struct psc_mparticles *mprts, int p, unsigned int n)
{
  assert(psc_mparticles_ops(mprts) == &psc_mparticles_PTYPE_ops);
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);
  struct psc_mparticles_PTYPE_patch *patch = &sub->patch[p];

  return psc_particle_PTYPE_buf_at_ptr(&patch->buf, n);
}

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_get_n_prts

static inline int
psc_mparticles_PTYPE_get_n_prts(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);

  return psc_particle_PTYPE_buf_size(&sub->patch[p].buf);
}

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch_reserve

static inline void
psc_mparticles_PTYPE_patch_reserve(struct psc_mparticles *mprts, int p,
				   unsigned int new_capacity)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);
  struct psc_mparticles_PTYPE_patch *patch = &sub->patch[p];

  unsigned int old_capacity = psc_particle_PTYPE_buf_capacity(&patch->buf);
  psc_particle_PTYPE_buf_reserve(&patch->buf, new_capacity);
  new_capacity = psc_particle_PTYPE_buf_capacity(&patch->buf);

  if (new_capacity == old_capacity) {
    return;
  }

#if PTYPE == PTYPE_SINGLE
  free(patch->prt_array_alt);
  patch->prt_array_alt = (particle_PTYPE_t *) malloc(new_capacity * sizeof(*patch->prt_array_alt));
  patch->b_idx = (unsigned int *) realloc(patch->b_idx, new_capacity * sizeof(*patch->b_idx));
  patch->b_ids = (unsigned int *) realloc(patch->b_ids, new_capacity * sizeof(*patch->b_ids));
#endif

#if PTYPE == PTYPE_SINGLE_BY_BLOCK
  free(patch->prt_array_alt);
  patch->prt_array_alt = malloc(new_capacity * sizeof(*patch->prt_array_alt));
  patch->b_idx = realloc(patch->b_idx, new_capacity * sizeof(*patch->b_idx));
  patch->b_ids = realloc(patch->b_ids, new_capacity * sizeof(*patch->b_ids));
#endif
}

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch_resize

static inline void
psc_mparticles_PTYPE_patch_resize(struct psc_mparticles *mprts, int p, int n_prts)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);
  struct psc_mparticles_PTYPE_patch *patch = &sub->patch[p];

  psc_particle_PTYPE_buf_resize(&patch->buf, n_prts);
}

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch_capacity

static inline unsigned int
psc_mparticles_PTYPE_patch_capacity(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);
  struct psc_mparticles_PTYPE_patch *patch = &sub->patch[p];

  return psc_particle_PTYPE_buf_capacity(&patch->buf);
}

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch_push_back

static inline void
psc_mparticles_PTYPE_patch_push_back(struct psc_mparticles *mprts, int p,
				     particle_PTYPE_t prt)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);
  struct psc_mparticles_PTYPE_patch *patch = &sub->patch[p];

  psc_particle_PTYPE_buf_push_back(&patch->buf, prt);
}

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch_get_b_mx

static inline const int *
psc_mparticles_PTYPE_patch_get_b_mx(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);
  struct psc_mparticles_PTYPE_patch *patch = &sub->patch[p];

  return patch->b_mx;
}

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch_get_b_dxi

static inline const particle_PTYPE_real_t *
psc_mparticles_PTYPE_patch_get_b_dxi(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);
  struct psc_mparticles_PTYPE_patch *patch = &sub->patch[p];

  return patch->b_dxi;
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
  psc_particle_PTYPE_iter_t rv;
  rv.n     = iter.n + 1;
  rv.p     = iter.p;
  rv.mprts = iter.mprts;

  return rv;
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
  psc_particle_PTYPE_range_t rv;
  rv.begin.n     = 0;
  rv.begin.p     = p;
  rv.begin.mprts = mprts;
  rv.end.n       = psc_mparticles_PTYPE_get_n_prts(mprts, p);
  rv.end.p       = p;
  rv.end.mprts   = mprts;

  return rv;
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_range_size

static inline unsigned int
psc_particle_PTYPE_range_size(psc_particle_PTYPE_range_t prts)
{
  return prts.end.n - prts.begin.n;
}

#include <math.h>

#endif // PTYPE_CUDA

#undef particle_PTYPE_real_t
#undef particle_PTYPE_t

#undef psc_particle_PTYPE_buf_t
#undef psc_particle_PTYPE_buf_ctor
#undef psc_particle_PTYPE_buf_dtor
#undef psc_particle_PTYPE_buf_size
#undef psc_particle_PTYPE_buf_resize
#undef psc_particle_PTYPE_buf_reserve
#undef psc_particle_PTYPE_buf_capacity
#undef psc_particle_PTYPE_buf_push_back
#undef psc_particle_PTYPE_buf_at_ptr

#undef psc_mparticles_PTYPE_patch
#undef psc_mparticles_PTYPE
#undef psc_mparticles_PTYPE_ops
#undef psc_mparticles_PTYPE_patch_get_buf
#undef psc_mparticles_PTYPE_get_one
#undef psc_mparticles_PTYPE_get_n_prts
#undef psc_mparticles_PTYPE_patch_reserve
#undef psc_mparticles_PTYPE_patch_push_back
#undef psc_mparticles_PTYPE_patch_resize
#undef psc_mparticles_PTYPE_patch_capacity
#undef psc_mparticles_PTYPE_patch_get_b_dxi
#undef psc_mparticles_PTYPE_patch_get_b_mx
#undef psc_particle_PTYPE_iter_t
#undef psc_particle_PTYPE_iter_equal
#undef psc_particle_PTYPE_iter_next 
#undef psc_particle_PTYPE_iter_deref 
#undef psc_particle_PTYPE_iter_at 
#undef psc_particle_PTYPE_range_t 
#undef psc_particle_PTYPE_range_mprts 
#undef psc_particle_PTYPE_range_size

