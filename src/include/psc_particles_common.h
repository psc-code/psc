

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
#define psc_mparticles_PTYPE_get_one psc_mparticles_single_get_one
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_single_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_single_patch_resize
#define psc_mparticles_PTYPE_patch_capacity psc_mparticles_single_patch_capacity
#define psc_particle_PTYPE_iter_t psc_particle_single_iter_t
#define psc_particle_PTYPE_range_t psc_particle_single_range_t

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
#define psc_mparticles_PTYPE_get_one psc_mparticles_double_get_one
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_double_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_double_patch_resize
#define psc_mparticles_PTYPE_patch_capacity psc_mparticles_double_patch_capacity
#define psc_particle_PTYPE_iter_t psc_particle_double_iter_t
#define psc_particle_PTYPE_range_t psc_particle_double_range_t 

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
#define psc_mparticles_PTYPE_get_one psc_mparticles_single_by_block_get_one
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_single_by_block_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_single_by_block_patch_resize
#define psc_mparticles_PTYPE_patch_capacity psc_mparticles_single_by_block_patch_capacity
#define psc_particle_PTYPE_iter_t psc_particle_single_by_block_iter_t
#define psc_particle_PTYPE_range_t psc_particle_single_by_block_range_t 

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
#define psc_mparticles_PTYPE_get_one psc_mparticles_fortran_get_one
#define psc_mparticles_PTYPE_patch_push_back psc_mparticles_fortran_patch_push_back
#define psc_mparticles_PTYPE_patch_resize psc_mparticles_fortran_patch_resize
#define psc_mparticles_PTYPE_patch_capacity psc_mparticles_fortran_patch_capacity
#define psc_particle_PTYPE_iter_t psc_particle_fortran_iter_t
#define psc_particle_PTYPE_range_t psc_particle_fortran_range_t 

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

#endif

// ======================================================================

#if PTYPE != PTYPE_CUDA

struct psc_particle_PTYPE_range_t;

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch

struct psc_mparticles_PTYPE_patch
{
  psc_particle_PTYPE_buf_t buf;

  int b_mx[3];
  particle_PTYPE_real_t b_dxi[3];

  struct psc_mparticles *mprts;
  int p;
  
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

  void reserve(unsigned int new_capacity);

  psc_particle_PTYPE_range_t range();

  unsigned int size() const
  {
    return psc_particle_PTYPE_buf_size(&buf);
  }

  void push_back(const particle_PTYPE_t& prt)
  {
    psc_particle_PTYPE_buf_push_back(&buf, prt);
  }

  psc_particle_PTYPE_buf_t& get_buf()
  {
    return buf;
  }

  const int* get_b_mx() const
  {
    return b_mx;
  }

  const particle_PTYPE_real_t* get_b_dxi() const
  {
    return b_dxi;
  }
};

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE

struct psc_mparticles_PTYPE
{
  using patch_t = psc_mparticles_PTYPE_patch;
  
  patch_t *patch;
};

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
// psc_mparticles_PTYPE_patch::reserve

inline void psc_mparticles_PTYPE_patch::reserve(unsigned int new_capacity)
{
  unsigned int old_capacity = psc_particle_PTYPE_buf_capacity(&buf);
  psc_particle_PTYPE_buf_reserve(&buf, new_capacity);
  new_capacity = psc_particle_PTYPE_buf_capacity(&buf);

  if (new_capacity == old_capacity) {
    return;
  }

#if PTYPE == PTYPE_SINGLE
  free(prt_array_alt);
  prt_array_alt = (particle_PTYPE_t *) malloc(new_capacity * sizeof(*prt_array_alt));
  b_idx = (unsigned int *) realloc(b_idx, new_capacity * sizeof(*b_idx));
  b_ids = (unsigned int *) realloc(b_ids, new_capacity * sizeof(*b_ids));
#endif

#if PTYPE == PTYPE_SINGLE_BY_BLOCK
  free(prt_array_alt);
  prt_array_alt = (particle_PTYPE_t *) malloc(new_capacity * sizeof(*prt_array_alt));
  b_idx = (unsigned int *) realloc(b_idx, new_capacity * sizeof(*b_idx));
  b_ids = (unsigned int *) realloc(b_ids, new_capacity * sizeof(*b_ids));
#endif
}

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch_resize

static inline void
psc_mparticles_PTYPE_patch_resize(struct psc_mparticles *mprts, int p, int n_prts)
{
  struct psc_mparticles_PTYPE *sub = psc_mparticles_PTYPE(mprts);
  struct psc_mparticles_PTYPE_patch *patch = &sub->patch[p];

  patch->buf.resize(n_prts);
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

  sub->patch[p].push_back(prt);
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_iter_t

struct psc_particle_PTYPE_iter_t
{
  particle_PTYPE_t& operator*()
  {
    return *psc_mparticles_PTYPE_get_one(const_cast<struct psc_mparticles *>(mprts), p, n);
  }

  particle_PTYPE_t& operator[](int m)
  {
    return *psc_mparticles_PTYPE_get_one(const_cast<struct psc_mparticles *>(mprts), p, n + m);
  }

  bool operator==(const psc_particle_PTYPE_iter_t& other)
  {
    assert(mprts == other.mprts && p == other.p);
    return n == other.n;
  }
  
  bool operator!=(const psc_particle_PTYPE_iter_t& other)
  {
    return !(*this == other);
  }

  psc_particle_PTYPE_iter_t operator++()
  {
    n++;
    return *this;;
  }
  
  int n;
  int p;
  const struct psc_mparticles *mprts;
};

// ----------------------------------------------------------------------
// psc_particle_PTYPE_range_t

struct psc_particle_PTYPE_range_t
{
  psc_particle_PTYPE_iter_t begin;
  psc_particle_PTYPE_iter_t end;

  unsigned int size() { return end.n - begin.n; }

  particle_PTYPE_t& operator[](int m)
  {
    assert(begin.n == 0);
    return *psc_mparticles_PTYPE_get_one(const_cast<struct psc_mparticles *>(begin.mprts), begin.p, m);
  }
};

// ----------------------------------------------------------------------
// psc_mparticles_PTYPE_patch::range

inline psc_particle_PTYPE_range_t psc_mparticles_PTYPE_patch::range()
{
  psc_particle_PTYPE_range_t rv;
  rv.begin.n     = 0;
  rv.begin.p     = p;
  rv.begin.mprts = mprts;
  rv.end.n       = size();
  rv.end.p       = p;
  rv.end.mprts   = mprts;

  return rv;
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
#undef psc_mparticles_PTYPE_get_one
#undef psc_mparticles_PTYPE_patch_push_back
#undef psc_mparticles_PTYPE_patch_resize
#undef psc_mparticles_PTYPE_patch_capacity
#undef psc_particle_PTYPE_iter_t
#undef psc_particle_PTYPE_range_t 

