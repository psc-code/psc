
#include "psc_particle_common.h"

#include <stdlib.h>
#include <assert.h>
#include <mrc_bits.h>
#include <algorithm>

#if PTYPE == PTYPE_SINGLE

#define particle_PTYPE_real_t particle_single_real_t
#define particle_PTYPE_t particle_single_t

#define psc_particle_PTYPE_buf_t psc_particle_single_buf_t
#define psc_particle_PTYPE_buf_ctor psc_particle_single_buf_ctor
#define psc_particle_PTYPE_buf_dtor psc_particle_single_buf_dtor
#define psc_particle_PTYPE_buf_reserve psc_particle_single_buf_reserve
#define psc_particle_PTYPE_buf_capacity psc_particle_single_buf_capacity
#define psc_particle_PTYPE_buf_push_back psc_particle_single_buf_push_back
#define psc_particle_PTYPE_buf_at_ptr psc_particle_single_buf_at_ptr
#define psc_particle_PTYPE_range_t psc_particle_single_range_t

#elif PTYPE == PTYPE_DOUBLE

#define particle_PTYPE_real_t particle_double_real_t
#define particle_PTYPE_t particle_double_t

#define psc_particle_PTYPE_buf_t psc_particle_double_buf_t
#define psc_particle_PTYPE_buf_ctor psc_particle_double_buf_ctor
#define psc_particle_PTYPE_buf_dtor psc_particle_double_buf_dtor
#define psc_particle_PTYPE_buf_reserve psc_particle_double_buf_reserve
#define psc_particle_PTYPE_buf_capacity psc_particle_double_buf_capacity
#define psc_particle_PTYPE_buf_push_back psc_particle_double_buf_push_back
#define psc_particle_PTYPE_buf_at_ptr psc_particle_double_buf_at_ptr
#define psc_particle_PTYPE_range_t psc_particle_double_range_t

#elif PTYPE == PTYPE_SINGLE_BY_BLOCK

#define particle_PTYPE_real_t particle_single_by_block_real_t
#define particle_PTYPE_t particle_single_by_block_t

#define psc_particle_PTYPE_buf_t psc_particle_single_by_block_buf_t
#define psc_particle_PTYPE_buf_ctor psc_particle_single_by_block_buf_ctor
#define psc_particle_PTYPE_buf_dtor psc_particle_single_by_block_buf_dtor
#define psc_particle_PTYPE_buf_reserve psc_particle_single_by_block_buf_reserve
#define psc_particle_PTYPE_buf_capacity psc_particle_single_by_block_buf_capacity
#define psc_particle_PTYPE_buf_push_back psc_particle_single_by_block_buf_push_back
#define psc_particle_PTYPE_buf_at_ptr psc_particle_single_by_block_buf_at_ptr
#define psc_particle_PTYPE_range_t psc_particle_single_by_block_range_t

#elif PTYPE == PTYPE_FORTRAN

#define particle_PTYPE_real_t particle_fortran_real_t
#define particle_PTYPE_t particle_fortran_t

#define psc_particle_PTYPE_buf_t psc_particle_fortran_buf_t
#define psc_particle_PTYPE_buf_ctor psc_particle_fortran_buf_ctor
#define psc_particle_PTYPE_buf_dtor psc_particle_fortran_buf_dtor
#define psc_particle_PTYPE_buf_reserve psc_particle_fortran_buf_reserve
#define psc_particle_PTYPE_buf_capacity psc_particle_fortran_buf_capacity
#define psc_particle_PTYPE_buf_push_back psc_particle_fortran_buf_push_back
#define psc_particle_PTYPE_buf_at_ptr psc_particle_fortran_buf_at_ptr

#elif PTYPE == PTYPE_CUDA

#define particle_PTYPE_real_t particle_cuda_real_t
#define particle_PTYPE_t particle_cuda_t

#define psc_particle_PTYPE_buf_t psc_particle_cuda_buf_t
#define psc_particle_PTYPE_buf_ctor psc_particle_cuda_buf_ctor
#define psc_particle_PTYPE_buf_dtor psc_particle_cuda_buf_dtor
#define psc_particle_PTYPE_buf_reserve psc_particle_cuda_buf_reserve
#define psc_particle_PTYPE_buf_capacity psc_particle_cuda_buf_capacity
#define psc_particle_PTYPE_buf_push_back psc_particle_cuda_buf_push_back
#define psc_particle_PTYPE_buf_at_ptr psc_particle_cuda_buf_at_ptr

#endif

// ======================================================================
// psc_particle_PTYPE_buf_t

struct psc_particle_PTYPE_range_t;

struct psc_particle_PTYPE_buf_t
{
  using particle_t = particle_PTYPE_t;
#ifdef psc_particle_PTYPE_range_t
  using range_t = psc_particle_PTYPE_range_t;
#endif
  
  particle_PTYPE_t *m_data;
  unsigned int m_size;
  unsigned int m_capacity;

  unsigned int size() const { return m_size; }
  unsigned int capacity() const { return m_capacity; }

  void resize(unsigned int new_size)
  {
    assert(new_size <= m_capacity);
    m_size = new_size;
  }

  void reserve(unsigned int new_capacity)
  {
    if (new_capacity <= m_capacity)
      return;

    new_capacity = std::max(new_capacity, m_capacity * 2);
    
    m_data = (particle_PTYPE_t *) realloc(m_data, new_capacity * sizeof(*m_data));
    m_capacity = new_capacity;
  }

  void push_back(const particle_PTYPE_t& prt)
  {
    unsigned int n = m_size;
    if (n >= m_capacity) {
      reserve(n + 1);
    }
    m_data[n++] = prt;
    m_size = n;
  }

  particle_PTYPE_t& operator[](int n)
  {
    return m_data[n];
  }
};

// ----------------------------------------------------------------------
// psc_particle_PTYPE_buf_ctor

static inline void
psc_particle_PTYPE_buf_ctor(psc_particle_PTYPE_buf_t *buf)
{
  buf->m_data = NULL;
  buf->m_size = 0;
  buf->m_capacity = 0;
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_buf_dtor

static inline void
psc_particle_PTYPE_buf_dtor(psc_particle_PTYPE_buf_t *buf)
{
  free(buf->m_data);
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_buf_reserve

static inline void
psc_particle_PTYPE_buf_reserve(psc_particle_PTYPE_buf_t *buf, unsigned int new_capacity)
{
  if (new_capacity <= buf->m_capacity)
    return;

  new_capacity = MAX(new_capacity, buf->m_capacity * 2);
  
  buf->m_data = (particle_PTYPE_t *) realloc(buf->m_data, new_capacity * sizeof(*buf->m_data));
  buf->m_capacity = new_capacity;
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_buf_capacity

static inline unsigned int
psc_particle_PTYPE_buf_capacity(const psc_particle_PTYPE_buf_t *buf)
{
  return buf->m_capacity;
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_buf_push_back

static inline void
psc_particle_PTYPE_buf_push_back(psc_particle_PTYPE_buf_t *buf, particle_PTYPE_t prt)
{
  unsigned int n = buf->m_size;
  if (n >= buf->m_capacity) {
    psc_particle_PTYPE_buf_reserve(buf, n + 1);
  }
  buf->m_data[n++] = prt;
  buf->m_size = n;
}

// ----------------------------------------------------------------------
// psc_particle_PTYPE_buf_at_ptr

static inline particle_PTYPE_t *
psc_particle_PTYPE_buf_at_ptr(psc_particle_PTYPE_buf_t *buf, unsigned int n)
{
  // FIXME? could do bounds check here...
  return &buf->m_data[n];
}


#undef particle_PTYPE_real_t
#undef particle_PTYPE_t

#undef psc_particle_PTYPE_buf_t
#undef psc_particle_PTYPE_buf_ctor
#undef psc_particle_PTYPE_buf_dtor
#undef psc_particle_PTYPE_buf_reserve
#undef psc_particle_PTYPE_buf_capacity
#undef psc_particle_PTYPE_buf_push_back
#undef psc_particle_PTYPE_buf_at_ptr
#undef psc_particle_PTYPE_range_t

