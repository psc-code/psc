
#include "psc_cuda2.h"
#include "psc_particles_as_cuda2.h"
#include "psc_fields_as_cuda2.h"

#define THREADS_PER_BLOCK (512)

#define F3_CURR F3_DEV
typedef fields_t::real_t *flds_curr_t;
#define F3_EM F3_DEV
typedef fields_t::real_t *flds_em_t;

#include "../psc_push_particles/inc_params.c"
#include "../psc_push_particles/inc_cache.c"
#include "../psc_push_particles/inc_interpolate.c"
#include "../psc_push_particles/inc_push.c"
#include "../psc_push_particles/inc_curr.c"
#include "../psc_push_particles/inc_step.c"

// ----------------------------------------------------------------------
// find_patch

CUDA_DEVICE static int
find_patch()
{
  return blockIdx.z / prm.b_mx[2];
}

// ----------------------------------------------------------------------
// find_ci0
//
// cell index of lower left cell of current block

CUDA_DEVICE static void
find_ci0(int *ci0)
{
  ci0[0] = blockIdx.x * BLOCKSIZE_X;
  ci0[1] = blockIdx.y * BLOCKSIZE_Y;
  ci0[2] = (blockIdx.z % prm.b_mx[2]) * BLOCKSIZE_Z;
}

// ----------------------------------------------------------------------
// find_bid

CUDA_DEVICE static int
find_bid()
{
  return (blockIdx.z * prm.b_mx[1] + blockIdx.y) * prm.b_mx[0] + blockIdx.x;
}

// ----------------------------------------------------------------------
// push_mprts_ab

CUDA_GLOBAL static void CUDA_LAUNCH_BOUNDS(THREADS_PER_BLOCK, 3)
push_mprts_ab(mprts_array_t mprts_arr,
	      unsigned int *b_off,
	      float *d_flds0, unsigned int size)
{
  CUDA_SHARED fields_t::real_t *flds_em;
  CUDA_SHARED curr_cache_t curr_cache;
  {
    int ci0[3]; find_ci0(ci0);
    int p = find_patch();
    fields_t::real_t *d_flds = d_flds0 + p * size;
    flds_em = em_cache_create(d_flds, ci0);
    curr_cache = curr_cache_create(d_flds, ci0);
  }

  int block_begin;
  CUDA_SHARED int block_end;
  {
    int bid = find_bid();
    block_begin = b_off[bid];
    block_end = b_off[bid + 1];
  }

  CUDA_SYNCTHREADS();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    push_one(mprts_arr, n, flds_em, curr_cache);
  }

  {
    int ci0[3]; find_ci0(ci0);
    int p = find_patch();
    fields_t::real_t *d_flds = d_flds0 + p * size;
    curr_cache_destroy(curr_cache, d_flds, ci0);
  }
}

#ifdef __CUDACC__

// ----------------------------------------------------------------------
// zero_currents

static void
zero_currents(struct psc_mfields *mflds)
{
  struct psc_mfields_cuda2 *mflds_sub = psc_mfields_cuda2(mflds);

  unsigned int size = mflds_sub->im[0] * mflds_sub->im[1] * mflds_sub->im[2];

  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_cuda2_real_t *d_flds = mflds_sub->d_flds + p * size * mflds->nr_fields;
    check(cudaMemset(d_flds + JXI * size, 0, 3 * size * sizeof(*d_flds)));
  }
}

#endif

// ----------------------------------------------------------------------
// push_mprts_loop

static void
push_mprts_loop(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);
  struct psc_mfields_cuda2 *mflds_sub = psc_mfields_cuda2(mflds);

  int *bs = mprts_sub->bs;
  assert(bs[0] == BLOCKSIZE_X && bs[1] == BLOCKSIZE_Y && bs[2] == BLOCKSIZE_Z);

  unsigned int fld_size = mflds->nr_fields *
    mflds_sub->im[0] * mflds_sub->im[1] * mflds_sub->im[2];

#ifdef __CUDACC__
  dim3 dimGrid(mprts_sub->b_mx[0],
	       mprts_sub->b_mx[1],
	       mprts_sub->b_mx[2] * mprts->nr_patches);

  mprts_array_t mprts_arr = { .xi4 = mprts_sub->d_xi4, .pxi4 = mprts_sub->d_pxi4, };

  push_mprts_ab<<<dimGrid, THREADS_PER_BLOCK>>>
    (mprts_arr, mprts_sub->d_b_off,
     mflds_sub->d_flds, fld_size);

  cuda_sync_if_enabled();
#else
  int dimGrid[3] = { mprts_sub->b_mx[0],
		     mprts_sub->b_mx[1],
		     mprts_sub->b_mx[2] * mprts->nr_patches };

  mprts_array_t mprts_arr = { .xi4 = mprts_sub->h_xi4, .pxi4 = mprts_sub->h_pxi4, };

  for (blockIdx.z = 0; blockIdx.z < dimGrid[2]; blockIdx.z++) {
    for (blockIdx.y = 0; blockIdx.y < dimGrid[1]; blockIdx.y++) {
      for (blockIdx.x = 0; blockIdx.x < dimGrid[0]; blockIdx.x++) {
	for (threadIdx.x = 0; threadIdx.x < THREADS_PER_BLOCK; threadIdx.x++) {
	  push_mprts_ab(mprts_arr, mprts_sub->h_b_off,
			mflds_sub->h_flds, fld_size);
	}
      }
    }
  }
#endif
}

// ----------------------------------------------------------------------
// cuda2_1vbec_push_mprts

void
SFX(cuda2_1vbec_push_mprts)(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  c_prm_set(ppsc);
  params_1vb_set(ppsc, mprts, mflds);

#ifdef __CUDACC__
  zero_currents(mflds);
#else
  psc_mfields_zero_range(mflds, JXI, JXI + 3);
#endif

  push_mprts_loop(mprts, mflds);
}

