
#include "psc_cuda2.h"
#include "psc_particles_as_cuda2.h"
#include "psc_fields_cuda2.h"

#define THREADS_PER_BLOCK (512)

#include "../psc_push_particles/inc_params.c"
#include "../psc_push_particles/inc_cache.c"
#include "../psc_push_particles/inc_interpolate.c"
#include "../psc_push_particles/inc_push.c"
#include "../psc_push_particles/inc_curr.c"
#include "../psc_push_particles/inc_step.c"

// ----------------------------------------------------------------------
// find_block_pos_patch

CUDA_DEVICE static int
find_block_pos_patch(int *ci0)
{
#if EM_CACHE == EM_CACHE_CUDA
  ci0[0] = block_Idx.x * BLOCKSIZE_X;
  ci0[1] = block_Idx.y * BLOCKSIZE_Y;
  ci0[2] = (blockIdx.z % prm.b_mx[2]) * BLOCKSIZE_Z;
#endif

  return blockIdx.z / prm.b_mx[2];
}

// ----------------------------------------------------------------------
// find_bid

CUDA_DEVICE static int
find_bid()
{
  return (blockIdx.z * prm.b_mx[1] + blockIdx.y) * prm.b_mx[0] + blockIdx.x;
}

#ifdef __CUDACC__

// ----------------------------------------------------------------------
// push_mprts_ab

__global__ static void __launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_ab(float4 *d_xi4, float4 *d_pxi4,
	      unsigned int *d_off,
	      float *d_flds0, unsigned int size)
{
  int ci0[3];
  int p, bid;
  p = find_block_pos_patch(ci0);
  real *d_flds = d_flds0 + p * size;

  DECLARE_EM_CACHE(flds_em, d_flds, size, ci0);
  real *flds_curr = d_flds;

  bid = find_bid();
  __shared__ int block_begin, block_end;
  block_begin = d_off[bid];
  block_end = d_off[bid + 1];

  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    push_one_mprts(d_xi4, d_pxi4, n, flds_em, flds_curr);
  }
}

// ----------------------------------------------------------------------
// push_mprts_a

__global__ static void __launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_a(float4 *d_xi4, float4 *d_pxi4,
	     unsigned int *d_off,
	     float *d_flds0, unsigned int size)
{
  int ci0[3];
  int p, bid;
  p = find_block_pos_patch(ci0);
  real *d_flds = d_flds0 + p * size;

  DECLARE_EM_CACHE(flds_em, d_flds, size, ci0);
  __shared__ real *flds_curr;
  flds_curr = d_flds;

  bid = find_bid();
  __shared__ int block_begin, block_end;
  block_begin = d_off[bid];
  block_end = d_off[bid + 1];

  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    push_one_mprts_a(d_xi4, d_pxi4, n, flds_em, flds_curr);
  }
}

// ----------------------------------------------------------------------
// push_mprts_b

__global__ static void __launch_bounds__(THREADS_PER_BLOCK, 3)
push_mprts_b(float4 *d_xi4, float4 *d_pxi4,
	     unsigned int *d_off,
	     float *d_flds0, unsigned int size)
{
  int ci0[3];
  int p, bid;
  p = find_block_pos_patch(ci0);
  real *d_flds = d_flds0 + p * size;

  real *flds_curr = d_flds;

  bid = find_bid();
  __shared__ int block_begin, block_end;
  block_begin = d_off[bid];
  block_end = d_off[bid + 1];

  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    push_one_mprts_b(d_xi4, d_pxi4, n, flds_curr);
  }
}

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
#ifdef __CUDACC__
  struct psc_mfields_cuda2 *mflds_sub = psc_mfields_cuda2(mflds);

  int *bs = mprts_sub->bs;
  assert(bs[0] == BLOCKSIZE_X && bs[1] == BLOCKSIZE_Y && bs[2] == BLOCKSIZE_Z);

  unsigned int fld_size = mflds->nr_fields *
    mflds_sub->im[0] * mflds_sub->im[1] * mflds_sub->im[2];

  dim3 dimGrid(mprts_sub->b_mx[0],
	       mprts_sub->b_mx[1],
	       mprts_sub->b_mx[2] * mprts->nr_patches);

  push_mprts_ab<<<dimGrid, THREADS_PER_BLOCK>>>
    (mprts_sub->d_xi4, mprts_sub->d_pxi4, mprts_sub->d_b_off,
     mflds_sub->d_flds, fld_size);

  cuda_sync_if_enabled();
#else
  int dimGrid[3] = { mprts_sub->b_mx[0],
		     mprts_sub->b_mx[1],
		     mprts_sub->b_mx[2] * mprts->nr_patches };

  for (blockIdx.z = 0; blockIdx.z < dimGrid[2]; blockIdx.z++) {
    for (blockIdx.y = 0; blockIdx.y < dimGrid[1]; blockIdx.y++) {
      for (blockIdx.x = 0; blockIdx.x < dimGrid[0]; blockIdx.x++) {
	for (threadIdx.x = 0; threadIdx.x < THREADS_PER_BLOCK; threadIdx.x++) {
	  int ci0[3];
	  int p = find_block_pos_patch(ci0);
	  int bid = find_bid();
	  int block_begin, block_end;
	  unsigned int *b_off = mprts_sub->h_b_off;
	  block_begin = b_off[bid];
	  block_end = b_off[bid + 1];
	  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
	    if (n < block_begin) {
	      continue;
	    }
	    //	for (int n = block_begin; n < block_end; n++) {
	    push_one_mprts(mprts, mflds, n, p);
	  }
	}
      }
    }
  }
#endif
}

#ifdef __CUDACC__
// ----------------------------------------------------------------------
// push_mprts_loop_b12

static void
push_mprts_loop_b1(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);
  struct psc_mfields_cuda2 *mflds_sub = psc_mfields_cuda2(mflds);

  int *bs = mprts_sub->bs;
  assert(bs[0] == BLOCKSIZE_X && bs[1] == BLOCKSIZE_Y && bs[2] == BLOCKSIZE_Z);

  unsigned int fld_size = mflds->nr_fields *
    mflds_sub->im[0] * mflds_sub->im[1] * mflds_sub->im[2];

  dim3 dimGrid(mprts_sub->b_mx[0],
	       mprts_sub->b_mx[1],
	       mprts_sub->b_mx[2] * mprts->nr_patches);

  push_mprts_a<<<dimGrid, THREADS_PER_BLOCK>>>
    (mprts_sub->d_xi4, mprts_sub->d_pxi4, mprts_sub->d_b_off,
     mflds_sub->d_flds, fld_size);

  cuda_sync_if_enabled();
}

static void
push_mprts_loop_b2(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);
  struct psc_mfields_cuda2 *mflds_sub = psc_mfields_cuda2(mflds);

  int *bs = mprts_sub->bs;
  assert(bs[0] == BLOCKSIZE_X && bs[1] == BLOCKSIZE_Y && bs[2] == BLOCKSIZE_Z);

  unsigned int fld_size = mflds->nr_fields *
    mflds_sub->im[0] * mflds_sub->im[1] * mflds_sub->im[2];

  dim3 dimGrid(mprts_sub->b_mx[0],
	       mprts_sub->b_mx[1],
	       mprts_sub->b_mx[2] * mprts->nr_patches);

  push_mprts_b<<<dimGrid, THREADS_PER_BLOCK>>>
    (mprts_sub->d_xi4, mprts_sub->d_pxi4, mprts_sub->d_b_off,
     mflds_sub->d_flds, fld_size);

  cuda_sync_if_enabled();
}

#endif

// ----------------------------------------------------------------------
// cuda2_1vbec_push_mprts

void
SFX(cuda2_1vbec_push_mprts)(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  params_1vb_set(ppsc, mprts, mflds);

#ifdef __CUDACC__
  zero_currents(mflds);
#else
  psc_mfields_zero_range(mflds, JXI, JXI + 3);
#endif

  push_mprts_loop(mprts, mflds);
}

#ifdef __CUDACC__

void
SFX(cuda2_1vbec_push_mprts_a)(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  params_1vb_set(ppsc, mprts, mflds);

  zero_currents(mflds);
}

void
SFX(cuda2_1vbec_push_mprts_b)(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  push_mprts_loop_b1(mprts, mflds);
  push_mprts_loop_b2(mprts, mflds);
}

void
SFX(cuda2_1vbec_push_mprts_b1)(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  push_mprts_loop_b1(mprts, mflds);
}

void
SFX(cuda2_1vbec_push_mprts_b2)(struct psc_mparticles *mprts, struct psc_mfields *mflds)
{
  push_mprts_loop_b2(mprts, mflds);
}

#endif
