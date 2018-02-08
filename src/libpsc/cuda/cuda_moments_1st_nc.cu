
#include "cuda_mparticles.h"
#include "cuda_mfields.h"
#include "cuda_mparticles_const.h"
#include "cuda_mfields_const.h"

#define THREADS_PER_BLOCK (512)

// FIXME/TODO: we could do this w/o prior reordering, but currently the
// generic moment calculation code first reorders anyway (which it shouldn't)

// ======================================================================
// GCurr

class GCurr {
public:
  DFields d_flds;

  __device__ GCurr(DFields _d_flds) :
    d_flds(_d_flds)
  {
  }

  __device__ void add(int m, int jy, int jz, float val)
  {
    float *addr = &d_flds(m, 0,jy,jz);
    atomicAdd(addr, val);
  }
};

// ======================================================================

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static int
find_block_pos_patch(int *block_pos)
{
  block_pos[1] = blockIdx.x;
  block_pos[2] = blockIdx.y % d_cmprts_const.b_mx[2];

  return blockIdx.y / d_cmprts_const.b_mx[2];
}

// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
rho_1st_nc_cuda_run(int block_start,
		    float4 *d_xi4, float4 *d_pxi4,
		    uint *d_off, int nr_total_blocks, uint *d_ids,
		    DMFields d_flds0)
{
  int block_pos[3];
  int p = find_block_pos_patch<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(block_pos);
  int bid = find_bid();
  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  GCurr scurr(d_flds0[p]);

  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    if (REORDER) {
      uint id = d_ids[n];
      LOAD_PARTICLE_POS(prt, d_xi4, id);
      LOAD_PARTICLE_MOM(prt, d_pxi4, id);
    } else {
      LOAD_PARTICLE_POS(prt, d_xi4, n);
      LOAD_PARTICLE_MOM(prt, d_pxi4, n);
    }

    float fnq = prt.qni_wni * d_cmprts_const.fnqs;
    
    int lf[3];
    float of[3];
    find_idx_off_1st(prt.xi, lf, of, float(0.));

    scurr.add(0, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnq);
    scurr.add(0, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnq);
    scurr.add(0, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnq);
    scurr.add(0, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnq);
  }
}

// ----------------------------------------------------------------------
// n_1st_cuda_run

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
n_1st_cuda_run(int block_start,
	       float4 *d_xi4, float4 *d_pxi4,
	       uint *d_off, int nr_total_blocks, uint *d_ids,
	       DMFields d_flds0)
{
  int block_pos[3];
  int p = find_block_pos_patch<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>(block_pos);
  int bid = find_bid();
  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  GCurr scurr(d_flds0[p]);

  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    if (REORDER) {
      uint id = d_ids[n];
      LOAD_PARTICLE_POS(prt, d_xi4, id);
      LOAD_PARTICLE_MOM(prt, d_pxi4, id);
    } else {
      LOAD_PARTICLE_POS(prt, d_xi4, n);
      LOAD_PARTICLE_MOM(prt, d_pxi4, n);
    }

    int kind = __float_as_int(prt.kind_as_float);
    float wni = prt.qni_wni * d_cmprts_const.q_inv[kind];
    float fnq = wni * d_cmprts_const.fnqs;
    
    int lf[3];
    float of[3];
    find_idx_off_1st(prt.xi, lf, of, float(-.5));

    scurr.add(kind, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnq);
    scurr.add(kind, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnq);
    scurr.add(kind, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnq);
    scurr.add(kind, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnq);
  }
}

// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run_patches_no_reorder

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER>
static void
rho_1st_nc_cuda_run_patches_no_reorder(struct cuda_mparticles *cmprts, struct cuda_mfields *cmres)
{
  cuda_mparticles_const_set(cmprts);
  cuda_mfields_const_set(cmres);
  
  dim3 dimGrid(cmprts->b_mx_[1], cmprts->b_mx_[2] * cmprts->n_patches);

  rho_1st_nc_cuda_run<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, REORDER>
    <<<dimGrid, THREADS_PER_BLOCK>>>
    (0, cmprts->d_xi4.data().get(), cmprts->d_pxi4.data().get(),
     cmprts->d_off.data().get(),
     cmprts->n_blocks, cmprts->d_id.data().get(), *cmres);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// n_1st_cuda_run_patches_no_reorder

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER>
static void
n_1st_cuda_run_patches_no_reorder(struct cuda_mparticles *cmprts, struct cuda_mfields *cmres)
{
  cuda_mparticles_const_set(cmprts);
  cuda_mfields_const_set(cmres);

  dim3 dimGrid(cmprts->b_mx_[1], cmprts->b_mx_[2] * cmprts->n_patches);

  n_1st_cuda_run<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, REORDER>
    <<<dimGrid, THREADS_PER_BLOCK>>>
    (0, cmprts->d_xi4.data().get(), cmprts->d_pxi4.data().get(), cmprts->d_off.data().get(),
     cmprts->n_blocks, cmprts->d_id.data().get(), *cmres);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run_patches

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
rho_1st_nc_cuda_run_patches(struct cuda_mparticles *cmprts, struct cuda_mfields *cmres)
{
  cmprts->reorder(); // FIXME/OPT?

  if (!cmprts->need_reorder) {
    rho_1st_nc_cuda_run_patches_no_reorder<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, false>(cmprts, cmres);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// n_1st_cuda_run_patches

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
n_1st_cuda_run_patches(struct cuda_mparticles *cmprts, struct cuda_mfields *cmres)
{
  cmprts->reorder(); // FIXME/OPT?

  if (!cmprts->need_reorder) {
    n_1st_cuda_run_patches_no_reorder<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, false>(cmprts, cmres);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// cuda_moments_yz_rho_1st_nc

void
cuda_moments_yz_rho_1st_nc(struct cuda_mparticles *cmprts, struct cuda_mfields *cmres)
{
  rho_1st_nc_cuda_run_patches<1, 4, 4>(cmprts, cmres);
}

// ----------------------------------------------------------------------
// cuda_moments_yz_n_1st

void
cuda_moments_yz_n_1st(struct cuda_mparticles *cmprts, struct cuda_mfields *cmres)
{
  n_1st_cuda_run_patches<1, 4, 4>(cmprts, cmres);
}

