
#include "cuda_mparticles.h"
#include "cuda_mfields.h"

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

// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run

template<typename BS, bool REORDER>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
rho_1st_nc_cuda_run(DMparticlesCuda<BS> dmprts,
		    float4 *d_xi4, float4 *d_pxi4,
		    uint *d_off, int nr_total_blocks, uint *d_ids,
		    DMFields d_flds0)
{
  int block_pos[3];
  int p = dmprts.find_block_pos_patch(block_pos);
  int bid = dmprts.find_bid();
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

    float fnq = prt.qni_wni * dmprts.fnqs();
    
    int lf[3];
    float of[3];
    dmprts.find_idx_off_1st(prt.xi, lf, of, float(0.));

    scurr.add(0, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnq);
    scurr.add(0, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnq);
    scurr.add(0, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnq);
    scurr.add(0, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnq);
  }
}

// ----------------------------------------------------------------------
// n_1st_cuda_run

template<typename BS, bool REORDER>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
n_1st_cuda_run(DMparticlesCuda<BS> dmprts,
	       float4 *d_xi4, float4 *d_pxi4,
	       uint *d_off, int nr_total_blocks, uint *d_ids,
	       DMFields d_flds0)
{
  int block_pos[3];
  int p = dmprts.find_block_pos_patch(block_pos);
  int bid = dmprts.find_bid();
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
    float wni = prt.qni_wni * dmprts.q_inv(kind);
    float fnq = wni * dmprts.fnqs();
    
    int lf[3];
    float of[3];
    dmprts.find_idx_off_1st(prt.xi, lf, of, float(-.5));

    scurr.add(kind, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnq);
    scurr.add(kind, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnq);
    scurr.add(kind, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnq);
    scurr.add(kind, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnq);
  }
}

// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run_patches_no_reorder

template<typename BS, bool REORDER>
static void
rho_1st_nc_cuda_run_patches_no_reorder(cuda_mparticles<BS>* cmprts, struct cuda_mfields *cmres)
{
  dim3 dimGrid(cmprts->b_mx()[1], cmprts->b_mx()[2] * cmprts->n_patches);

  rho_1st_nc_cuda_run<BS, REORDER>
    <<<dimGrid, THREADS_PER_BLOCK>>>
    (*cmprts, cmprts->d_xi4.data().get(), cmprts->d_pxi4.data().get(),
     cmprts->d_off.data().get(),
     cmprts->n_blocks, cmprts->d_id.data().get(), *cmres);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// n_1st_cuda_run_patches_no_reorder

template<typename BS, bool REORDER>
static void
n_1st_cuda_run_patches_no_reorder(cuda_mparticles<BS>* cmprts, struct cuda_mfields *cmres)
{
  dim3 dimGrid(cmprts->b_mx()[1], cmprts->b_mx()[2] * cmprts->n_patches);

  n_1st_cuda_run<BS, REORDER>
    <<<dimGrid, THREADS_PER_BLOCK>>>
    (*cmprts, cmprts->d_xi4.data().get(), cmprts->d_pxi4.data().get(), cmprts->d_off.data().get(),
     cmprts->n_blocks, cmprts->d_id.data().get(), *cmres);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// cuda_moments_yz_rho_1st_nc

template<typename BS>
void cuda_moments_yz_rho_1st_nc(cuda_mparticles<BS>* cmprts, struct cuda_mfields *cmres)
{
  cmprts->reorder(); // FIXME/OPT?

  if (!cmprts->need_reorder) {
    rho_1st_nc_cuda_run_patches_no_reorder<BS, false>(cmprts, cmres);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// cuda_moments_yz_n_1st

template<typename BS>
void cuda_moments_yz_n_1st(cuda_mparticles<BS>* cmprts, struct cuda_mfields *cmres)
{
  cmprts->reorder(); // FIXME/OPT?

  if (!cmprts->need_reorder) {
    n_1st_cuda_run_patches_no_reorder<BS, false>(cmprts, cmres);
  } else {
    assert(0);
  }
}

template void cuda_moments_yz_rho_1st_nc<BS144>(cuda_mparticles<BS144>* cmprts, struct cuda_mfields *cmres);
template void cuda_moments_yz_n_1st<BS144>(cuda_mparticles<BS144>* cmprts, struct cuda_mfields *cmres);
