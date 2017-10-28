
#include "cuda_mparticles.h"
#include "cuda_mfields.h"
#include "cuda_mfields_const.h"

#include "psc_cuda.h"
#include "particles_cuda.h"

#undef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK (512)

#define NO_CHECKERBOARD
//#define DEBUG

#include "cuda_common.h"

#define MAX_KINDS (4)

__constant__ float c_q_inv[MAX_KINDS];

static void
set_constants()
{
  struct psc *psc = ppsc;
  
  assert(psc->nr_kinds <= MAX_KINDS);
  float q_inv[MAX_KINDS];
  for (int k = 0; k < psc->nr_kinds; k++) {
    q_inv[k] = 1.f / psc->kinds[k].q;
  }

  cudaError_t ierr = cudaMemcpyToSymbol(c_q_inv, q_inv, MAX_KINDS * sizeof(*q_inv));
  cudaCheck(ierr);
}
  
// FIXME/TODO: we could do this w/o prior reordering, but currently the
// generic moment calculation code first reorders anyway (which it shouldn't)

// ======================================================================
// GCurr

class GCurr {
public:
  real *d_flds;

  __device__ GCurr(real *_d_flds) :
    d_flds(_d_flds)
  {
  }

  __device__ void add(int m, int jy, int jz, float val)
  {
    float *addr = &D_F3(d_flds, m, 0,jy,jz);
    atomicAdd(addr, val);
  }
};

#define LOAD_PARTICLE_POS_(pp, d_xi4, n) do {				\
    float4 _xi4 = d_xi4[n];						\
    (pp).xi[0]         = _xi4.x;					\
    (pp).xi[1]         = _xi4.y;					\
    (pp).xi[2]         = _xi4.z;					\
    (pp).kind_as_float = _xi4.w;					\
} while (0)

#define LOAD_PARTICLE_MOM_(pp, d_pxi4, n) do {				\
    float4 _pxi4 = d_pxi4[n];						\
    (pp).pxi[0]        = _pxi4.x;					\
    (pp).pxi[1]        = _pxi4.y;					\
    (pp).pxi[2]        = _pxi4.z;					\
    (pp).qni_wni       = _pxi4.w;					\
} while (0)

__device__ static void
find_idx_off_1st(const real xi[3], int j[3], real h[3], real shift,
		 struct cuda_mparticles_params mprts_prm)
{
  for (int d = 0; d < 3; d++) {
    real pos = xi[d] * mprts_prm.dxi[d] + shift;
    j[d] = __float2int_rd(pos);
    h[d] = pos - j[d];
  }
}

// ======================================================================

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static int
find_block_pos_patch(struct cuda_mparticles_params mprts_prm, int *block_pos)
{
  block_pos[1] = blockIdx.x;
  block_pos[2] = blockIdx.y % mprts_prm.b_mx[2];

  return blockIdx.y / mprts_prm.b_mx[2];
}

__device__ static int
find_bid(struct cuda_mparticles_params mprts_prm)
{
  return blockIdx.y * mprts_prm.b_mx[1] + blockIdx.x;
}

// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
rho_1st_nc_cuda_run(int block_start, struct cuda_mparticles_params mprts_prm,
		    float4 *d_xi4, float4 *d_pxi4,
		    unsigned int *d_off, int nr_total_blocks, unsigned int *d_ids,
		    float *d_flds0, unsigned int size)
{
  int block_pos[3];
  int p = find_block_pos_patch<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    (mprts_prm, block_pos);
  int bid = find_bid(mprts_prm);
  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  GCurr scurr(d_flds0 + p * size);

  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    if (REORDER) {
      unsigned int id = d_ids[n];
      LOAD_PARTICLE_POS_(prt, d_xi4, id);
      LOAD_PARTICLE_MOM_(prt, d_pxi4, id);
    } else {
      LOAD_PARTICLE_POS_(prt, d_xi4, n);
      LOAD_PARTICLE_MOM_(prt, d_pxi4, n);
    }

    real fnq = prt.qni_wni * mprts_prm.fnqs;
    
    int lf[3];
    real of[3];
    find_idx_off_1st(prt.xi, lf, of, real(0.), mprts_prm);

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
n_1st_cuda_run(int block_start, struct cuda_mparticles_params mprts_prm,
	       float4 *d_xi4, float4 *d_pxi4,
	       unsigned int *d_off, int nr_total_blocks, unsigned int *d_ids,
	       float *d_flds0, unsigned int size)
{
  int block_pos[3];
  int p = find_block_pos_patch<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    (mprts_prm, block_pos);
  int bid = find_bid(mprts_prm);
  int block_begin = d_off[bid];
  int block_end = d_off[bid + 1];

  GCurr scurr(d_flds0 + p * size);

  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    if (REORDER) {
      unsigned int id = d_ids[n];
      LOAD_PARTICLE_POS_(prt, d_xi4, id);
      LOAD_PARTICLE_MOM_(prt, d_pxi4, id);
    } else {
      LOAD_PARTICLE_POS_(prt, d_xi4, n);
      LOAD_PARTICLE_MOM_(prt, d_pxi4, n);
    }

    int kind = __float_as_int(prt.kind_as_float);
    real wni = prt.qni_wni * c_q_inv[kind];
    real fnq = wni * mprts_prm.fnqs;
    
    int lf[3];
    real of[3];
    find_idx_off_1st(prt.xi, lf, of, real(-.5), mprts_prm);

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
rho_1st_nc_cuda_run_patches_no_reorder(struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
  struct psc_mfields_cuda *mres_cuda = psc_mfields_cuda(mres);
  struct cuda_mfields *cmres = mres_cuda->cmflds;

  cuda_mfields_const_set(cmres);
  
  struct cuda_mparticles_params mprts_prm;
  cuda_mparticles_params_set(&mprts_prm, cmprts);

  unsigned int fld_size = mres->nr_fields * cmres->im[0] * cmres->im[1] * cmres->im[2];

  int gx = mprts_prm.b_mx[1];
  int gy = mprts_prm.b_mx[2] * mprts->nr_patches;
  dim3 dimGrid(gx, gy);

  rho_1st_nc_cuda_run<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, REORDER>
    <<<dimGrid, THREADS_PER_BLOCK>>>
    (0, mprts_prm,
     cmprts->d_xi4, cmprts->d_pxi4,
     cmprts->d_off,
     cmprts->n_blocks, cmprts->d_id,
     cmres->d_flds, fld_size);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// n_1st_cuda_run_patches_no_reorder

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER>
static void
n_1st_cuda_run_patches_no_reorder(struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
  struct psc_mfields_cuda *mres_cuda = psc_mfields_cuda(mres);
  struct cuda_mfields *cmres = mres_cuda->cmflds;

  cuda_mfields_const_set(cmres);

  struct cuda_mparticles_params mprts_prm;
  cuda_mparticles_params_set(&mprts_prm, cmprts);

  unsigned int fld_size = mres->nr_fields *
    cmres->im[0] * cmres->im[1] * cmres->im[2];

  int gx = mprts_prm.b_mx[1];
  int gy = mprts_prm.b_mx[2] * mprts->nr_patches;
  dim3 dimGrid(gx, gy);

  n_1st_cuda_run<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, REORDER>
    <<<dimGrid, THREADS_PER_BLOCK>>>
    (0, mprts_prm,
     cmprts->d_xi4, cmprts->d_pxi4, cmprts->d_off,
     cmprts->n_blocks, cmprts->d_id,
     cmres->d_flds, fld_size);
  cuda_sync_if_enabled();
}

// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run_patches

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
rho_1st_nc_cuda_run_patches(struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
    
  cuda_mparticles_reorder(cmprts); // FIXME/OPT?

  if (!cmprts->need_reorder) {
    rho_1st_nc_cuda_run_patches_no_reorder<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, false>(mprts, mres);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// n_1st_cuda_run_patches

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
static void
n_1st_cuda_run_patches(struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  struct cuda_mparticles *cmprts = psc_mparticles_cuda(mprts)->cmprts;
    
  cuda_mparticles_reorder(cmprts); // FIXME/OPT?

  if (!cmprts->need_reorder) {
    n_1st_cuda_run_patches_no_reorder<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, false>(mprts, mres);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// yz_moments_rho_1st_nc_cuda_run_patches

void
yz_moments_rho_1st_nc_cuda_run_patches(struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  rho_1st_nc_cuda_run_patches<1, 4, 4>(mprts, mres);
}

// ----------------------------------------------------------------------
// yz_moments_n_1st_cuda_run_patches

void
yz_moments_n_1st_cuda_run_patches(struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  static bool first_time = true;
  if (first_time) {
    set_constants();
    first_time = false;
  }
  n_1st_cuda_run_patches<1, 4, 4>(mprts, mres);
}

