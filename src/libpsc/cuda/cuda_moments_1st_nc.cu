
#include "psc_cuda.h"
#include "particles_cuda.h"

#undef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK (512)

enum IP {
  IP_STD, // standard interpolation
  IP_EC,  // energy-conserving interpolation
};

enum DEPOSIT {
  DEPOSIT_VB_2D,
  DEPOSIT_VB_3D,
};

enum CURRMEM {
  CURRMEM_SHARED,
  CURRMEM_GLOBAL,
};

#define NO_CHECKERBOARD
//#define DEBUG

#include "cuda_common.h"

static __constant__ __device__ float c_dqs[4]; // FIXME hardcoded

// ----------------------------------------------------------------------
// set_consts

static void
set_consts(struct cuda_params *prm)
{
  check(cudaMemcpyToSymbol(c_dqs, prm->dq, sizeof(c_dqs)));
}

// ======================================================================
// GCurr

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
class GCurr {
public:
  static const int shared_size = 1;

  real *scurr;
  real *d_flds;

  __device__ GCurr(real *_scurr, real *_d_flds) :
    scurr(_scurr), d_flds(_d_flds)
  {
  }

  __device__ void add_to_fld(struct cuda_params prm, int *ci0)
  {
  }

  __device__ void add(int m, int jy, int jz, float val, struct cuda_params prm, int *ci0)
  {
    float *addr = &F3_DEV_YZ(JXI+m, jy+ci0[1],jz+ci0[2]);
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

__device__ static void
find_idx_off_1st(const real xi[3], int j[3], real h[3], real shift,
		 struct cuda_params prm)
{
  for (int d = 0; d < 3; d++) {
    real pos = xi[d] * prm.dxi[d] + shift;
    j[d] = __float2int_rd(pos);
    h[d] = pos - j[d];
  }
}

// ----------------------------------------------------------------------
// push_part_one

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER, enum IP IP>
__device__ static void
push_part_one(struct d_particle *prt, int n, unsigned int *d_ids,
	      float4 *d_xi4, float4 *d_pxi4,
	      int ci0[3], struct cuda_params prm)
{
  unsigned int id;
  if (REORDER) {
    id = d_ids[n];
    LOAD_PARTICLE_POS_(*prt, d_xi4, id);
  } else {
    LOAD_PARTICLE_POS_(*prt, d_xi4, n);
  }
  // here we have x^{n+.5}
}

// ----------------------------------------------------------------------
// yz_calc_j

template<enum DEPOSIT DEPOSIT, class CURR>
__device__ static void
yz_calc_j(struct d_particle *prt, int n, float4 *d_xi4, float4 *d_pxi4,
	  CURR &scurr,
	  struct cuda_params prm, int nr_total_blocks, int p_nr,
	  unsigned int *d_bidx, int bid, int *ci0)
{
  real fnq = prt->qni_wni * prm.fnqs;
  
  int lf[3];
  real of[3];
  find_idx_off_1st(prt->xi, lf, of, real(0.), prm);
  lf[1] -= ci0[1];
  lf[2] -= ci0[2];

  scurr.add(0, lf[1]  , lf[2]  , (1.f - of[1]) * (1.f - of[2]) * fnq, prm, ci0);
  scurr.add(0, lf[1]+1, lf[2]  , (      of[1]) * (1.f - of[2]) * fnq, prm, ci0);
  scurr.add(0, lf[1]  , lf[2]+1, (1.f - of[1]) * (      of[2]) * fnq, prm, ci0);
  scurr.add(0, lf[1]+1, lf[2]+1, (      of[1]) * (      of[2]) * fnq, prm, ci0);
}

// ======================================================================

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static int
find_block_pos_patch(struct cuda_params prm, int *block_pos, int *ci0)
{
  block_pos[1] = blockIdx.x;
  block_pos[2] = blockIdx.y % prm.b_mx[2];

  ci0[0] = 0;
  ci0[1] = block_pos[1] * BLOCKSIZE_Y;
  ci0[2] = block_pos[2] * BLOCKSIZE_Z;

  return blockIdx.y / prm.b_mx[2];
}

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static int
find_block_pos_patch_q(struct cuda_params prm, int *block_pos, int *ci0, int block_start)
{
  int grid_dim_y = (prm.b_mx[2] + 1) / 2;
  block_pos[1] = blockIdx.x * 2;
  block_pos[2] = (blockIdx.y % grid_dim_y) * 2;
  block_pos[1] += block_start & 1;
  block_pos[2] += block_start >> 1;
  if (block_pos[1] >= prm.b_mx[1] ||
      block_pos[2] >= prm.b_mx[2])
    return -1;

  ci0[0] = 0;
  ci0[1] = block_pos[1] * BLOCKSIZE_Y;
  ci0[2] = block_pos[2] * BLOCKSIZE_Z;

  return blockIdx.y / grid_dim_y;
}

__device__ static int
find_bid(struct cuda_params prm)
{
  return blockIdx.y * prm.b_mx[1] + blockIdx.x;
}

__device__ static int
find_bid_q(struct cuda_params prm, int p, int *block_pos)
{
  // FIXME won't work if b_mx[1,2] not even (?)
  return block_pos_to_block_idx(block_pos, prm.b_mx) + p * prm.b_mx[1] * prm.b_mx[2];
}

#define DECLARE_AND_ZERO_SCURR						\
  __shared__ real _scurr[CURR::shared_size];				\
  CURR scurr(_scurr, d_flds0 + p * size)				\

#define FIND_BLOCK_RANGE_CURRMEM(CURRMEM)				\
  int block_pos[3], ci0[3];						\
  int p, bid;								\
  if (CURRMEM == CURRMEM_SHARED) {					\
    p = find_block_pos_patch_q<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>	\
      (prm, block_pos, ci0, block_start);				\
    if (p < 0)								\
      return;								\
    									\
    bid = find_bid_q(prm, p, block_pos);				\
  } else if (CURRMEM == CURRMEM_GLOBAL) {				\
    p = find_block_pos_patch<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>	\
      (prm, block_pos, ci0);						\
    bid = find_bid(prm);						\
  }									\
  int block_begin = d_off[bid];						\
  int block_end = d_off[bid + 1]					\
    
// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER,
	 enum IP IP, enum DEPOSIT DEPOSIT, enum CURRMEM CURRMEM, class CURR>
__global__ static void
__launch_bounds__(THREADS_PER_BLOCK, 3)
rho_1st_nc_cuda_run(int block_start, struct cuda_params prm, float4 *d_xi4, float4 *d_pxi4,
	      float4 *d_alt_xi4, float4 *d_alt_pxi4,
	      unsigned int *d_off, int nr_total_blocks, unsigned int *d_ids, unsigned int *d_bidx,
	      float *d_flds0, unsigned int size)
{
  FIND_BLOCK_RANGE_CURRMEM(CURRMEM);
  //  DECLARE_AND_CACHE_FIELDS;
  DECLARE_AND_ZERO_SCURR;
  
  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    push_part_one<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, REORDER, IP>
      (&prt, n, d_ids, d_xi4, d_pxi4, ci0, prm);

    yz_calc_j<DEPOSIT_VB_2D, CURR>
      (&prt, n, d_xi4, d_pxi4, scurr, prm, nr_total_blocks, p, d_bidx, bid, ci0);
  }
  
  scurr.add_to_fld(prm, ci0);
}

// ----------------------------------------------------------------------

#define CUDA_PUSH_MPRTS_TOP(CURRMEM)					\
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);	\
  struct psc_mfields_cuda *mres_cuda = psc_mfields_cuda(mres);		\
									\
  struct cuda_params prm;						\
  set_params(&prm, ppsc, mprts, mres);					\
  set_consts(&prm);							\
									\
  unsigned int fld_size = mres->nr_fields *				\
    mres_cuda->im[0] * mres_cuda->im[1] * mres_cuda->im[2];		\
									\
  int gx, gy;								\
  if (CURRMEM == CURRMEM_GLOBAL) {					\
    gx = prm.b_mx[1];							\
    gy = prm.b_mx[2] * mprts->nr_patches;				\
  }									\
  dim3 dimGrid(gx, gy);							\

// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run_patches_no_reorder

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, bool REORDER,
	 enum IP IP, enum DEPOSIT DEPOSIT, enum CURRMEM CURRMEM>
static void
rho_1st_nc_cuda_run_patches_no_reorder(struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  CUDA_PUSH_MPRTS_TOP(CURRMEM);

  if (CURRMEM == CURRMEM_GLOBAL) {
    rho_1st_nc_cuda_run<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, REORDER, IP, DEPOSIT, CURRMEM,
			GCurr<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z> >
      <<<dimGrid, THREADS_PER_BLOCK>>>
      (0, prm, mprts_cuda->d_xi4, mprts_cuda->d_pxi4,
       mprts_cuda->d_alt_xi4, mprts_cuda->d_alt_pxi4, mprts_cuda->d_off,
       mprts_cuda->nr_total_blocks, mprts_cuda->d_ids, mprts_cuda->d_bidx,
       mres_cuda->d_flds, fld_size);
    cuda_sync_if_enabled();
  } else {
    assert(0);
  }

  free_params(&prm);
}

// ----------------------------------------------------------------------
// rho_1st_nc_cuda_run_patches

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z, enum IP IP, enum DEPOSIT DEPOSIT,
	 enum CURRMEM CURRMEM>
static void
rho_1st_nc_cuda_run_patches(struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  struct psc_mparticles_cuda *mprts_cuda = psc_mparticles_cuda(mprts);
    
  psc_mparticles_cuda_copy_to_dev(mprts);
  
  if (!mprts_cuda->need_reorder) {
    rho_1st_nc_cuda_run_patches_no_reorder<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, false, IP, DEPOSIT, CURRMEM>(mprts, mres);
  } else {
    assert(0);
#if 0
    cuda_push_mprts_ab<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z, true, IP, DEPOSIT, CURRMEM>(mprts, mres);
    mprts_cuda->need_reorder = false;
#endif
  }
}

// ----------------------------------------------------------------------
// yz_moments_rho_1st_nc_cuda_run_patches

void
yz_moments_rho_1st_nc_cuda_run_patches(struct psc_mparticles *mprts, struct psc_mfields *mres)
{
  rho_1st_nc_cuda_run_patches<1, 4, 4, IP_EC, DEPOSIT_VB_3D, CURRMEM_GLOBAL>(mprts, mres);
#if 0
  // FIXME, make sure is reordered -- or handle if not
  for (int p = 0; p < mres->nr_patches; p++) {
    do_rho_run(p, psc_mfields_get_patch(mres, p), psc_mparticles_get_patch(mprts, p));
  }
#endif
}

