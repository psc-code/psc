
#include "cuda_mparticles.h"
#include "cuda_mfields.h"

#include "cuda_mparticles_const.h"
#include "cuda_mfields_const.h"

#include "psc.h" // FIXME

#define THREADS_PER_BLOCK (512)

// ======================================================================
// field caching

#define F3_CACHE(fld_cache, m, jy, jz)					\
  ((fld_cache)[(((m-EX)							\
		 *(BLOCKSIZE_Z + 4) + ((jz)-(-2)))			\
		*(BLOCKSIZE_Y + 4) + ((jy)-(-2)))])

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static void
cache_fields(float *fld_cache, DMFields d_flds0, int *ci0, int p)
{
  DFields d_flds = d_flds0[p];

  int ti = threadIdx.x;
  int n = BLOCKSIZE_X * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4);
  while (ti < n) {
    int tmp = ti;
    int jy = tmp % (BLOCKSIZE_Y + 4) - 2;
    tmp /= BLOCKSIZE_Y + 4;
    int jz = tmp % (BLOCKSIZE_Z + 4) - 2;
    // OPT? currently it seems faster to do the loop rather than do m by threadidx
    for (int m = EX; m <= HZ; m++) {
      F3_CACHE(fld_cache, m, jy, jz) = d_flds(m, 0,jy+ci0[1],jz+ci0[2]);
    }
    ti += THREADS_PER_BLOCK;
  }
}

// ----------------------------------------------------------------------
// push_xi
//
// advance position using velocity

__device__ static void
push_xi(struct d_particle *p, const float vxi[3], float dt)
{
  int d;
  for (d = 1; d < 3; d++) {
    p->xi[d] += dt * vxi[d];
  }
}

// ----------------------------------------------------------------------
// calc_vxi
//
// calculate velocity from momentum

__device__ static void
calc_vxi(float vxi[3], struct d_particle p)
{
  float root = rsqrtf(float(1.) + sqr(p.pxi[0]) + sqr(p.pxi[1]) + sqr(p.pxi[2]));

  int d;
  for (d = 0; d < 3; d++) {
    vxi[d] = p.pxi[d] * root;
  }
}

// ----------------------------------------------------------------------
// push_pxi_dt
//
// advance moments according to EM fields

__device__ static void
push_pxi_dt(struct d_particle *p,
	    float exq, float eyq, float ezq, float hxq, float hyq, float hzq)
{
  int kind = __float_as_int(p->kind_as_float);
  float dq = d_cmprts_const.dq[kind];
  float pxm = p->pxi[0] + dq*exq;
  float pym = p->pxi[1] + dq*eyq;
  float pzm = p->pxi[2] + dq*ezq;
  
  float root = dq * rsqrtf(float(1.) + sqr(pxm) + sqr(pym) + sqr(pzm));
  float taux = hxq * root, tauy = hyq * root, tauz = hzq * root;
  
  float tau = float(1.) / (float(1.) + sqr(taux) + sqr(tauy) + sqr(tauz));
  float pxp = ( (float(1.) + sqr(taux) - sqr(tauy) - sqr(tauz)) * pxm
	       +(float(2.)*taux*tauy + float(2.)*tauz)*pym
	       +(float(2.)*taux*tauz - float(2.)*tauy)*pzm)*tau;
  float pyp = ( (float(2.)*taux*tauy - float(2.)*tauz)*pxm
	       +(float(1.) - sqr(taux) + sqr(tauy) - sqr(tauz)) * pym
	       +(float(2.)*tauy*tauz + float(2.)*taux)*pzm)*tau;
  float pzp = ( (float(2.)*taux*tauz + float(2.)*tauy)*pxm
	       +(float(2.)*tauy*tauz - float(2.)*taux)*pym
	       +(float(1.) - sqr(taux) - sqr(tauy) + sqr(tauz))*pzm)*tau;
  
  p->pxi[0] = pxp + dq * exq;
  p->pxi[1] = pyp + dq * eyq;
  p->pxi[2] = pzp + dq * ezq;
}

// ----------------------------------------------------------------------
// push_part_one

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__device__ static void
push_part_one(struct d_particle *prt, int n, uint *d_ids, float4 *d_xi4, float4 *d_pxi4,
	      float4 *d_alt_xi4, float4 *d_alt_pxi4,
	      float *fld_cache, int ci0[3])
{
  LOAD_PARTICLE_POS(*prt, d_xi4, n);

  // here we have x^{n+.5}, p^n

  // field interpolation
  float exq, eyq, ezq, hxq, hyq, hzq;
  int lg[3];
  float og[3];
  find_idx_off_1st(prt->xi, lg, og, float(0.));
  lg[1] -= ci0[1];
  lg[2] -= ci0[2];
  
#if 0
  if (lg[1] < -2 || lg[1] >= BLOCKSIZE_Y + 1) {
    printf("lg[1] %d\n", lg[1]);
  }
  if (lg[2] < -2 || lg[2] >= BLOCKSIZE_Z + 1) {
    printf("lg[2] %d\n", lg[2]);
  }
#endif
  exq = ((1.f - og[1]) * (1.f - og[2]) * F3_CACHE(fld_cache, EX, lg[1]+0, lg[2]+0) +
	 (      og[1]) * (1.f - og[2]) * F3_CACHE(fld_cache, EX, lg[1]+1, lg[2]+0) +
	 (1.f - og[1]) * (      og[2]) * F3_CACHE(fld_cache, EX, lg[1]+0, lg[2]+1) +
	 (      og[1]) * (      og[2]) * F3_CACHE(fld_cache, EX, lg[1]+1, lg[2]+1));
  eyq = ((1.f - og[2]) * F3_CACHE(fld_cache, EY, lg[1]  , lg[2]+0) +
	 (      og[2]) * F3_CACHE(fld_cache, EY, lg[1]  , lg[2]+1));
  ezq = ((1.f - og[1]) * F3_CACHE(fld_cache, EZ, lg[1]+0, lg[2]  ) +
	 (      og[1]) * F3_CACHE(fld_cache, EZ, lg[1]+1, lg[2]  ));
  hxq = (F3_CACHE(fld_cache, HX, lg[1]  , lg[2]  ));
  hyq = ((1.f - og[1]) * F3_CACHE(fld_cache, HY, lg[1]+0, lg[2]  ) +
	 (      og[1]) * F3_CACHE(fld_cache, HY, lg[1]+1, lg[2]  ));
  hzq = ((1.f - og[2]) * F3_CACHE(fld_cache, HZ, lg[1]  , lg[2]+0) +
	 (      og[2]) * F3_CACHE(fld_cache, HZ, lg[1]  , lg[2]+1));

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  LOAD_PARTICLE_MOM(*prt, d_pxi4, n);
  push_pxi_dt(prt, exq, eyq, ezq, hxq, hyq, hzq);
  STORE_PARTICLE_MOM(*prt, d_pxi4, n);

  float vxi[3];
  calc_vxi(vxi, *prt);

  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
  push_xi(prt, vxi, d_cmprts_const.dt);
  STORE_PARTICLE_POS(*prt, d_xi4, n);
}

// ----------------------------------------------------------------------
// push_mprts_ab

template<int BLOCKSIZE_X, int BLOCKSIZE_Y, int BLOCKSIZE_Z>
__global__ static void
push_mprts_ab(int block_start, float4 *d_xi4, float4 *d_pxi4,
	      float4 *d_alt_xi4, float4 *d_alt_pxi4,
	      uint *d_off, int nr_total_blocks, uint *d_ids, uint *d_bidx,
	      DMFields d_flds0, uint size)
{
  int block_pos[3], ci0[3];						\
  int p, bid;								\
  p = find_block_pos_patch_q<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>	\
    (block_pos, ci0, block_start);					\
  if (p < 0)								\
    return;								\
    									\
  bid = find_bid_q(p, block_pos);					\
  int block_begin = d_off[bid];						\
  int block_end = d_off[bid + 1];					\
    
  __shared__ float fld_cache[6 * 1 * (BLOCKSIZE_Y + 4) * (BLOCKSIZE_Z + 4)];
  cache_fields<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
    (fld_cache, d_flds0, ci0, p);
    
  __syncthreads();
  for (int n = (block_begin & ~31) + threadIdx.x; n < block_end; n += THREADS_PER_BLOCK) {
    if (n < block_begin) {
      continue;
    }
    struct d_particle prt;
    push_part_one<BLOCKSIZE_X, BLOCKSIZE_Y, BLOCKSIZE_Z>
      (&prt, n, d_ids, d_xi4, d_pxi4, d_alt_xi4, d_alt_pxi4, fld_cache, ci0);
  }
}

// ----------------------------------------------------------------------
// cuda_push_mprts_ab

static void
cuda_push_mprts_ab(struct cuda_mparticles *cmprts, struct cuda_mfields *cmflds)
{
  cuda_mfields_const_set(cmflds);
  cuda_mparticles_const_set(cmprts);

  uint fld_size = cmflds->n_fields * cmflds->n_cells_per_patch;

  int gx, gy;
  gx = cmprts->b_mx_[1];
  gy = cmprts->b_mx_[2] * cmprts->n_patches;
  dim3 dimGrid(gx, gy);

  for (int block_start = 0; block_start < 4; block_start++) {
    push_mprts_ab<1, 1, 1>
      <<<dimGrid, THREADS_PER_BLOCK>>>
      (block_start, cmprts->d_xi4.data().get(), cmprts->d_pxi4.data().get(),
       cmprts->d_alt_xi4.data().get(), cmprts->d_alt_pxi4.data().get(), cmprts->d_off.data().get(),
       cmprts->n_blocks, cmprts->d_id.data().get(), cmprts->d_bidx.data().get(),
       DMFields(cmflds), fld_size);
    cuda_sync_if_enabled();
  }
}

// ----------------------------------------------------------------------
// xyz_cuda_push_mprts

static void
xyz_cuda_push_mprts(struct cuda_mparticles *cmprts, struct cuda_mfields *cmflds)
{
  assert(!cmprts->need_reorder);
  cuda_push_mprts_ab(cmprts, cmflds);
}

// ----------------------------------------------------------------------
// cuda_push_mprts_xyz

void cuda_push_mprts_xyz(cuda_mparticles *cmprts, cuda_mfields *cmflds)
{
  return xyz_cuda_push_mprts(cmprts, cmflds);
}

