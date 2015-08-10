
#include "psc_debug.h"

// ======================================================================
// EXT_PREPARE_SORT
//
// if enabled, calculate the new block position for each particle once
// it's moved, and also append particles that left the block to an extra
// list at the end of all local particles (hopefully there's enough room...)

#ifdef EXT_PREPARE_SORT

#ifdef PSC_PARTICLES_AS_SINGLE
static struct psc_particles_single *prts_sub;
#pragma omp threadprivate(prts_sub)
#endif

static inline void
ext_prepare_sort_before(struct psc_particles *prts)
{
  prts_sub = psc_particles_single(prts);
  memset(prts_sub->b_cnt, 0,
	 (prts_sub->nr_blocks + 1) * sizeof(*prts_sub->b_cnt));
}

static inline void
ext_prepare_sort(struct psc_particles *prts, int n, particle_t *prt,
		 int *b_pos)
{
  /* FIXME, only if blocksize == 1! */
  int *b_mx = prts_sub->b_mx;
  if (b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
      b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
    prts_sub->b_idx[n] = b_pos[2] * b_mx[1] + b_pos[1];
  } else { /* out of bounds */
    prts_sub->b_idx[n] = prts_sub->nr_blocks;
    assert(prts_sub->b_cnt[prts_sub->nr_blocks] < prts_sub->n_alloced);
    /* append to back */
    *particles_get_one(prts, prts->n_part + prts_sub->b_cnt[prts_sub->nr_blocks]) = *prt;
  }
  prts_sub->b_cnt[prts_sub->b_idx[n]]++;
}

#else

static inline void
ext_prepare_sort_before(struct psc_particles *prts)
{
}

static inline void
ext_prepare_sort(struct psc_particles *prts, int n, particle_t *prt,
		 int *b_pos)
{
}

#endif

// ======================================================================

// ----------------------------------------------------------------------
// push_one

#ifdef __CUDACC__

__device__ static void
push_one(particle_t *prt, int n, float4 *d_xi4, float4 *d_pxi4,
	 real *flds_em, flds_curr_t flds_curr, int ci0[3])
{
  PARTICLE_CUDA2_LOAD_POS(*prt, d_xi4, n);

  // here we have x^{n+.5}, p^n

  // field interpolation
  real exq, eyq, ezq, hxq, hyq, hzq;
  int lg[3];
  real og[3];
  find_idx_off_1st_rel(prt->xi, lg, og, real(0.));
  INTERPOLATE_1ST_EC(flds_em, exq, eyq, ezq, hxq, hyq, hzq);

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  PARTICLE_CUDA2_LOAD_MOM(*prt, d_pxi4, n);
  int kind = particle_kind(prt);
  real dq = prm.dq_kind[kind];
  push_pxi(prt, exq, eyq, ezq, hxq, hyq, hzq, dq);
  PARTICLE_CUDA2_STORE_MOM(*prt, d_pxi4, n);

  real vxi[3];
  calc_vxi(vxi, prt);

  particle_real_t xm[3], xp[3];
  int lf[3];

  // position xm at x^(n+.5)
  real h0[3];
  find_idx_off_pos_1st_rel(prt->xi, lg, h0, xm, real(0.));

  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
  push_xi(prt, vxi, prm.dt);
  PARTICLE_CUDA2_STORE_POS(*prt, d_xi4, n);

  // position xp at x^(n+.5)
  real h1[3];
  find_idx_off_pos_1st_rel(prt->xi, lf, h1, xp, real(0.));

  calc_j(flds_curr, xm, xp, lf, lg, prt, vxi);
}

__device__ static void
push_one_a(particle_t *prt, int n, float4 *d_xi4, float4 *d_pxi4,
	   real *flds_em, flds_curr_t flds_curr, int ci0[3])
{
  PARTICLE_CUDA2_LOAD_POS(*prt, d_xi4, n);

  // here we have x^{n+.5}, p^n

  // field interpolation
  real exq, eyq, ezq, hxq, hyq, hzq;
  int lg[3];
  real og[3];
  find_idx_off_1st_rel(prt->xi, lg, og, real(0.));

  INTERPOLATE_1ST_EC(flds_em, exq, eyq, ezq, hxq, hyq, hzq);

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  PARTICLE_CUDA2_LOAD_MOM(*prt, d_pxi4, n);
  int kind = particle_kind(prt);
  real dq = prm.dq_kind[kind];
  push_pxi(prt, exq, eyq, ezq, hxq, hyq, hzq, dq);
  PARTICLE_CUDA2_STORE_MOM(*prt, d_pxi4, n);
}

__device__ static void
push_one_b(particle_t *prt, int n, float4 *d_xi4, float4 *d_pxi4,
	   real *flds_em, flds_curr_t flds_curr, int ci0[3])
{
  PARTICLE_CUDA2_LOAD_POS(*prt, d_xi4, n);
  PARTICLE_CUDA2_LOAD_MOM(*prt, d_pxi4, n);

  real vxi[3];
  calc_vxi(vxi, prt);

  particle_real_t xm[3], xp[3];
  int lg[3], lf[3];

  // position xm at x^(n+.5)
  real h0[3];
  find_idx_off_pos_1st_rel(prt->xi, lg, h0, xm, real(0.));

  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
  push_xi(prt, vxi, prm.dt);
  PARTICLE_CUDA2_STORE_POS(*prt, d_xi4, n);

  // position xp at x^(n+.5)
  real h1[3];
  find_idx_off_pos_1st_rel(prt->xi, lf, h1, xp, real(0.));

  calc_j(flds_curr, xm, xp, lf, lg, prt, vxi);
}

#else

static inline void
push_one(particle_t *prt, struct psc_fields *flds, struct psc_particles *prts, int n)
{
  // field interpolation
  int lg[3], lh[3];
  particle_real_t og[3], oh[3], xm[3];
  find_idx_off_pos_1st_rel(&particle_x(prt), lg, og, xm, 0.f); // FIXME passing xi hack
  find_idx_off_1st_rel(&particle_x(prt), lh, oh, -.5f);

  particle_real_t exq, eyq, ezq, hxq, hyq, hzq;
  INTERPOLATE_1ST(flds, exq, eyq, ezq, hxq, hyq, hzq);
  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
  int kind = particle_kind(prt);
  particle_real_t dq = prm.dq_kind[kind];
  push_pxi(prt, exq, eyq, ezq, hxq, hyq, hzq, dq);

  particle_real_t vxi[3];
  calc_vxi(vxi, prt);
#if CALC_J == CALC_J_1VB_2D
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0)
  push_xi(prt, vxi, .5f * prm.dt);
  
  // OUT OF PLANE CURRENT DENSITY AT (n+1.0)*dt
  calc_j_oop(flds, prt, vxi);
  
  // x^(n+1), p^(n+1) -> x^(n+1.5), p^(n+1)
  push_xi(prt, vxi, .5f * prm.dt);
#else
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
  push_xi(prt, vxi, prm.dt);
#endif
  
  int lf[3];
  particle_real_t of[3], xp[3];
  find_idx_off_pos_1st_rel(&particle_x(prt), lf, of, xp, 0.f);
  ext_prepare_sort(prts, n, prt, lf);

  // CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
  calc_j(flds, xm, xp, lf, lg, prt, vxi);
}

#endif

// ----------------------------------------------------------------------
// push_one_mprts

#if PSC_PARTICLES_AS_CUDA2

#ifdef __CUDACC__

CUDA_DEVICE static void
push_one_mprts(float4 *d_xi4, float4 *d_pxi4, int n,
	       real *flds_em, flds_curr_t flds_curr, int ci0[3])
{
  particle_t prt;

  push_one(&prt, n, d_xi4, d_pxi4, flds_em, flds_curr, ci0);
}

CUDA_DEVICE static void
push_one_mprts_a(float4 *d_xi4, float4 *d_pxi4, int n,
		 real *flds_em, flds_curr_t flds_curr, int ci0[3])
{
  particle_t prt;

  push_one_a(&prt, n, d_xi4, d_pxi4, flds_em, flds_curr, ci0);
}

CUDA_DEVICE static void
push_one_mprts_b(float4 *d_xi4, float4 *d_pxi4, int n,
		 real *flds_em, flds_curr_t flds_curr, int ci0[3])
{
  particle_t prt;

  push_one_b(&prt, n, d_xi4, d_pxi4, flds_em, flds_curr, ci0);
}

#else

static inline void
push_one_mprts(struct psc_mparticles *mprts, struct psc_mfields *mflds, int n, int p)
{
  struct psc_mparticles_cuda2 *mprts_sub = psc_mparticles_cuda2(mprts);

  struct psc_fields *flds = psc_mfields_get_patch(mflds, p);

  particle_t prt;
  PARTICLE_CUDA2_LOAD_POS(prt, mprts_sub->h_xi4, n);
  PARTICLE_CUDA2_LOAD_MOM(prt, mprts_sub->h_pxi4, n);
  
  push_one(&prt, flds, NULL, n);

  PARTICLE_CUDA2_STORE_POS(prt, mprts_sub->h_xi4, n);
  PARTICLE_CUDA2_STORE_MOM(prt, mprts_sub->h_pxi4, n);
}

#endif

#else

static void
push_one_prts(struct psc_particles *prts, struct psc_fields *flds, int n)
{
  particle_t *prt = particles_get_one(prts, n);

  push_one(prt, flds, prts, n);
}

#endif
