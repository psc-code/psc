
#include "psc_debug.h"

#ifndef SFX

#if DIM == DIM_1
#define SFX(s) s ## _1
#elif DIM == DIM_YZ
#define SFX(s) s ## _yz
#elif DIM == DIM_XYZ
#define SFX(s) s ## _xyz
#endif

#endif // !SFX

// ----------------------------------------------------------------------
// the following macros are convoluted, as they handle two quite different
// cases for now:
// - cuda2/acc particles provide the mprts-based array, so n
//   includes the patch offset
// - otherwise, we provide per-patch particles, so n is the patch-local
//   particle number

#if PSC_PARTICLES_AS_CUDA2 || PSC_PARTICLES_AS_ACC

typedef struct { float4 *xi4; float4 *pxi4; } mprts_array_t;

#if PSC_PARTICLES_AS_CUDA2

#define PARTICLE_LOAD(prt, mprts_arr, n)		\
  particle_t _prt;					\
  prt = &_prt;						\
  PARTICLE_CUDA2_LOAD_POS(*prt, mprts_arr.xi4, n);	\
  PARTICLE_CUDA2_LOAD_MOM(*prt, mprts_arr.pxi4, n)

#define PARTICLE_STORE(prt, mprts_arr, n)		\
  PARTICLE_CUDA2_STORE_POS(*prt, mprts_arr.xi4, n);	\
  PARTICLE_CUDA2_STORE_MOM(*prt, mprts_arr.pxi4, n)	\

#elif PSC_PARTICLES_AS_ACC

#define PARTICLE_LOAD(prt, mprts_arr, n)		\
  particle_t _prt;					\
  prt = &_prt;						\
  PARTICLE_ACC_LOAD_POS(*prt, mprts_arr.xi4, n);	\
  PARTICLE_ACC_LOAD_MOM(*prt, mprts_arr.pxi4, n)

#define PARTICLE_STORE(prt, mprts_arr, n)		\
  PARTICLE_ACC_STORE_POS(*prt, mprts_arr.xi4, n);	\
  PARTICLE_ACC_STORE_MOM(*prt, mprts_arr.pxi4, n)

#endif

#else

typedef particle_iter_t mprts_array_t;

#define PARTICLE_LOAD(prt, mprts_arr, n)	\
  prt = particle_iter_at(mprts_arr, n)

#define PARTICLE_STORE(prt, mprts_arr, n) do {} while (0)

#endif

// ======================================================================
// EXT_PREPARE_SORT
//
// if enabled, calculate the new block position for each particle once
// it's moved, and also append particles that left the block to an extra
// list at the end of all local particles (hopefully there's enough room...)

#ifdef EXT_PREPARE_SORT

static inline void
ext_prepare_sort_before(struct psc_mparticles *mprts, int p)
{
  struct psc_mparticles_single *sub = psc_mparticles_single(mprts);
  struct psc_mparticles_single_patch *patch = &sub->patch[p];
  
  memset(patch->b_cnt, 0, (patch->nr_blocks + 1) * sizeof(*patch->b_cnt));
}

static inline void
ext_prepare_sort(struct psc_mparticles *mprts, int p, int n, particle_t *prt,
		 int *b_pos)
{
  struct psc_mparticles_single *sub = psc_mparticles_single(mprts);
  struct psc_mparticles_single_patch *patch = &sub->patch[p];
  particle_range_t prts = particle_range_mprts(mprts, p);
  unsigned int n_prts = particle_range_size(prts);
  /* FIXME, only if blocksize == 1! */
  int *b_mx = patch->b_mx;
  if (b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
      b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
    patch->b_idx[n] = b_pos[2] * b_mx[1] + b_pos[1];
  } else { /* out of bounds */
    patch->b_idx[n] = patch->nr_blocks;
    /* append to back */
    *particle_iter_at(prts.begin, n_prts + patch->b_cnt[patch->nr_blocks]) = *prt;
  }
  patch->b_cnt[patch->b_idx[n]]++;
}

#else

static inline void
ext_prepare_sort_before(struct psc_mparticles *mprts, int p)
{
}

static inline void
ext_prepare_sort(struct psc_mparticles *mprts, int p, int n, particle_t *prt,
		 int *b_pos)
{
}

#endif

// ======================================================================

// ----------------------------------------------------------------------
// push_one
//
// as opposed to what the name implies, mprts_arr may actually be the
// per-patch particles, in which case n needs to be the patch-local
// patch number (see also above)
//
// the cuda2 version is very similar to the generic one, except
// - it does spread the particles loads/stores as necessary
// - it doesn't handle 1VB_2D

CUDA_DEVICE static void
push_one(mprts_array_t mprts_arr, int n,
	 em_cache_t flds_em, curr_cache_t curr_cache)
{
  FieldsEM EM(flds_em);
#if PSC_PARTICLES_AS_CUDA2
  particle_t _prt, *prt = &_prt;
  PARTICLE_CUDA2_LOAD_POS(*prt, mprts_arr.xi4, n);

  // here we have x^{n+.5}, p^n

  // field interpolation
  int lg[3];
  real og[3];
  find_idx_off_1st_rel(prt->xi, lg, og, real(0.));

  INTERPOLATE_FIELDS(flds_em);

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0) 
  PARTICLE_CUDA2_LOAD_MOM(*prt, mprts_arr.pxi4, n);
  int kind = particle_kind(prt);
  real dq = prm.dq_kind[kind];
  push_pxi(prt, exq, eyq, ezq, hxq, hyq, hzq, dq);
  PARTICLE_CUDA2_STORE_MOM(*prt, mprts_arr.pxi4, n);

  real vxi[3];
  calc_vxi(vxi, prt);

  particle_real_t xm[3], xp[3];
  int lf[3];

  // position xm at x^(n+.5)
  real h0[3];
  find_idx_off_pos_1st_rel(prt->xi, lg, h0, xm, real(0.));

  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0) 
  push_xi(prt, vxi, prm.dt);
  PARTICLE_CUDA2_STORE_POS(*prt, mprts_arr.xi4, n);

  // position xp at x^(n+.5)
  real h1[3];
  find_idx_off_pos_1st_rel(prt->xi, lf, h1, xp, real(0.));

  calc_j(curr_cache, xm, xp, lf, lg, prt, vxi);

#else

  particle_t *prt;
  PARTICLE_LOAD(prt, mprts_arr, n);
  
  // field interpolation
  particle_real_t *xi = &particle_x(prt);

  particle_real_t xm[3];
  for (int d = 0; d < 3; d++) {
    xm[d] = xi[d] * c_prm.dxi[d];
  }

  // FIELD INTERPOLATION

  IP<NNN> ip;
  INTERPOLATE_FIELDS(flds_em);

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
  int kind = particle_kind(prt);
  particle_real_t dq = prm.dq_kind[kind];
  push_p(&particle_px(prt), ip.E, ip.H, dq);

  particle_real_t vxi[3];
  calc_v(vxi, &particle_px(prt));
#if CALC_J == CALC_J_1VB_2D
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0)
  push_x(&prt->xi, vxi, .5f * c_prm.dt);
  
  // OUT OF PLANE CURRENT DENSITY AT (n+1.0)*dt
  calc_j_oop(curr_cache, prt, vxi);
  
  // x^(n+1), p^(n+1) -> x^(n+1.5), p^(n+1)
  push_x(&prt->xi, vxi, .5f * c_prm.dt);
#else
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
  push_x(&prt->xi, vxi, c_prm.dt);
#endif
  
  int lf[3];
  particle_real_t of[3], xp[3];
  find_idx_off_pos_1st_rel(&particle_x(prt), lf, of, xp, 0.f);
  //  ext_prepare_sort(prts, n, prt, lf);

  // CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
  int lg[3];
  IF_DIM_X( lg[0] = ip.lg1; );
  IF_DIM_Y( lg[1] = ip.lg2; );
  IF_DIM_Z( lg[2] = ip.lg3; );
  calc_j(curr_cache, xm, xp, lf, lg, prt, vxi);

#if !(PUSH_DIM & DIM_X)
  prt->xi = 0.f;
#endif
#if !(PUSH_DIM & DIM_Y)
  prt->yi = 0.f;
#endif
#if !(PUSH_DIM & DIM_Z)
  prt->zi = 0.f;
#endif

  PARTICLE_STORE(prt, mprts_arr, n);
#endif
}

// ----------------------------------------------------------------------
// stagger_one

CUDA_DEVICE static void
stagger_one(mprts_array_t mprts_arr, int n,
	    em_cache_t flds_em)
{
  FieldsEM EM(flds_em);
  particle_t *prt;
  PARTICLE_LOAD(prt, mprts_arr, n);
  
  // field interpolation
  particle_real_t *xi = &particle_x(prt);

  particle_real_t xm[3];
  for (int d = 0; d < 3; d++) {
    xm[d] = xi[d] * c_prm.dxi[d];
  }

  // FIELD INTERPOLATION

  IP<NNN> ip;
  INTERPOLATE_FIELDS(flds_em);

  // x^(n+1/2), p^{n+1/2} -> x^(n+1/2), p^{n}
  int kind = particle_kind(prt);
  particle_real_t dq = prm.dq_kind[kind];
  push_p(&particle_px(prt), ip.E, ip.H, -.5f * dq);

  PARTICLE_STORE(prt, mprts_arr, n);
}
