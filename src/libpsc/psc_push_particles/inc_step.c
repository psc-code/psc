
#include "psc_debug.h"

#define PARTICLE_LOAD(prt, mprts_arr, n)	\
  prt = &mprts_arr[n]

#define PARTICLE_STORE(prt, mprts_arr, n) do {} while (0)

// ======================================================================
// EXT_PREPARE_SORT
//
// if enabled, calculate the new block position for each particle once
// it's moved, and also append particles that left the block to an extra
// list at the end of all local particles (hopefully there's enough room...)

template<typename mparticles_t, typename OPT_EXT>
struct ExtPrepareSort
{
  using particle_t = typename mparticles_t::particle_t;

  static void before(typename mparticles_t::patch_t& prtss)
  {}

  static void sort(typename mparticles_t::patch_t prts, int n, particle_t *prt, int *b_pos)
  {}
};

template<typename mparticles_t>
struct ExtPrepareSort<mparticles_t, opt_ext_prepare_sort>
{
  using particle_t = typename mparticles_t::particle_t;
  
  static void before(typename mparticles_t::patch_t& prts)
  {
    memset(prts.b_cnt, 0, (prts.nr_blocks + 1) * sizeof(*prts.b_cnt));
  }

  static void sort(typename mparticles_t::patch_t& prts, int n, particle_t *prt,
		   int *b_pos)
  {
    unsigned int n_prts = prts.size();
    /* FIXME, only if blocksize == 1! */
    int *b_mx = prts.pi_.b_mx_;
    if (b_pos[1] >= 0 && b_pos[1] < b_mx[1] &&
	b_pos[2] >= 0 && b_pos[2] < b_mx[2]) {
      prts.b_idx[n] = b_pos[2] * b_mx[1] + b_pos[1];
    } else { /* out of bounds */
      prts.b_idx[n] = prts.nr_blocks;
      /* append to back */
      prts[n_prts + prts.b_cnt[prts.nr_blocks]] = *prt;
    }
    prts.b_cnt[prts.b_idx[n]]++;
  }
};

// ======================================================================

// ----------------------------------------------------------------------
// push_one

template<typename C, typename particles_t>
CUDA_DEVICE static void
push_one(particles_t& prts, int n,
	 typename C::Mfields::fields_t flds_em, curr_cache_t curr_cache)
{
  using dim = typename C::dim;
  using FieldsEM = typename C::FieldsEM;
  using IP = InterpolateEM<FieldsEM, typename C::ip, dim>;
  using AdvanceParticle_t = AdvanceParticle<real_t, dim>;

  AdvanceParticle_t advance(prts.grid().dt);
  FieldsEM EM(flds_em);

  particle_t *prt;
  PARTICLE_LOAD(prt, prts, n);
  
  // field interpolation
  real_t *xi = &prt->xi;

  real_t xm[3];
  for (int d = 0; d < 3; d++) {
    xm[d] = xi[d] * c_prm.dxi[d];
  }

  // FIELD INTERPOLATION

  IP ip;
  ip.set_coeffs(xm);
  // FIXME, we're not using EM instead flds_em
  real_t E[3] = { ip.ex(flds_em), ip.ey(flds_em), ip.ez(flds_em) };
  real_t H[3] = { ip.hx(flds_em), ip.hy(flds_em), ip.hz(flds_em) };

  // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
  int kind = prt->kind();
  real_t dq = prm.dq_kind[kind];
  advance.push_p(&prt->pxi, E, H, dq);

  real_t vxi[3];
  advance.calc_v(vxi, &prt->pxi);
#if CALC_J == CALC_J_1VB_2D
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0)
  advance.push_x(&prt->xi, vxi, .5f);
  
  // OUT OF PLANE CURRENT DENSITY AT (n+1.0)*dt
  calc_j_oop(curr_cache, prt, vxi);
  
  // x^(n+1), p^(n+1) -> x^(n+1.5), p^(n+1)
  advance.push_x(&prt->xi, vxi, .5f);
#else
  // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
  advance.push_x(&prt->xi, vxi);
#endif
  
  int lf[3];
  real_t of[3], xp[3];
  find_idx_off_pos_1st_rel(&prt->xi, lf, of, xp, real_t(0.));
  //  ext_prepare_sort(prts, n, prt, lf);

  // CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
  int lg[3];
  IF_DIM_X( lg[0] = ip.cx.g.l; );
  IF_DIM_Y( lg[1] = ip.cy.g.l; );
  IF_DIM_Z( lg[2] = ip.cz.g.l; );
  calc_j(curr_cache, xm, xp, lf, lg, prt, vxi);

#ifdef PUSH_DIM
#if !(PUSH_DIM & DIM_X)
  prt->xi = 0.f;
#endif
#if !(PUSH_DIM & DIM_Y)
  prt->yi = 0.f;
#endif
#if !(PUSH_DIM & DIM_Z)
  prt->zi = 0.f;
#endif
#endif

  PARTICLE_STORE(prt, mprts_arr, n);
}

// ----------------------------------------------------------------------
// stagger_one

template<typename C>
CUDA_DEVICE static void
stagger_one(mparticles_t::patch_t& prts, int n, typename C::Mfields::fields_t flds_em)
{
  using dim = typename C::dim;
  using FieldsEM = typename C::FieldsEM;
  using IP = InterpolateEM<FieldsEM, typename C::ip, dim>;
  using AdvanceParticle_t = AdvanceParticle<real_t, dim>;

  AdvanceParticle_t advance(prts.grid().dt);

  FieldsEM EM(flds_em);
  particle_t *prt;
  PARTICLE_LOAD(prt, prts, n);
  
  // field interpolation
  real_t *xi = &prt->xi;

  real_t xm[3];
  for (int d = 0; d < 3; d++) {
    xm[d] = xi[d] * c_prm.dxi[d];
  }

  // FIELD INTERPOLATION

  IP ip;
  ip.set_coeffs(xm);
  // FIXME, we're not using EM instead flds_em
  real_t E[3] = { ip.ex(flds_em), ip.ey(flds_em), ip.ez(flds_em) };
  real_t H[3] = { ip.hx(flds_em), ip.hy(flds_em), ip.hz(flds_em) };

  // x^(n+1/2), p^{n+1/2} -> x^(n+1/2), p^{n}
  int kind = prt->kind();
  real_t dq = prm.dq_kind[kind];
  advance.push_p(&prt->pxi, E, H, -.5f * dq);

  PARTICLE_STORE(prt, mprts_arr, n);
}
