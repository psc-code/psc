
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "1vb/psc_push_particles_1vb.h"

#include "fields.hxx"

#include "interpolate.hxx"
#include "inc_push.c"
#include "inc_curr.c"

#define MAX_NR_KINDS (10)

// ======================================================================
// PushParticles1vb

template<typename C>
struct PushParticles1vb
{
  using Mparticles = typename C::Mparticles;
  using Mfields = typename C::Mfields;
  using mparticles_t = PscMparticles<Mparticles>;
  using real_t = typename Mparticles::real_t;
  using Real3 = Vec3<real_t>;
  using dim = typename C::dim;
  using FieldsEM = typename C::FieldsEM;
  using IP = InterpolateEM<FieldsEM, typename C::ip, dim>;
  using AdvanceParticle_t = AdvanceParticle<real_t, dim>;
  using Current = typename C::Current_t;
  
  // ----------------------------------------------------------------------
  // push

  static void push(typename C::Mparticles::patch_t& prts, int n,
		   typename C::Mfields::fields_t flds_em, typename C::curr_cache_t curr_cache)
  {
    AdvanceParticle_t advance(prts.grid().dt);
    FieldsEM EM(flds_em);
    Current current(prts.grid());
    PI<real_t> pi(prts.grid());
    Real3 dxi = Real3{ 1., 1., 1. } / Real3(prts.grid().dx);
    real_t dq_kind[MAX_NR_KINDS];
    auto& kinds = prts.grid().kinds;
    assert(kinds.size() <= MAX_NR_KINDS);
    for (int k = 0; k < kinds.size(); k++) {
      dq_kind[k] = .5f * prts.grid().eta * prts.grid().dt * kinds[k].q / kinds[k].m;
    }

    particle_t *prt = &prts[n];
  
    // field interpolation
    real_t *xi = &prt->xi;

    real_t xm[3];
    for (int d = 0; d < 3; d++) {
      xm[d] = xi[d] * dxi[d];
    }

    // FIELD INTERPOLATION

    IP ip;
    ip.set_coeffs(xm);
    // FIXME, we're not using EM instead flds_em
    real_t E[3] = { ip.ex(flds_em), ip.ey(flds_em), ip.ez(flds_em) };
    real_t H[3] = { ip.hx(flds_em), ip.hy(flds_em), ip.hz(flds_em) };

    // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
    int kind = prt->kind();
    real_t dq = dq_kind[kind];
    advance.push_p(&prt->pxi, E, H, dq);

    real_t vxi[3];
    advance.calc_v(vxi, &prt->pxi);

    int lf[3];
    real_t of[3], xp[3];
#if CALC_J == CALC_J_1VB_2D
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0)
    advance.push_x(&prt->xi, vxi, .5f);
  
    // OUT OF PLANE CURRENT DENSITY AT (n+1.0)*dt
    pi.find_idx_off_1st_rel(&prt->xi, lf, of, real_t(0.));
    current.calc_j_oop(curr_cache, particle_qni_wni(prt), vxi, lf, of);
  
    // x^(n+1), p^(n+1) -> x^(n+1.5), p^(n+1)
    advance.push_x(&prt->xi, vxi, .5f);
#else
    // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
    advance.push_x(&prt->xi, vxi);
#endif
  
    pi.find_idx_off_pos_1st_rel(&prt->xi, lf, of, xp, real_t(0.));

    // CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
    int lg[3];
    if (!dim::InvarX::value) { lg[0] = ip.cx.g.l; }
    if (!dim::InvarY::value) { lg[1] = ip.cy.g.l; }
    if (!dim::InvarZ::value) { lg[2] = ip.cz.g.l; }
    current.calc_j(curr_cache, xm, xp, lf, lg, particle_qni_wni(prt), vxi);

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
  }

  // ----------------------------------------------------------------------
  // stagger

  static void stagger(typename mparticles_t::patch_t& prts, int n, typename C::Mfields::fields_t flds_em)
  {
    AdvanceParticle_t advance(prts.grid().dt);
    FieldsEM EM(flds_em);
    Real3 dxi = Real3{ 1., 1., 1. } / Real3(prts.grid().dx);
    real_t dq_kind[MAX_NR_KINDS];
    auto& kinds = prts.grid().kinds;
    assert(kinds.size() <= MAX_NR_KINDS);
    for (int k = 0; k < kinds.size(); k++) {
      dq_kind[k] = .5f * prts.grid().eta * prts.grid().dt * kinds[k].q / kinds[k].m;
    }

    particle_t *prt = &prts[n];
  
    // field interpolation
    real_t *xi = &prt->xi;

    real_t xm[3];
    for (int d = 0; d < 3; d++) {
      xm[d] = xi[d] * dxi[d];
    }

    // FIELD INTERPOLATION

    IP ip;
    ip.set_coeffs(xm);
    // FIXME, we're not using EM instead flds_em
    real_t E[3] = { ip.ex(flds_em), ip.ey(flds_em), ip.ez(flds_em) };
    real_t H[3] = { ip.hx(flds_em), ip.hy(flds_em), ip.hz(flds_em) };

    // x^(n+1/2), p^{n+1/2} -> x^(n+1/2), p^{n}
    int kind = prt->kind();
    real_t dq = dq_kind[kind];
    advance.push_p(&prt->pxi, E, H, -.5f * dq);
  }

  // ----------------------------------------------------------------------
  // push_mprts

  static void push_mprts(Mparticles& mprts, Mfields& mflds)
  {
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      auto& prts = mprts[p];

      flds.zero(JXI, JXI + 3);

      typename C::curr_cache_t curr_cache(flds);

      unsigned int n_prts = prts.size();
      for (int n = 0; n < n_prts; n++) {
	push(prts, n, flds, curr_cache);
      }
    }
  }

  // ----------------------------------------------------------------------
  // stagger_mprts
  
  static void stagger_mprts(Mparticles& mprts, Mfields& mflds)
  {
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto flds = mflds[p];
      auto& prts = mprts[p];

      flds.zero(JXI, JXI + 3);

      unsigned int n_prts = prts.size();
      for (int n = 0; n < n_prts; n++) {
	stagger(prts, n, flds);
      }
    }
  }
};

// ======================================================================

template<typename C>
void push_p_ops<C>::push_mprts(typename C::Mparticles& mprts,
			       typename C::Mfields& mflds)
{
  PushParticles1vb<C>::push_mprts(mprts, mflds);
}

template<typename C>
void push_p_ops<C>::push_mprts(struct psc_mparticles *mprts,
			       struct psc_mfields *mflds_base)
{
  auto mf = mflds_base->get_as<mfields_t>(EX, EX + 6);
  auto mp = mparticles_t(mprts);
  PushParticles1vb<C>::push_mprts(*mp.sub(), *mf.sub());
  mf.put_as(mflds_base, JXI, JXI+3);
}

template<typename C>
void push_p_ops<C>::stagger_mprts(struct psc_mparticles *mprts,
				  struct psc_mfields *mflds_base)
{
  auto mf = mflds_base->get_as<mfields_t>(EX, EX + 6);
  auto mp = mparticles_t(mprts);
  PushParticles1vb<C>::stagger_mprts(*mp.sub(), *mf.sub());
  mf.put_as(mflds_base, JXI, JXI+3);
}

template struct push_p_ops<push_p_conf>;
