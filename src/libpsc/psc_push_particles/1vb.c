
#include <mrc_profile.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "1vb/psc_push_particles_1vb.h"
#include "../libpsc/psc_checks/checks_impl.hxx"

#define MAX_NR_KINDS (10)

// ======================================================================
// PushParticles1vb

template<typename C>
struct PushParticles1vb
{
  using Mparticles = typename C::Mparticles;
  using MfieldsState = typename C::MfieldsState;
  using real_t = typename Mparticles::real_t;
  using Real3 = Vec3<real_t>;
  using dim = typename C::dim;
  using InterpolateEM_t = typename C::InterpolateEM_t;
  using AdvanceParticle_t = typename C::AdvanceParticle_t;
  using Current = typename C::Current_t;
  using particle_t = typename Mparticles::particle_t;
  using curr_cache_t = typename Current::curr_cache_t;
  using checks_order = checks_order_1st;
  
  // ----------------------------------------------------------------------
  // push_mprts

  static void push_mprts(Mparticles& mprts, MfieldsState& mflds)
  {
    for (int p = 0; p < mprts.n_patches(); p++) {
      mflds[p].zero(JXI, JXI + 3);
      push_mprts_patch(mflds[p], mprts[p]);
    }
  }

  // ----------------------------------------------------------------------
  // stagger_mprts
  
  static void stagger_mprts(Mparticles& mprts, MfieldsState& mflds)
  {
    for (int p = 0; p < mprts.n_patches(); p++) {
      stagger_mprts_patch(mflds[p], mprts[p]);
    }
  }

private:

  // ----------------------------------------------------------------------
  // push_mprts_patch
  
  static void push_mprts_patch(typename MfieldsState::fields_t flds, typename Mparticles::patch_t& prts)
  {
    typename InterpolateEM_t::fields_t EM(flds);
    InterpolateEM_t ip;
    AdvanceParticle_t advance(prts.grid().dt);
    curr_cache_t curr_cache(flds);
    Current current(prts.grid());

    PI<real_t> pi(prts.grid());
    Real3 dxi = Real3{ 1., 1., 1. } / Real3(prts.grid().domain.dx);
    real_t dq_kind[MAX_NR_KINDS];
    auto& kinds = prts.grid().kinds;
    assert(kinds.size() <= MAX_NR_KINDS);
    for (int k = 0; k < kinds.size(); k++) {
      dq_kind[k] = .5f * prts.grid().norm.eta * prts.grid().dt * kinds[k].q / kinds[k].m;
    }

    unsigned int n_prts = prts.size();
    for (int n = 0; n < n_prts; n++) {
      particle_t& prt = prts[n];
  
      // field interpolation
      real_t *xi = prt.x;

      real_t xm[3];
      for (int d = 0; d < 3; d++) {
	xm[d] = xi[d] * dxi[d];
      }

      // FIELD INTERPOLATION

      ip.set_coeffs(xm);
      // FIXME, we're not using EM instead flds_em
      real_t E[3] = { ip.ex(EM), ip.ey(EM), ip.ez(EM) };
      real_t H[3] = { ip.hx(EM), ip.hy(EM), ip.hz(EM) };

      // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
      real_t dq = dq_kind[prt.kind];
      advance.push_p(prt.p, E, H, dq);

      real_t vxi[3];
      advance.calc_v(vxi, prt.p);

      int lf[3];
      real_t of[3], xp[3];
#if CALC_J == CALC_J_1VB_2D
      // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0)
      advance.push_x(prt.x, vxi, .5f);
  
      // OUT OF PLANE CURRENT DENSITY AT (n+1.0)*dt
      pi.find_idx_off_1st_rel(prt.x, lf, of, real_t(0.));
      current.calc_j_oop(curr_cache, prts.prt_qni_wni(prt), vxi, lf, of);
  
      // x^(n+1), p^(n+1) -> x^(n+1.5), p^(n+1)
      advance.push_x(prt.x, vxi, .5f);
#else
      // x^(n+0.5), p^(n+1.0) -> x^(n+1.5), p^(n+1.0)
      advance.push_x(prt.x, vxi);
#endif
  
      pi.find_idx_off_pos_1st_rel(prt.x, lf, of, xp, real_t(0.));

      // CURRENT DENSITY BETWEEN (n+.5)*dt and (n+1.5)*dt
      int lg[3];
      if (!dim::InvarX::value) { lg[0] = ip.cx.g.l; }
      if (!dim::InvarY::value) { lg[1] = ip.cy.g.l; }
      if (!dim::InvarZ::value) { lg[2] = ip.cz.g.l; }
      current.calc_j(curr_cache, xm, xp, lf, lg, prts.prt_qni_wni(prt), vxi);

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
  }

  // ----------------------------------------------------------------------
  // stagger_mprts_patch
  
  static void stagger_mprts_patch(typename MfieldsState::fields_t flds, typename Mparticles::patch_t& prts)
  {
    typename InterpolateEM_t::fields_t EM(flds);
    InterpolateEM_t ip;
        
    AdvanceParticle_t advance(prts.grid().dt);

    Real3 dxi = Real3{ 1., 1., 1. } / Real3(prts.grid().domain.dx);
    real_t dq_kind[MAX_NR_KINDS];
    auto& kinds = prts.grid().kinds;
    assert(kinds.size() <= MAX_NR_KINDS);
    for (int k = 0; k < kinds.size(); k++) {
      dq_kind[k] = .5f * prts.grid().eta * prts.grid().dt * kinds[k].q / kinds[k].m;
    }
      
    unsigned int n_prts = prts.size();
    for (int n = 0; n < n_prts; n++) {
      particle_t *prt = &prts[n];
      
      // field interpolation
      real_t *xi = &prt->xi;
      
      real_t xm[3];
      for (int d = 0; d < 3; d++) {
	xm[d] = xi[d] * dxi[d];
      }
      
      // FIELD INTERPOLATION

      ip.set_coeffs(xm);
      // FIXME, we're not using EM instead flds_em
      real_t E[3] = { ip.ex(EM), ip.ey(EM), ip.ez(EM) };
      real_t H[3] = { ip.hx(EM), ip.hy(EM), ip.hz(EM) };
      
      // x^(n+1/2), p^{n+1/2} -> x^(n+1/2), p^{n}
      int kind = prt->kind();
      real_t dq = dq_kind[kind];
      advance.push_p(&prt->pxi, E, H, -.5f * dq);
    }
  }
};
