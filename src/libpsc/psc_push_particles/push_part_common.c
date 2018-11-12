
#include <mrc_profile.h>

#include "../libpsc/psc_checks/checks_impl.hxx"

template<typename C>
struct PushParticlesCommon
{
  using Mparticles = typename C::Mparticles;
  using MfieldsState = typename C::MfieldsState;
  using particle_t = typename Mparticles::particle_t;
  using real_t = typename Mparticles::real_t;
  using Real3 = Vec3<real_t>;
  using AdvanceParticle_t = typename C::AdvanceParticle_t;
  using InterpolateEM_t = typename C::InterpolateEM_t;
};

// ======================================================================
// PushParticles__

template<typename C>
struct PushParticles__ : PushParticlesCommon<C>
{
  using Base = PushParticlesCommon<C>;
  using typename Base::Mparticles;
  using typename Base::MfieldsState;
  using typename Base::particle_t;
  using typename Base::real_t;
  using typename Base::Real3;
  using typename Base::AdvanceParticle_t;
  using typename Base::InterpolateEM_t;
  
  using Dim = typename C::dim;
  using Current = typename C::Current_t;
  using checks_order = checks_order_2nd; // FIXME, sometimes 1st even with Esirkepov

  static void push_mprts(Mparticles& mprts, MfieldsState& mflds)
  {
    for (int p = 0; p < mprts.n_patches(); p++) {
      mflds[p].zero(JXI, JXI + 3);
      push_mprts_patch(mflds[p], mprts[p]);
    }
  }

private:

  // ----------------------------------------------------------------------
  // push_mprts_patch
  
  static void push_mprts_patch(typename MfieldsState::fields_t flds, typename Mparticles::patch_t& prts)
  {
    typename InterpolateEM_t::fields_t EM(flds);
    typename Current::fields_t J(flds);
    InterpolateEM_t ip;
    AdvanceParticle_t advance(prts.grid().dt);
    Current current(prts.grid());
  
    real_t dqs = .5f * prts.grid().norm.eta * prts.grid().dt;
    Real3 dxi = Real3{ 1., 1., 1. } / Real3(prts.grid().domain.dx);
  
    for (auto& prt: prts) {
      real_t *x = prt.x;
      real_t vv[3];

      // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt
      // FIELD INTERPOLATION

      real_t xm[3];
      for (int d = 0; d < 3; d++) {
	xm[d] = x[d] * dxi[d];
      }
      ip.set_coeffs(xm);

      current.charge_before(ip);

      real_t E[3] = { ip.ex(EM), ip.ey(EM), ip.ez(EM) };
      real_t H[3] = { ip.hx(EM), ip.hy(EM), ip.hz(EM) };

      // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
      real_t dq = dqs * prts.prt_qni(prt) / prts.prt_mni(prt);
      advance.push_p(prt.p, E, H, dq);

      // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0)
      advance.calc_v(vv, prt.p);
      advance.push_x(x, vv);

      // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt
      current.charge_after(x);

      // CURRENT DENSITY AT (n+1.0)*dt
      current.prep(prts.prt_qni_wni(prt), vv);
      current.calc(J);
    }
  }

};

