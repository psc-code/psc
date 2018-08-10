
#include <mrc_profile.h>

#include "../libpsc/psc_checks/checks_impl.hxx"

template<typename C>
struct PushParticles__
{
  using Mparticles = typename C::Mparticles;
  using MfieldsState = typename C::MfieldsState;
  using fields_t = typename MfieldsState::fields_t;
  using AdvanceParticle_t = typename C::AdvanceParticle_t;
  using real_t = typename Mparticles::real_t;
  using Real3 = Vec3<real_t>;
  using checks_order = checks_order_2nd;

  static void push_mprts(Mparticles& mprts, MfieldsState& mflds)
  {
    static int pr;
    if (!pr) {
      pr = prof_register(__func__, 1., 0, 0);
    }

    prof_start(pr);
    for (int p = 0; p < mprts.n_patches(); p++) {
      mflds[p].zero(JXI, JXI + 3);
      push_mprts_patch(mflds[p], mprts[p]);
    }
    prof_stop(pr);
  }

private:
  static void push_mprts_patch(fields_t flds, typename Mparticles::patch_t& prts)
  {
    using InterpolateEM_t = typename C::InterpolateEM_t;
    using CurrentE_t = typename C::CurrentE_t;
    using particle_t = typename Mparticles::particle_t;

    real_t dqs = .5f * prts.grid().norm.eta * prts.grid().dt;
    Real3 dxi = Real3{ 1., 1., 1. } / Real3(prts.grid().domain.dx);
  
    AdvanceParticle_t advance(prts.grid().dt);
    InterpolateEM_t ip;
    CurrentE_t c(prts.grid());

    Fields3d<fields_t> EM(flds); // FIXME, EM and J are identical here
    Fields3d<fields_t> J(flds);

    for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
      particle_t& prt = *prt_iter;
      real_t *x = prt.x;
      real_t vv[3];

      // CHARGE DENSITY FORM FACTOR AT (n+.5)*dt
      // FIELD INTERPOLATION

      real_t xm[3];
      for (int d = 0; d < 3; d++) {
	xm[d] = x[d] * dxi[d];
      }
      ip.set_coeffs(xm);

      c.charge_before(ip);

      real_t E[3] = { ip.ex(EM), ip.ey(EM), ip.ez(EM) };
      real_t H[3] = { ip.hx(EM), ip.hy(EM), ip.hz(EM) };

      // x^(n+0.5), p^n -> x^(n+0.5), p^(n+1.0)
      real_t dq = dqs * prts.prt_qni(prt) / prts.prt_mni(prt);
      advance.push_p(prt.p, E, H, dq);

      // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0)
      advance.calc_v(vv, prt.p);
      advance.push_x(x, vv);

      // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt
      c.charge_after(x);

      // CURRENT DENSITY AT (n+1.0)*dt
      c.prep(prts.prt_qni_wni(prt), vv);
      c.calc(J);
    }
  }

};

