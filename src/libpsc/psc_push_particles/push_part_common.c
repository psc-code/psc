
#include <mrc_profile.h>

#include <particle_iter.h>

template<typename C>
struct PushParticles__
{
  using Mparticles = typename C::Mparticles;
  using Mfields = typename C::Mfields;
  using fields_t = typename Mfields::fields_t;
  using AdvanceParticle_t = typename C::AdvanceParticle_t;
  using real_t = typename Mparticles::real_t;
  using Real3 = Vec3<real_t>;

  static void push_mprts(Mparticles& mprts, Mfields& mflds)
  {
    static int pr;
    if (!pr) {
      pr = prof_register(__func__, 1., 0, 0);
    }

    prof_start(pr);
    for (int p = 0; p < mprts.n_patches(); p++) {
      mflds[p].zero(JXI, JXI + 3);
      do_push_part(mflds[p], mprts[p]);
    }
    prof_stop(pr);
  }

private:
  static void do_push_part(fields_t flds, typename Mparticles::patch_t& prts)
  {
    using dim = typename C::dim;
    using InterpolateEM_t = typename C::InterpolateEM_t;
    using CurrentE_t = typename C::CurrentE_t;
    using particle_t = typename Mparticles::particle_t;

    real_t dqs = .5f * prts.grid().eta * prts.grid().dt;
    Real3 dxi = Real3{ 1., 1., 1. } / Real3(prts.grid().dx);
  
    AdvanceParticle_t advance(prts.grid().dt);
    InterpolateEM_t ip;
    CurrentE_t c(prts.grid());

    Fields3d<fields_t> EM(flds); // FIXME, EM and J are identical here
    Fields3d<fields_t> J(flds);

    PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
      particle_t *part = &*prt_iter;
      real_t *x = &part->xi;
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
      real_t dq = dqs * prts.prt_qni(*part) / prts.prt_mni(*part);
      advance.push_p(&part->pxi, E, H, dq);

      // x^(n+0.5), p^(n+1.0) -> x^(n+1.0), p^(n+1.0)
      advance.calc_v(vv, &part->pxi);
      advance.push_x(x, vv);

      // CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt
      c.charge_after(x);

      // CURRENT DENSITY AT (n+1.0)*dt
      c.prep(prts.prt_qni_wni(*part), vv);
      c.calc(J);
    }
  }

};

// ----------------------------------------------------------------------
// PscPushParticles_

template<typename PushParticles_t>
struct PscPushParticles_
{
  using Mparticles = typename PushParticles_t::Mparticles;
  using Mfields = typename PushParticles_t::Mfields;
  using mparticles_t = PscMparticles<Mparticles>;
  using mfields_t = PscMfields<Mfields>;

  static void push_mprts(struct psc_mparticles *mprts,
			 struct psc_mfields *mflds_base)
  {
    mfields_t mf = mflds_base->get_as<mfields_t>(EX, EX + 6);
    mparticles_t mp = mparticles_t(mprts);

    PushParticles_t::push_mprts(*mp.sub(), *mf.sub());

    mf.put_as(mflds_base, JXI, JXI+3);
  }
};

#ifdef CONFIG
template struct PscPushParticles_<PushParticles__<CONFIG>>;
#endif
