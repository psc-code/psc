
#include "psc_debug.h"

#include "interpolate.hxx"

// ======================================================================
// InterpolateEM

template<typename F, typename OPT_IP, typename OPT_DIM>
struct InterpolateEM
{
  using IP = InterpolateEM<F, OPT_IP, OPT_DIM>;
  using real_t = typename F::real_t;
  using ip_coeffs_t = ip_coeffs<real_t, OPT_IP>;
  using ip_coeff_t = typename ip_coeffs_t::ip_coeff_t;

  void set_coeffs(real_t xm[3])
  {
    IF_DIM_X( cx.set(xm[0]); );
    IF_DIM_Y( cy.set(xm[1]); );
    IF_DIM_Z( cz.set(xm[2]); );
  }

  using Helper = InterpolateEM_Helper<F, IP, OPT_IP, OPT_DIM>;
  real_t ex(F EM) { return Helper::ex(*this, EM); }
  real_t ey(F EM) { return Helper::ey(*this, EM); }
  real_t ez(F EM) { return Helper::ez(*this, EM); }
  real_t hx(F EM) { return Helper::hx(*this, EM); }
  real_t hy(F EM) { return Helper::hy(*this, EM); }
  real_t hz(F EM) { return Helper::hz(*this, EM); }

  ip_coeffs_t cx, cy, cz;
};

using IP = InterpolateEM<Fields3d<fields_t>, opt_ip, opt_dim>;

