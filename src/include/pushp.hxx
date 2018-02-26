
#ifndef PUSHP_HXX
#define PUSHP_HXX

#include "cuda_compat.h"

template<typename real_t, typename dim>
struct AdvanceParticle
{
  __host__ __device__ AdvanceParticle(real_t dt)
    : dt_(dt)
  {}
  
  // ----------------------------------------------------------------------
  // push_x
  
  __host__ __device__ inline void push_x(real_t x[3], const real_t v[3], real_t dt_fac=real_t{1.})
  {
    if (!dim::InvarX::value) { x[0] += dt_fac * dt_ * v[0]; }
    if (!dim::InvarY::value) { x[1] += dt_fac * dt_ * v[1]; }
    if (!dim::InvarZ::value) { x[2] += dt_fac * dt_ * v[2]; }
  }

  // ----------------------------------------------------------------------
  // push_p
  //
  // advance momentum based on Lorentz force from EM fields
  
  __host__ __device__ inline void push_p(real_t p[3], const real_t E[3], const real_t H[3], real_t dq)
  {
    real_t pxm = p[0] + dq * E[0];
    real_t pym = p[1] + dq * E[1];
    real_t pzm = p[2] + dq * E[2];
    
    real_t root = dq * rsqrt(real_t(1.) + sqr(pxm) + sqr(pym) + sqr(pzm));
    real_t taux = H[0] * root, tauy = H[1] * root, tauz = H[2] * root;
    
    real_t tau = real_t(1.) / (real_t(1.) + sqr(taux) + sqr(tauy) + sqr(tauz));
    real_t pxp = ( (real_t(1.) + sqr(taux) - sqr(tauy) - sqr(tauz)) * pxm
	      +(real_t(2.)*taux*tauy + real_t(2.)*tauz)*pym
	      +(real_t(2.)*taux*tauz - real_t(2.)*tauy)*pzm)*tau;
    real_t pyp = ( (real_t(2.)*taux*tauy - real_t(2.)*tauz)*pxm
	      +(real_t(1.) - sqr(taux) + sqr(tauy) - sqr(tauz)) * pym
	      +(real_t(2.)*tauy*tauz + real_t(2.)*taux)*pzm)*tau;
    real_t pzp = ( (real_t(2.)*taux*tauz + real_t(2.)*tauy)*pxm
	      +(real_t(2.)*tauy*tauz - real_t(2.)*taux)*pym
	      +(real_t(1.) - sqr(taux) - sqr(tauy) + sqr(tauz))*pzm)*tau;
    
    p[0] = pxp + dq * E[0];
    p[1] = pyp + dq * E[1];
    p[2] = pzp + dq * E[2];
  }

// ----------------------------------------------------------------------
// calc_v

  __host__ __device__ inline void calc_v(real_t v[3], const real_t p[3])
  {
    real_t root = rsqrt(1.f + sqr(p[0]) + sqr(p[1]) + sqr(p[2]));
    for (int d = 0; d < 3; d++) {
      v[d] = p[d] * root;
    }
  }
  
private:
  real_t dt_;
};

#endif
