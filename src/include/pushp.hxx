
#ifndef PUSHP_HXX
#define PUSHP_HXX

#include "cuda_compat.h"

template<typename real_t, typename dim>
struct AdvanceParticle
{
  // ----------------------------------------------------------------------
  // push_x
  
  __host__ __device__ inline void push_x(real_t x[3], const real_t v[3], real_t dt)
  {
    if (!dim::InvarX::value) { x[0] += dt * v[0]; }
    if (!dim::InvarY::value) { x[1] += dt * v[1]; }
    if (!dim::InvarZ::value) { x[2] += dt * v[2]; }
  }
};

// ----------------------------------------------------------------------
// push_p
//
// advance momentum based on Lorentz force from EM fields

template<typename R>
__host__ __device__ inline void push_p(R p[3], const R E[3], const R H[3], R dq)
{
  R pxm = p[0] + dq * E[0];
  R pym = p[1] + dq * E[1];
  R pzm = p[2] + dq * E[2];
  
  R root = dq * rsqrt(R(1.) + sqr(pxm) + sqr(pym) + sqr(pzm));
  R taux = H[0] * root, tauy = H[1] * root, tauz = H[2] * root;
  
  R tau = R(1.) / (R(1.) + sqr(taux) + sqr(tauy) + sqr(tauz));
  R pxp = ( (R(1.) + sqr(taux) - sqr(tauy) - sqr(tauz)) * pxm
           +(R(2.)*taux*tauy + R(2.)*tauz)*pym
	   +(R(2.)*taux*tauz - R(2.)*tauy)*pzm)*tau;
  R pyp = ( (R(2.)*taux*tauy - R(2.)*tauz)*pxm
	   +(R(1.) - sqr(taux) + sqr(tauy) - sqr(tauz)) * pym
	   +(R(2.)*tauy*tauz + R(2.)*taux)*pzm)*tau;
  R pzp = ( (R(2.)*taux*tauz + R(2.)*tauy)*pxm
	   +(R(2.)*tauy*tauz - R(2.)*taux)*pym
	   +(R(1.) - sqr(taux) - sqr(tauy) + sqr(tauz))*pzm)*tau;
  
  p[0] = pxp + dq * E[0];
  p[1] = pyp + dq * E[1];
  p[2] = pzp + dq * E[2];
}

// ----------------------------------------------------------------------
// calc_v

template<typename R>
__host__ __device__ inline void calc_v(R v[3], const R p[3])
{
  R root = rsqrt(1.f + sqr(p[0]) + sqr(p[1]) + sqr(p[2]));
  for (int d = 0; d < 3; d++) {
    v[d] = p[d] * root;
  }
}

#endif
