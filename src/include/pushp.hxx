
#ifndef PUSHP_HXX
#define PUSHP_HXX

#include "cuda_compat.h"

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

#endif
