
#ifndef PUSHP_HXX
#define PUSHP_HXX

#include "cuda_compat.h"

// ----------------------------------------------------------------------
// push_p

template<typename R>
inline void push_p(R p[3], const R E[3], const R H[3], R dq)
{
  R pxm = p[0] + dq * E[0];
  R pym = p[1] + dq * E[1];
  R pzm = p[2] + dq * E[2];

  R root = dq / std::sqrt(1.f + pxm*pxm + pym*pym + pzm*pzm);
  R taux = H[0] * root;
  R tauy = H[1] * root;
  R tauz = H[2] * root;
  
  R tau = 1.f / (1.f + taux*taux + tauy*tauy + tauz*tauz);
  R pxp = ((1.f+taux*taux-tauy*tauy-tauz*tauz)*pxm + 
			 (2.f*taux*tauy+2.f*tauz)*pym + 
			 (2.f*taux*tauz-2.f*tauy)*pzm)*tau;
  R pyp = ((2.f*taux*tauy-2.f*tauz)*pxm +
			 (1.f-taux*taux+tauy*tauy-tauz*tauz)*pym +
			 (2.f*tauy*tauz+2.f*taux)*pzm)*tau;
  R pzp = ((2.f*taux*tauz+2.f*tauy)*pxm +
			 (2.f*tauy*tauz-2.f*taux)*pym +
			 (1.f-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau;
  
  p[0] = pxp + dq * E[0];
  p[1] = pyp + dq * E[1];
  p[2] = pzp + dq * E[2];
}

#endif
