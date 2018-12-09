
#ifndef PUSHP_HXX
#define PUSHP_HXX

#include "cuda_compat.h"

template<typename real_t, typename dim>
struct AdvanceParticle
{
  using Real3 = Vec3<real_t>;
  
  __host__ __device__ AdvanceParticle(real_t dt)
    : dt_(dt)
  {}
  
  // ----------------------------------------------------------------------
  // push_x
  
  __host__ __device__ inline void push_x(Real3& x, const Real3& v, real_t dt_fac=real_t{1.})
  {
    if (!dim::InvarX::value) { x[0] += dt_fac * dt_ * v[0]; }
    if (!dim::InvarY::value) { x[1] += dt_fac * dt_ * v[1]; }
    if (!dim::InvarZ::value) { x[2] += dt_fac * dt_ * v[2]; }
  }

  // ----------------------------------------------------------------------
  // push_p
  //
  // advance momentum based on Lorentz force from EM fields
  
  __host__ __device__ inline void push_p(Real3& p, const Real3& E, const Real3& H, real_t dq)
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

  __host__ __device__ inline Real3 calc_v(const Real3& u)
  {
    real_t root = rsqrt(1.f + sqr(u[0]) + sqr(u[1]) + sqr(u[2]));
    return u * Real3{root, root, root};
  }
  
private:
  real_t dt_;
};

#endif
