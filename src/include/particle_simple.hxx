
#pragma once

// ======================================================================
// ParticleSimple

template<class R>
struct ParticleSimple
{
  using real_t = R;
  using Real3 = Vec3<real_t>;

  ParticleSimple()
  {}

  ParticleSimple(Real3 x, Real3 u, real_t qni_wni, int kind)
    : x_{x}, u_{u}, qni_wni_{qni_wni}, kind_{kind}
  {}

  Real3  x() const { return x_; }
  Real3& x(  )     { return x_; }
  Real3  u() const { return u_; }
  Real3& u(  )     { return u_; }
  int kind() const { return kind_; }
  real_t qni_wni() const { return qni_wni_; }
  
  Real3 x_;
  real_t qni_wni_;
  Real3 u_;
  int kind_;
};

