
#pragma once

#include "cuda_compat.h"

// ======================================================================
// DParticleSimple

template<typename _Real>
struct DParticleSimple
{
  using real_t = _Real;
  using Real3 = Vec3<real_t>;

  __host__ __device__
  DParticleSimple() = default;

  __host__ __device__
  DParticleSimple(const DParticleSimple&) = default;

  __host__ __device__
  DParticleSimple(Real3 x, Real3 u, real_t qni_wni, int kind)
    : xi_{x},
      pxi_{u},
      kind_{kind},
      qni_wni_{qni_wni}
  {}

  __host__ __device__
  bool operator==(const DParticleSimple& other) const
  {
    return (xi_ == other.xi_ && qni_wni_ == other.qni_wni_ &&
	    pxi_ == other.pxi_ && kind_ == other.kind_);
  }

  __host__ __device__
  bool operator!=(const DParticleSimple& other) const { return !(*this == other); }

  __host__ __device__ Real3  x() const { return xi_; }
  __host__ __device__ Real3& x()       { return xi_; }
  __host__ __device__ Real3  u() const { return pxi_; }
  __host__ __device__ Real3& u()       { return pxi_; }
  __host__ __device__ int kind() const { return kind_; }
  __host__ __device__ real_t qni_wni() const { return qni_wni_; }
  
private:
  Real3 xi_;
  Real3 pxi_;
  int kind_;
  real_t qni_wni_;
};

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

