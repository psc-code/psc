
#pragma once

#include "cuda_compat.h"

// ======================================================================
// ParticleSimple

template<typename _Real>
struct ParticleSimple
{
  using real_t = _Real;
  using Real3 = Vec3<real_t>;

  ParticleSimple() = default;

  ParticleSimple(const ParticleSimple&) = default;

  __host__ __device__
  ParticleSimple(Real3 x, Real3 u, real_t qni_wni, int kind)
    : xi_{x},
      pxi_{u},
      kind_{kind},
      qni_wni_{qni_wni}
  {}

  __host__ __device__
  bool operator==(const ParticleSimple& other) const
  {
    return (xi_ == other.xi_ && qni_wni_ == other.qni_wni_ &&
	    pxi_ == other.pxi_ && kind_ == other.kind_);
  }

  __host__ __device__
  bool operator!=(const ParticleSimple& other) const { return !(*this == other); }

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

