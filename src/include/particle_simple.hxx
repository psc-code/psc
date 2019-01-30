
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
    : x_{x},
      u_{u},
      kind{kind},
      qni_wni{qni_wni}
  {}

  __host__ __device__
  bool operator==(const ParticleSimple& other) const
  {
    return (x_ == other.x_ && qni_wni == other.qni_wni &&
	    u_ == other.u_ && kind == other.kind);
  }

  __host__ __device__
  bool operator!=(const ParticleSimple& other) const { return !(*this == other); }

  __host__ __device__ Real3  x() const { return x_; }
  __host__ __device__ Real3& x()       { return x_; }
  __host__ __device__ Real3  u() const { return u_; }
  __host__ __device__ Real3& u()       { return u_; }
  
private:
  Real3 x_;
  Real3 u_;
public:
  int kind;
  real_t qni_wni;
};

