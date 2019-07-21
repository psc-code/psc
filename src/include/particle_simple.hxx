
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

  KG_INLINE ParticleSimple(Real3 x, Real3 u, real_t qni_wni, int kind)
    : x{x},
      u{u},
      kind{kind},
      qni_wni{qni_wni}
  {}

  KG_INLINE bool operator==(const ParticleSimple& other) const
  {
    return (x == other.x && qni_wni == other.qni_wni &&
	    u == other.u && kind == other.kind);
  }

  KG_INLINE bool operator!=(const ParticleSimple& other) const { return !(*this == other); }

public:
  Real3 x;
  Real3 u;
  int kind;
  real_t qni_wni;
};

