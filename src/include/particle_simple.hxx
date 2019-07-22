
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

template <typename Particle>
class DoComponents;

template <typename R>
class DoComponents<ParticleSimple<R>>
{
public:
  using Particle = ParticleSimple<R>;
  
  template <typename FUNC>
  static void run(FUNC&& put_component)
  {
    put_component("x", [](const Particle& prt) { return prt.x[0]; });
    put_component("y", [](const Particle& prt) { return prt.x[1]; });
    put_component("z", [](const Particle& prt) { return prt.x[2]; });
    put_component("ux", [](const Particle& prt) { return prt.u[0]; });
    put_component("uy", [](const Particle& prt) { return prt.u[1]; });
    put_component("uz", [](const Particle& prt) { return prt.u[2]; });
    put_component("kind", [](const Particle& prt) { return prt.kind; });
    put_component("qni_wni", [](const Particle& prt) { return prt.qni_wni; });
  }
};

