
#pragma once

#include "cuda_compat.h"
#include "particles.hxx"

// ======================================================================
// ParticleSimple

template <typename _Real>
struct ParticleSimple
{
  using real_t = _Real;
  using Real3 = Vec3<real_t>;

  ParticleSimple() = default;

  KG_INLINE ParticleSimple(Real3 x, Real3 u, real_t qni_wni, int kind,
                           psc::particle::Id id, psc::particle::Tag tag)
    : x{x}, u{u}, kind{kind}, qni_wni{qni_wni}
  {}

  KG_INLINE bool operator==(const ParticleSimple& other) const
  {
    return (x == other.x && qni_wni == other.qni_wni && u == other.u &&
            kind == other.kind);
  }

  KG_INLINE bool operator!=(const ParticleSimple& other) const
  {
    return !(*this == other);
  }

  KG_INLINE psc::particle::Id id() const { return 0; }
  KG_INLINE psc::particle::Tag tag() const { return 0; }

public:
  Real3 x;
  Real3 u;
  int kind;
  real_t qni_wni;
};

template <typename R>
class ForComponents<ParticleSimple<R>>
{
public:
  using Particle = ParticleSimple<R>;

  template <typename FUNC>
  static void run(FUNC&& func)
  {
    func("x", [](Particle& prt) { return &prt.x[0]; });
    func("y", [](Particle& prt) { return &prt.x[1]; });
    func("z", [](Particle& prt) { return &prt.x[2]; });
    func("ux", [](Particle& prt) { return &prt.u[0]; });
    func("uy", [](Particle& prt) { return &prt.u[1]; });
    func("uz", [](Particle& prt) { return &prt.u[2]; });
    func("kind", [](Particle& prt) { return &prt.kind; });
    func("qni_wni", [](Particle& prt) { return &prt.qni_wni; });
  }
};
