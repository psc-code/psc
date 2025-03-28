
#pragma once

#include <kg/Vec3.h>

#include <cstdint>

namespace psc
{
namespace particle
{

using Id = uint64_t;
using Tag = int;

// ======================================================================
// Inject
//
// standard type used to pass new particles to a particle injector class
// (which will take care of converting it the specific actual storage type)

struct Inject
{
  using Real = double;
  using Real3 = Vec3<Real>;

  Inject(const Real3& x, const Real3& u, Real w, int kind, Tag tag = {})
    : x{x}, u{u}, w{w}, kind{kind}, tag{tag}
  {}

  Real3 x;
  Real3 u;
  Real w;
  int kind;
  Tag tag;
};

} // namespace particle
} // namespace psc

inline std::ostream& operator<<(std::ostream& os,
                                const psc::particle::Inject& inj)
{
  os << "Inject{x=" << inj.x << ", u=" << inj.u << ", w=" << inj.w
     << ", kind=" << inj.kind << ", tag=" << inj.tag << "}";
  return os;
}