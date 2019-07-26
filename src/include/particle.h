
#pragma once

#include <kg/Vec3.h>

#include <cstdint>

namespace psc
{
namespace particle
{

using Id = uint64_t;

// ======================================================================
// Inject
//
// standard type used to pass new particles to a particle injector class
// (which will take care of converting it the specific actual storage type)

struct Inject
{
  using Real = double;
  using Real3 = Vec3<Real>;

  Real3 x;
  Real3 u;
  Real w;
  int kind;
  Id id;
};

} // namespace particle
} // namespace psc
