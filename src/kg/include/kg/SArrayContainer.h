
#ifndef KG_SARRAY_CONTAINER_H
#define KG_SARRAY_CONTAINER_H

#include <kg/Array3d.h>
#include <kg/Macros.h>
#include <kg/Vec3.h>

#include <gtensor/gtensor.h>

#include <cstring>

using namespace gt::placeholders;

// FIXME, do noexcept?
// FIXME, use size_t instead of int, at least for 1d offsets?

namespace kg
{

// ======================================================================
// SArrayContainer

template <typename C>
struct SArrayContainerInnerTypes;

template <typename D>
class SArrayContainer
{
public:
  using Derived = D;

  using InnerTypes = SArrayContainerInnerTypes<D>;
  using Storage = typename InnerTypes::Storage;
  using value_type = typename Storage::value_type;
  using reference = typename Storage::reference;
  using const_reference = typename Storage::const_reference;
  using pointer = typename Storage::pointer;
  using const_pointer = typename Storage::const_pointer;
};

} // namespace kg

#endif
