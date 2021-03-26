
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

  KG_INLINE SArrayContainer(const Box3& box) : box_{box} {}

  KG_INLINE const Int3& ib() const { return box_.ib(); }

  KG_INLINE const_reference operator()(int m, int i, int j, int k) const
  {
    return storage()(i - ib()[0], j - ib()[1], k - ib()[2], m);
  }

  KG_INLINE reference operator()(int m, int i, int j, int k)
  {
    return storage()(i - ib()[0], j - ib()[1], k - ib()[2], m);
  }

public:
  KG_INLINE Storage& storage() { return derived().storageImpl(); }
  KG_INLINE const Storage& storage() const { return derived().storageImpl(); }

protected:
  KG_INLINE Derived& derived() { return *static_cast<Derived*>(this); }
  KG_INLINE const Derived& derived() const
  {
    return *static_cast<const Derived*>(this);
  }

private:
  Box3 box_; //> lower bounds and length per direction
};

} // namespace kg

#endif
