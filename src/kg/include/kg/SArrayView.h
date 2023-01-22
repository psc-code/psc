
#ifndef KG_SARRAY_VIEW_H
#define KG_SARRAY_VIEW_H

#include <kg/Array3d.h>
#include <kg/Macros.h>
#include <kg/Vec3.h>

#include <psc/gtensor.h>

using namespace gt::placeholders;

namespace kg
{

// ======================================================================
// SArrayView

namespace detail
{
template <typename L>
GT_INLINE gt::shape_type<4> strides(const Box3& box, int n_comps);

template <>
GT_INLINE gt::shape_type<4> strides<LayoutSOA>(const Box3& box, int n_comps)
{
  return gt::shape(1, box.im(0), box.im(0) * box.im(1),
                   box.im(0) * box.im(1) * box.im(2));
}

template <>
GT_INLINE gt::shape_type<4> strides<LayoutAOS>(const Box3& box, int n_comps)
{
  return gt::shape(n_comps, n_comps * box.im(0),
                   n_comps * box.im(0) * box.im(1), 1);
}

} // namespace detail

template <typename T, typename L = kg::LayoutSOA>
struct SArrayView
{
  using Storage = gt::gtensor_span<T, 4>;
  using value_type = typename Storage::value_type;
  using reference = typename Storage::reference;
  using const_reference = typename Storage::const_reference;
  using pointer = typename Storage::pointer;
  using const_pointer = typename Storage::const_pointer;

  KG_INLINE SArrayView(const Box3& box, int n_comps, pointer data)
    : ib_(box.ib()),
      storage_(data, gt::shape(box.im(0), box.im(1), box.im(2), n_comps),
               detail::strides<L>(box, n_comps))
  {}

  KG_INLINE const Int3& ib() const { return ib_; }
  KG_INLINE Storage& storage() { return storage_; }
  KG_INLINE const Storage& storage() const { return storage_; }

private:
  Storage storage_;
  Int3 ib_; //> lower bounds per direction
};

} // namespace kg

#endif
