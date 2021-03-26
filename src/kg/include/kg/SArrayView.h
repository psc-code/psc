
#ifndef KG_SARRAY_VIEW_H
#define KG_SARRAY_VIEW_H

#include <kg/SArrayContainer.h>

#include <gtensor/gtensor.h>

namespace kg
{

// ======================================================================
// SArrayView

template <typename T, typename L = kg::LayoutSOA>
struct SArrayView;

template <typename T, typename L>
struct SArrayContainerInnerTypes<SArrayView<T, L>>
{
  using Storage = gt::gtensor_span<T, 4>;
};

namespace detail
{
template <typename L>
gt::shape_type<4> strides(const Box3& box, int n_comps);

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

template <typename T, typename L>
struct SArrayView : kg::SArrayContainer<SArrayView<T, L>>
{
  using Base = kg::SArrayContainer<SArrayView<T, L>>;
  using Storage = typename Base::Storage;
  using real_t = typename Base::value_type;

  KG_INLINE SArrayView(const Box3& box, int n_comps, real_t* data)
    : Base(box),
      storage_(data, gt::shape(box.im(0), box.im(1), box.im(2), n_comps),
               detail::strides<L>(box, n_comps))
  {}

private:
  Storage storage_;

  KG_INLINE Storage& storageImpl() { return storage_; }
  KG_INLINE const Storage& storageImpl() const { return storage_; }

  friend class kg::SArrayContainer<SArrayView<T, L>>;
};

} // namespace kg

#endif
