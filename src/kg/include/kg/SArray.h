
#ifndef KG_SARRAY_H
#define KG_SARRAY_H

#include <kg/SArrayContainer.h>

#include <gtensor/gtensor.h>

namespace kg
{

// ======================================================================
// SArray

template <typename T, typename L = LayoutSOA>
struct SArray;

template <typename T, typename L>
struct SArrayContainerInnerTypes<SArray<T, L>>
{
  using Storage = gt::gtensor<T, 4>;
};

template <typename T, typename L>
struct SArray : SArrayContainer<SArray<T, L>>
{
  using Base = SArrayContainer<SArray<T, L>>;
  using Storage = typename Base::Storage;
  using real_t = typename Base::value_type;

  SArray(const Box3& box, int n_comps)
    : Base(box.ib()),
      storage_(gt::shape(box.im(0), box.im(1), box.im(2), n_comps))
  {
    static_assert(std::is_same<L, LayoutSOA>::value,
                  "FIXME need to support LayoutAOS");
  }

private:
  Storage storage_;

  KG_INLINE Storage& storageImpl() { return storage_; }
  KG_INLINE const Storage& storageImpl() const { return storage_; }

  friend class SArrayContainer<SArray<T, L>>;
};

} // namespace kg

#endif
