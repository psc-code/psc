
#ifndef KG_SARRAY_VIEW_H
#define KG_SARRAY_VIEW_H

#include <kg/SArrayContainer.h>
#include <kg/Storage.h>

namespace kg
{

// ======================================================================
// SArrayView

template <typename T, typename L = kg::LayoutSOA>
struct SArrayView;

template <typename T, typename L>
struct SArrayContainerInnerTypes<SArrayView<T, L>>
{
  using Layout = L;
  using Storage = StorageNoOwnership<T>;
};

template <typename T, typename L>
struct SArrayView : kg::SArrayContainer<SArrayView<T, L>>
{
  using Base = kg::SArrayContainer<SArrayView<T, L>>;
  using Storage = typename Base::Storage;
  using real_t = typename Base::value_type;

  KG_INLINE SArrayView(const Box3 &box, int n_comps, real_t* data)
    : Base{box, n_comps}, storage_{data}
  {}

private:
  Storage storage_;

  KG_INLINE Storage& storageImpl() { return storage_; }
  KG_INLINE const Storage& storageImpl() const { return storage_; }

  friend class kg::SArrayContainer<SArrayView<T, L>>;
};

} // namespace kg

#endif
