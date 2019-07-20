
#ifndef KG_SARRAY_H
#define KG_SARRAY_H

#include <kg/SArrayContainer.h>

namespace kg
{

// ======================================================================
// StorageUniquePtr
//
// Mostly, std::vector would be an equivalent choice, though this
// Storagle class is move-only, preventing accidental copies

template <typename T>
class StorageUniquePtr
{
public:
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;

  StorageUniquePtr(size_t size) : data_{new T[size]} {}

  const_reference operator[](int offset) const { return data_[offset]; }
  reference operator[](int offset) { return data_[offset]; }

  // FIXME access to underlying storage might better be avoided?
  // use of this makes assumption that storage is contiguous
  const_pointer data() const { return data_.get(); }
  pointer data() { return data_.get(); }

private:
  std::unique_ptr<value_type[]> data_;
};

// ======================================================================
// SArray

template <typename T, typename L = LayoutSOA>
struct SArray;

template <typename T, typename L>
struct SArrayContainerInnerTypes<SArray<T, L>>
{
  using Layout = L;
  using Storage = StorageUniquePtr<T>;
  //using Storage = std::vector<T>;
};

template <typename T, typename L>
struct SArray : SArrayContainer<SArray<T, L>>
{
  using Base = SArrayContainer<SArray<T, L>>;
  using Storage = typename Base::Storage;
  using real_t = typename Base::value_type;

  SArray(Int3 ib, Int3 im, int n_comps)
    : Base{ib, im, n_comps},
      storage_(Base::size())
  {}

private:
  Storage storage_;

  Storage& storageImpl() { return storage_; }
  const Storage& storageImpl() const { return storage_; }

  friend class SArrayContainer<SArray<T, L>>;
};

} // namespace kg

#endif
