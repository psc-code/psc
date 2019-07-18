
#ifndef KG_SARRAY_H
#define KG_SARRAY_H

#include <kg/SArrayContainer.h>

namespace kg
{

// ======================================================================
// StorageRaw

template <typename T>
class StorageRaw
{
public:
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;

  StorageRaw(pointer data) : data_{data} {}
  ~StorageRaw()
  {
    ::free(data_);
  }


  const_reference operator[](int offset) const { return data_[offset]; }
  reference operator[](int offset) { return data_[offset]; }

  // FIXME access to underlying storage might better be avoided?
  // use of this makes assumption that storage is contiguous
  const_pointer data() const { return data_; }
  pointer data() { return data_; }

private:
  pointer data_;
};

// ======================================================================
// SArray

template <typename T, typename L = LayoutSOA>
struct SArray;

template <typename T, typename L>
struct SArrayContainerInnerTypes<SArray<T, L>>
{
  using Layout = L;
  using Storage = StorageRaw<T>;
};

template <typename T, typename L>
struct SArray : SArrayContainer<SArray<T, L>>
{
  using Base = SArrayContainer<SArray<T, L>>;
  using Storage = typename Base::Storage;
  using real_t = typename Base::value_type;

  SArray(Int3 ib, Int3 im, int n_comps)
    : Base{ib, im, n_comps},
      storage_{(real_t*)calloc(Base::size(), sizeof(real_t))}
  {}

private:
  Storage storage_;

  Storage& storageImpl() { return storage_; }
  const Storage& storageImpl() const { return storage_; }

  friend class SArrayContainer<SArray<T, L>>;
};

} // namespace kg

#endif
