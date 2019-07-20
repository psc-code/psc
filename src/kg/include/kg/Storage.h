
#ifndef KG_STORAGE_H
#define KG_STORAGE_H

namespace kg
{

// ======================================================================
// StorageNoOwnership

template <typename T>
class StorageNoOwnership
{
public:
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;

  KG_INLINE StorageNoOwnership(pointer data) : data_{data} {}

  KG_INLINE const_reference operator[](int offset) const { return data_[offset]; }
  KG_INLINE reference operator[](int offset) { return data_[offset]; }

  // FIXME access to underlying storage might better be avoided?
  // use of this makes assumption that storage is contiguous
  const_pointer data() const { return data_; }
  pointer data() { return data_; }

private:
  pointer data_;
};

} // namespace kg

#endif
