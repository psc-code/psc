
#ifndef KG_SARRAY_CONTAINER_H
#define KG_SARRAY_CONTAINER_H

#include <kg/Array3d.h>
#include <kg/Vec3.h>
#include <kg/Macros.h>

#include <cstring>

// FIXME, do noexcept?
// FIXME, use size_t instead of int, at least for 1d offsets?

namespace kg
{

namespace detail
{
template <typename Layout>
struct LayoutDataOffset;

template <>
struct LayoutDataOffset<LayoutSOA>
{
  KG_INLINE static int run(int n_comps, const Int3& im, int m, Int3 idx)
  {
    return (((((m)*im[2] + idx[2]) * im[1] + idx[1]) * im[0] + idx[0]));
  }
};

template <>
struct LayoutDataOffset<LayoutAOS>
{
  KG_INLINE static int run(int n_comps, const Int3& im, int m, Int3 idx)
  {
    return ((((idx[2]) * im[1] + idx[1]) * im[0] + idx[0]) * n_comps + m);
  }
};

} // namespace detail

template <typename Layout>
KG_INLINE static int layoutDataOffset(int n_comps, const Int3& im, int m, Int3 idx)
{
  return detail::LayoutDataOffset<Layout>::run(n_comps, im, m, idx);
}

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
  using Layout = typename InnerTypes::Layout;
  using Storage = typename InnerTypes::Storage;
  using value_type = typename Storage::value_type;
  using reference = typename Storage::reference;
  using const_reference = typename Storage::const_reference;
  using pointer = typename Storage::pointer;
  using const_pointer = typename Storage::const_pointer;

  KG_INLINE SArrayContainer(const Box3& box, int n_comps)
    : box_{box}, n_comps_{n_comps}
  {}

  KG_INLINE const Box3& box() const { return box_; }
  KG_INLINE const Int3& ib() const { return box().ib(); }
  KG_INLINE const Int3& im() const { return box().im(); }
  KG_INLINE int ib(int d) const { return box().ib(d); }
  KG_INLINE int im(int d) const { return box().im(d); }
  KG_INLINE int n_cells() const { return box().size(); }
  KG_INLINE int n_comps() const { return n_comps_; }
  KG_INLINE int size() const { return n_comps() * n_cells(); }

  const_pointer data() const { return storage().data(); }
  pointer data() { return storage().data(); }

  KG_INLINE const_reference operator()(int m, int i, int j, int k) const
  {
    return storage()[index(m, {i, j, k})];
  }

  KG_INLINE reference operator()(int m, int i, int j, int k)
  {
    return storage()[index(m, {i, j, k})];
  }

  KG_INLINE int index(int m, Int3 idx) const
  {
#ifdef BOUNDS_CHECK
    assert(m >= 0 && m < n_comps_);
    assert(idx[0] >= ib_[0] && idx[0] < ib_[0] + im_[0]);
    assert(idx[1] >= ib_[1] && idx[1] < ib_[1] + im_[1]);
    assert(idx[2] >= ib_[2] && idx[2] < ib_[2] + im_[2]);
#endif

    return layoutDataOffset<Layout>(n_comps_, im(), m, idx - ib());
  }

  void zero(int m)
  {
    // FIXME, only correct for SOA!!!
    std::memset(&(*this)(m, ib()[0], ib()[1], ib()[2]), 0,
           n_cells() * sizeof(value_type));
  }

  void zero(int mb, int me)
  {
    for (int m = mb; m < me; m++) {
      zero(m);
    }
  }

  void zero() { std::memset(storage().data(), 0, sizeof(value_type) * size()); }

  void set(int m, const_reference val)
  {
    for (int k = ib()[2]; k < ib()[2] + im()[2]; k++) {
      for (int j = ib()[1]; j < ib()[1] + im()[1]; j++) {
        for (int i = ib()[0]; i < ib()[0] + im()[0]; i++) {
          (*this)(m, i, j, k) = val;
        }
      }
    }
  }

  void scale(int m, const_reference val)
  {
    for (int k = ib()[2]; k < ib()[2] + im()[2]; k++) {
      for (int j = ib()[1]; j < ib()[1] + im()[1]; j++) {
        for (int i = ib()[0]; i < ib()[0] + im()[0]; i++) {
          (*this)(m, i, j, k) *= val;
        }
      }
    }
  }

  template <typename F>
  void copy_comp(int mto, const F& from, int mfrom)
  {
    for (int k = ib()[2]; k < ib()[2] + im()[2]; k++) {
      for (int j = ib()[1]; j < ib()[1] + im()[1]; j++) {
        for (int i = ib()[0]; i < ib()[0] + im()[0]; i++) {
          (*this)(mto, i, j, k) = from(mfrom, i, j, k);
        }
      }
    }
  }

  template <typename F>
  void axpy_comp(int m_y, const_reference alpha, const F& x, int m_x)
  {
    for (int k = ib()[2]; k < ib()[2] + im()[2]; k++) {
      for (int j = ib()[1]; j < ib()[1] + im()[1]; j++) {
        for (int i = ib()[0]; i < ib()[0] + im()[0]; i++) {
          (*this)(m_y, i, j, k) += alpha * x(m_x, i, j, k);
        }
      }
    }
  }

  value_type max_comp(int m)
  {
    value_type rv = std::numeric_limits<value_type>::lowest();
    for (int k = ib()[2]; k < ib()[2] + im()[2]; k++) {
      for (int j = ib()[1]; j < ib()[1] + im()[1]; j++) {
        for (int i = ib()[0]; i < ib()[0] + im()[0]; i++) {
          rv = std::max(rv, (*this)(m, i, j, k));
        }
      }
    }
    return rv;
  }

  void dump()
  {
    for (int k = ib()[2]; k < ib()[2] + im()[2]; k++) {
      for (int j = ib()[1]; j < ib()[1] + im()[1]; j++) {
        for (int i = ib()[0]; i < ib()[0] + im()[0]; i++) {
          for (int m = 0; m < n_comps_; m++) {
            mprintf("dump: ijk %d:%d:%d m %d: %g\n", i, j, k, m,
                    (*this)(m, i, j, k));
          }
        }
      }
    }
  }

protected:
  KG_INLINE Storage& storage() { return derived().storageImpl(); }
  KG_INLINE const Storage& storage() const { return derived().storageImpl(); }
  KG_INLINE Derived& derived() { return *static_cast<Derived*>(this); }
  KG_INLINE const Derived& derived() const { return *static_cast<const Derived*>(this); }

private:
  Box3 box_; //> lower bounds and length per direction
  int n_comps_;  // # of components
};

} // namespace kg

#endif
