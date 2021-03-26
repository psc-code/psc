
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
#ifdef BOUNDS_CHECK
    assert(m >= 0 && m < n_comps_);
    assert(i >= ib()[0] && i < ib()[0] + im()[0]);
    assert(j >= ib()[1] && j < ib()[1] + im()[1]);
    assert(k >= ib()[2] && k < ib()[2] + im()[2]);
#endif

    return storage()(i - ib()[0], j - ib()[1], k - ib()[2], m);
  }

  KG_INLINE reference operator()(int m, int i, int j, int k)
  {
#if defined(BOUNDS_CHECK) && !defined(__CUDACC__)
    assert(m >= 0 && m < n_comps_);
    assert(i >= ib()[0] && i < ib()[0] + im()[0]);
    assert(j >= ib()[1] && j < ib()[1] + im()[1]);
    assert(k >= ib()[2] && k < ib()[2] + im()[2]);
#endif

    return storage()(i - ib()[0], j - ib()[1], k - ib()[2], m);
  }

  void zero(int m) { storage().view(_all, _all, _all, m) = value_type(); }

  void zero(int mb, int me)
  {
    storage().view(_all, _all, _all, _s(mb, me)) = value_type();
  }

  void zero() { storage() = value_type(); }

  void set(int m, const_reference val)
  {
    storage().view(_all, _all, _all, m) = val;
  }

  void scale(int m, const_reference val)
  {
    storage().view(_all, _all, _all, m) =
      storage().view(_all, _all, _all, m) * val;
  }

  template <typename F>
  void copy_comp(int mto, const F& from, int mfrom)
  {
    storage().view(_all, _all, _all, mto) = from.view(_all, _all, _all, mfrom);
  }

  template <typename F>
  void axpy_comp(int m_y, const_reference alpha, const F& x, int m_x)
  {
    storage().view(_all, _all, _all, m_y) =
      storage().view(_all, _all, _all, m_y) +
      alpha * x.storage().view(_all, _all, _all, m_x);
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
  KG_INLINE const Derived& derived() const
  {
    return *static_cast<const Derived*>(this);
  }

private:
  Box3 box_;    //> lower bounds and length per direction
  int n_comps_; // # of components
};

} // namespace kg

#endif
