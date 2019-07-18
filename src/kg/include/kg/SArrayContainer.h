
#ifndef KG_SARRAY_CONTAINER_H
#define KG_SARRAY_CONTAINER_H

#include <kg/Array3d.h>

// FIXME, do noexcept?

namespace kg
{

template <typename C>
struct SArrayContainerInnerTypes;

// ======================================================================
// fields3d_container

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

  SArrayContainer(Int3 ib, Int3 im, int n_comps)
    : ib_{ib}, im_{im}, n_comps_{n_comps}
  {}

  Int3 ib() const { return ib_; }
  Int3 im() const { return im_; }
  int n_cells() const { return im_[0] * im_[1] * im_[2]; }
  int n_comps() const { return n_comps_; }
  int size() const { return n_comps() * n_cells(); }

  const_pointer data() const { return storage().data(); }
  pointer data() { return storage().data(); }

  const_reference operator()(int m, int i, int j, int k) const
  {
    return storage()[index(m, i, j, k)];
  }
  reference operator()(int m, int i, int j, int k)
  {
    return storage()[index(m, i, j, k)];
  }

  void zero(int m)
  {
    // FIXME, only correct for SOA!!!
    memset(&(*this)(m, ib_[0], ib_[1], ib_[2]), 0,
           n_cells() * sizeof(value_type));
  }

  void zero(int mb, int me)
  {
    for (int m = mb; m < me; m++) {
      zero(m);
    }
  }

  void zero() { memset(storage().data(), 0, sizeof(value_type) * size()); }

  int index(int m, int i, int j, int k) const
  {
#ifdef BOUNDS_CHECK
    assert(m >= 0 && m < n_comps_);
    assert(i >= ib_[0] && i < ib_[0] + im_[0]);
    assert(j >= ib_[1] && j < ib_[1] + im_[1]);
    assert(k >= ib_[2] && k < ib_[2] + im_[2]);
#endif

    if (Layout::isAOS::value) {
      return (
        ((((k - ib_[2])) * im_[1] + (j - ib_[1])) * im_[0] + (i - ib_[0])) *
          n_comps_ +
        m);
    } else {
      return (((((m)*im_[2] + (k - ib_[2])) * im_[1] + (j - ib_[1])) * im_[0] +
               (i - ib_[0])));
    }
  }

  void set(int m, const_reference val)
  {
    for (int k = ib_[2]; k < ib_[2] + im_[2]; k++) {
      for (int j = ib_[1]; j < ib_[1] + im_[1]; j++) {
        for (int i = ib_[0]; i < ib_[0] + im_[0]; i++) {
          (*this)(m, i, j, k) = val;
        }
      }
    }
  }

  void scale(int m, const_reference val)
  {
    for (int k = ib_[2]; k < ib_[2] + im_[2]; k++) {
      for (int j = ib_[1]; j < ib_[1] + im_[1]; j++) {
        for (int i = ib_[0]; i < ib_[0] + im_[0]; i++) {
          (*this)(m, i, j, k) *= val;
        }
      }
    }
  }

  template <typename F>
  void copy_comp(int mto, const F& from, int mfrom)
  {
    for (int k = ib_[2]; k < ib_[2] + im_[2]; k++) {
      for (int j = ib_[1]; j < ib_[1] + im_[1]; j++) {
        for (int i = ib_[0]; i < ib_[0] + im_[0]; i++) {
          (*this)(mto, i, j, k) = from(mfrom, i, j, k);
        }
      }
    }
  }

  template <typename F>
  void axpy_comp(int m_y, const_reference alpha, const F& x, int m_x)
  {
    for (int k = ib_[2]; k < ib_[2] + im_[2]; k++) {
      for (int j = ib_[1]; j < ib_[1] + im_[1]; j++) {
        for (int i = ib_[0]; i < ib_[0] + im_[0]; i++) {
          (*this)(m_y, i, j, k) += alpha * x(m_x, i, j, k);
        }
      }
    }
  }

  value_type max_comp(int m)
  {
    value_type rv = std::numeric_limits<value_type>::lowest();
    for (int k = ib_[2]; k < ib_[2] + im_[2]; k++) {
      for (int j = ib_[1]; j < ib_[1] + im_[1]; j++) {
        for (int i = ib_[0]; i < ib_[0] + im_[0]; i++) {
          rv = std::max(rv, (*this)(m, i, j, k));
        }
      }
    }
    return rv;
  }

  void dump()
  {
    for (int k = ib_[2]; k < ib_[2] + im_[2]; k++) {
      for (int j = ib_[1]; j < ib_[1] + im_[1]; j++) {
        for (int i = ib_[0]; i < ib_[0] + im_[0]; i++) {
          for (int m = 0; m < n_comps_; m++) {
            mprintf("dump: ijk %d:%d:%d m %d: %g\n", i, j, k, m,
                    (*this)(m, i, j, k));
          }
        }
      }
    }
  }

protected:
  Storage& storage() { return derived().storageImpl(); }

  const Storage& storage() const { return derived().storageImpl(); }

  Derived& derived() { return *static_cast<Derived*>(this); }

  const Derived& derived() const { return *static_cast<const Derived*>(this); }

private:
  const Int3 ib_, im_; //> lower bounds and length per direction
  const int n_comps_;  // # of components
};

} // namespace kg

#endif
