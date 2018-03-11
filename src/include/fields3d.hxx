
#ifndef FIELDS3D_HXX
#define FIELDS3D_HXX

#include <type_traits>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <algorithm>

#include "grid.hxx"
#include "psc_fields.h"

template<bool AOS>
struct Layout
{
  using isAOS = std::integral_constant<bool, AOS>;
};

using LayoutAOS = Layout<true>;
using LayoutSOA = Layout<false>;

// ======================================================================
// fields3d

template<typename R, typename L=LayoutSOA>
struct fields3d {
  using real_t = R;
  using layout = L;

  fields3d(Int3 ib, Int3 im, int n_comps, real_t* data=nullptr)
    : ib_{ib}, im_{im},
      n_comps_{n_comps},
      data_{data}
  {
    if (!data_) {
      data_ = (real_t *) calloc(size(), sizeof(*data_));
    }
  }

  void dtor()
  {
    free(data_);
    data_ = NULL;
  }

  real_t  operator()(int m, int i, int j, int k) const { return data_[index(m, i, j, k)];  }
  real_t& operator()(int m, int i, int j, int k)       { return data_[index(m, i, j, k)];  }

  real_t* data() { return data_; }
  int index(int m, int i, int j, int k) const;
  int n_comps() const { return n_comps_; }
  int n_cells() const { return im_[0] * im_[1] * im_[2]; }
  int size()    const { return n_comps() * n_cells(); }
  Int3 ib()     const { return ib_; }
  Int3 im()     const { return im_; }

  void zero(int m)
  {
    memset(&(*this)(m, ib_[0], ib_[1], ib_[2]), 0, n_cells() * sizeof(real_t));
  }

  void zero(int mb, int me)
  {
    for (int m = mb; m < me; m++) {
      zero(m);
    }
  }

  void zero()
  {
    memset(data_, 0, sizeof(real_t) * size());
  }

  void set(int m, real_t val)
  {
    for (int k = ib_[2]; k < ib_[2] + im_[2]; k++) {
      for (int j = ib_[1]; j < ib_[1] + im_[1]; j++) {
	for (int i = ib_[0]; i < ib_[0] + im_[0]; i++) {
	  (*this)(m, i,j,k) = val;
	}
      }
    }
  }

  void scale(int m, real_t val)
  {
    for (int k = ib_[2]; k < ib_[2] + im_[2]; k++) {
      for (int j = ib_[1]; j < ib_[1] + im_[1]; j++) {
	for (int i = ib_[0]; i < ib_[0] + im_[0]; i++) {
	  (*this)(m, i,j,k) *= val;
	}
      }
    }
  }

  void copy_comp(int mto, const fields3d& from, int mfrom)
  {
    for (int k = ib_[2]; k < ib_[2] + im_[2]; k++) {
      for (int j = ib_[1]; j < ib_[1] + im_[1]; j++) {
	for (int i = ib_[0]; i < ib_[0] + im_[0]; i++) {
	  (*this)(mto, i,j,k) = from(mfrom, i,j,k);
	}
      }
    }
  }

  void axpy_comp(int m_y, real_t alpha, const fields3d& x, int m_x)
  {
    for (int k = ib_[2]; k < ib_[2] + im_[2]; k++) {
      for (int j = ib_[1]; j < ib_[1] + im_[1]; j++) {
	for (int i = ib_[0]; i < ib_[0] + im_[0]; i++) {
	  (*this)(m_y, i,j,k) += alpha * x(m_x, i,j,k);
	}
      }
    }
  }

  real_t max_comp(int m)
  {
    real_t rv = -std::numeric_limits<real_t>::max();
    for (int k = ib_[2]; k < ib_[2] + im_[2]; k++) {
      for (int j = ib_[1]; j < ib_[1] + im_[1]; j++) {
	for (int i = ib_[0]; i < ib_[0] + im_[0]; i++) {
	  rv = std::max(rv, (*this)(m, i,j,k));
	}
      }
    }
    return rv;
  }

  real_t* data_;
  Int3 ib_, im_; //> lower bounds and length per direction
  int n_comps_; // # of components
};

template<typename R, typename L>
int fields3d<R, L>::index(int m, int i, int j, int k) const
{
#ifdef BOUNDS_CHECK
  assert(m >= 0 && m < n_comps_);
  assert(i >= ib_[0] && i < ib_[0] + im_[0]);
  assert(j >= ib_[1] && j < ib_[1] + im_[1]);
  assert(k >= ib_[2] && k < ib_[2] + im_[2]);
#endif

  if (L::isAOS::value) {
    return (((((k - ib_[2])) * im_[1] +
	      (j - ib_[1])) * im_[0] +
	     (i - ib_[0])) * n_comps_ + m);
  } else {
    return (((((m) * im_[2] +
	       (k - ib_[2])) * im_[1] +
	      (j - ib_[1])) * im_[0] +
	     (i - ib_[0])));
  }
}

// ======================================================================
// MfieldsBase

struct MfieldsBase
{
  struct fields_t
  {
    struct real_t {};
  };
  
  MfieldsBase(const Grid_t& grid, int n_fields)
    : grid_(grid),
      n_fields_(n_fields)
  {}

  int n_patches() const { return grid_.n_patches(); }
  int n_comps() const { return n_fields_; }

  virtual void zero_comp(int m) = 0;
  virtual void set_comp(int m, double val) = 0;
  virtual void scale_comp(int m, double val) = 0;
  virtual void axpy_comp(int m_y, double alpha, MfieldsBase& x, int m_x) = 0;
  virtual double max_comp(int m) = 0;

  void zero()            { for (int m = 0; m < n_fields_; m++) zero_comp(m); }
  void scale(double val) { for (int m = 0; m < n_fields_; m++) scale_comp(m, val); }
  void axpy(double alpha, MfieldsBase& x)
  {
    for (int m = 0; m < n_fields_; m++) {
      axpy_comp(m, alpha, x, m);
    }
  }
  
protected:
  int n_fields_;
  const Grid_t& grid_;
};

// ======================================================================
// Mfields

template<typename F>
struct Mfields : MfieldsBase
{
  using fields_t = F;
  using real_t = typename fields_t::real_t;

  Mfields(const Grid_t& grid, int n_fields, int ibn[3])
    : MfieldsBase(grid, n_fields)
  {
    unsigned int size = 1;
    for (int d = 0; d < 3; d++) {
      ib[d] = -ibn[d];
      im[d] = grid_.ldims[d] + 2 * ibn[d];
      size *= im[d];
    }

    data = (real_t**) calloc(n_patches(), sizeof(*data));
    for (int p = 0; p < n_patches(); p++) {
      data[p] = (real_t *) calloc(n_fields * size, sizeof(real_t));
    }
  }

  ~Mfields()
  {
    if (data) { // FIXME, since this object exists without a constructor having been called, for now...
      for (int p = 0; p < n_patches(); p++) {
	free(data[p]);
      }
      free(data);
    }
  }

  fields_t operator[](int p)
  {
    return fields_t(ib, im, n_fields_, data[p]);
  }

  void zero_comp(int m) override
  {
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].zero(m);
    }
  }

  void set_comp(int m, double val) override
  {
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].set(m, val);
    }
  }
  
  void scale_comp(int m, double val) override
  {
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].scale(m, val);
    }
  }

  void copy_comp(int mto, Mfields& from, int mfrom)
  {
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].copy_comp(mto, from[p], mfrom);
    }
  }
  
  void axpy_comp(int m_y, double alpha, MfieldsBase& x_base, int m_x) override
  {
    // FIXME? dynamic_cast would actually be more appropriate
    Mfields& x = static_cast<Mfields&>(x_base);
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].axpy_comp(m_y, alpha, x[p], m_x);
    }
  }

  double max_comp(int m) override
  {
    double rv = -std::numeric_limits<double>::max();
    for (int p = 0; p < n_patches(); p++) {
      rv = std::max(rv, double((*this)[p].max_comp(m)));
    }
    return rv;
  }

  real_t **data;
  int ib[3]; //> lower left corner for each patch (incl. ghostpoints)
  int im[3]; //> extent for each patch (incl. ghostpoints)
};

// ======================================================================
// PscMfields

template<typename S>
struct PscMfields
{
  using sub_t = S;
  using fields_t = typename sub_t::fields_t;
  using real_t = typename fields_t::real_t;

  static_assert(std::is_convertible<sub_t*, MfieldsBase*>::value,
		"sub classes used in mfields_t must derive from psc_mfields_base");
  
  PscMfields(struct psc_mfields *mflds)
    : mflds_(mflds)
  {}

  unsigned int n_fields() const { return mflds_->nr_fields; }

  fields_t operator[](int p)
  {
    return (*sub())[p];
  }

  template<typename MF>
  MF get_as(int mb, int me)
  {
    return mflds_->get_as<MF>(mb, me);
  }
  
  template<typename MF>
  void put_as(MF mflds_base, int mb, int me)
  {
    return put_as(mflds_base.mflds(), mb, me);
  }

  void put_as(struct psc_mfields *mflds_base, int mb, int me)
  {
    psc_mfields_put_as(mflds_, mflds_base, mb, me);
  }

  struct psc_mfields *mflds() { return mflds_; }
  
  sub_t* operator->() { return sub(); }
  sub_t* sub() { return mrc_to_subobj(mflds_, sub_t); }
  
private:
  struct psc_mfields *mflds_;
};

using PscMfieldsBase = PscMfields<MfieldsBase>;

#endif

