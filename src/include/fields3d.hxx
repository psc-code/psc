
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

  real_t *data;
  int ib[3], im[3]; //> lower bounds and length per direction
  int nr_comp; //> nr of components
  int first_comp; // first component

  fields3d(const int _ib[3], const int _im[3], int _n_comps, int _first_comp=0,
	   real_t *_data=nullptr)
    : ib{ _ib[0], _ib[1], _ib[2] },
      im{ _im[0], _im[1], _im[2] },
      nr_comp{_n_comps},
      first_comp{_first_comp},
      data(_data)
  {
    if (!data) {
      data = (real_t *) calloc(size(), sizeof(*data));
    }
  }

  void dtor()
  {
    free(data);
    data = NULL;
  }

  real_t  operator()(int m, int i, int j, int k) const { return data[index(m, i, j, k)];  }
  real_t& operator()(int m, int i, int j, int k)       { return data[index(m, i, j, k)];  }

  int index(int m, int i, int j, int k) const;

  int size()
  {
    return nr_comp * im[0] * im[1] * im[2];
  }

  int n_cells()
  {
    return im[0] * im[1] * im[2];
  }

  void zero(int m)
  {
    memset(&(*this)(m, ib[0], ib[1], ib[2]), 0, n_cells() * sizeof(real_t));
  }

  void zero(int mb, int me)
  {
    for (int m = mb; m < me; m++) {
      zero(m);
    }
  }

  void zero()
  {
    memset(data, 0, sizeof(real_t) * size());
  }

  void set(int m, real_t val)
  {
    for (int k = ib[2]; k < ib[2] + im[2]; k++) {
      for (int j = ib[1]; j < ib[1] + im[1]; j++) {
	for (int i = ib[0]; i < ib[0] + im[0]; i++) {
	  (*this)(m, i,j,k) = val;
	}
      }
    }
  }

  void scale(int m, real_t val)
  {
    for (int k = ib[2]; k < ib[2] + im[2]; k++) {
      for (int j = ib[1]; j < ib[1] + im[1]; j++) {
	for (int i = ib[0]; i < ib[0] + im[0]; i++) {
	  (*this)(m, i,j,k) *= val;
	}
      }
    }
  }

  void copy_comp(int mto, const fields3d& from, int mfrom)
  {
    for (int k = ib[2]; k < ib[2] + im[2]; k++) {
      for (int j = ib[1]; j < ib[1] + im[1]; j++) {
	for (int i = ib[0]; i < ib[0] + im[0]; i++) {
	  (*this)(mto, i,j,k) = from(mfrom, i,j,k);
	}
      }
    }
  }

  void axpy_comp(int m_y, real_t alpha, const fields3d& x, int m_x)
  {
    for (int k = ib[2]; k < ib[2] + im[2]; k++) {
      for (int j = ib[1]; j < ib[1] + im[1]; j++) {
	for (int i = ib[0]; i < ib[0] + im[0]; i++) {
	  (*this)(m_y, i,j,k) += alpha * x(m_x, i,j,k);
	}
      }
    }
  }

  real_t max_comp(int m)
  {
    real_t rv = -std::numeric_limits<real_t>::max();
    for (int k = ib[2]; k < ib[2] + im[2]; k++) {
      for (int j = ib[1]; j < ib[1] + im[1]; j++) {
	for (int i = ib[0]; i < ib[0] + im[0]; i++) {
	  rv = std::max(rv, (*this)(m, i,j,k));
	}
      }
    }
    return rv;
  }
};

template<typename R, typename L>
int fields3d<R, L>::index(int m, int i, int j, int k) const
{
#ifdef BOUNDS_CHECK
  assert(m >= first_comp_ && m < n_comp_);
  assert(i >= ib[0] && i < ib[0] + im[0]);
  assert(j >= ib[1] && j < ib[1] + im[1]);
  assert(k >= ib[2] && k < ib[2] + im[2]);
#endif

  if (L::isAOS::value) {
    return (((((k - ib[2])) * im[1] +
	      (j - ib[1])) * im[0] +
	     (i - ib[0])) * nr_comp + m);
  } else {
    return (((((m - first_comp) * im[2] +
	       (k - ib[2])) * im[1] +
	      (j - ib[1])) * im[0] +
	     (i - ib[0])));
  }
}

// ======================================================================
// psc_mfields_base

struct psc_mfields_base
{
  psc_mfields_base(const Grid_t& grid, int n_fields)
    : grid_(grid),
      n_fields_(n_fields)
  {}

  int n_patches() const { return grid_.n_patches(); }

  virtual void zero_comp(int m) = 0;
  virtual void set_comp(int m, double val) = 0;
  virtual void scale_comp(int m, double val) = 0;
  
protected:
  int n_fields_;
  const Grid_t& grid_;
};

// ======================================================================
// psc_mfields_

template<typename F>
struct psc_mfields_ : psc_mfields_base
{
  using fields_t = F;
  using real_t = typename fields_t::real_t;

  psc_mfields_(const Grid_t& grid, int n_fields, int ibn[3], int first_comp)
    : psc_mfields_base(grid, n_fields),
      first_comp_(first_comp)
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

  ~psc_mfields_()
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
    return fields_t(ib, im, n_fields_, first_comp_);
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

  void copy_comp(int mto, psc_mfields_&from, int mfrom)
  {
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].copy_comp(mto, from[p], mfrom);
    }
  }
  
  void axpy_comp(int m_y, real_t alpha, psc_mfields_&x, int m_x)
  {
    for (int p = 0; p < n_patches(); p++) {
      (*this)[p].axpy_comp(m_y, alpha, x[p], m_x);
    }
  }

  double max_comp(int m)
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
private:
  int first_comp_;
};

// ======================================================================
// mfields_base

template<typename S>
struct mfields_base
{
  using sub_t = S;
  using fields_t = typename sub_t::fields_t;
  using real_t = typename fields_t::real_t;

  static_assert(std::is_convertible<sub_t*, psc_mfields_base*>::value,
		"sub classes used in mfields_t must derive from psc_mfields_base");
  
  mfields_base(struct psc_mfields *mflds)
    : mflds_(mflds),
      sub_(mrc_to_subobj(mflds, sub_t))
  {}

  unsigned int n_patches() const { return mflds_->nr_patches; }
  unsigned int n_fields() const { return mflds_->nr_fields; }

  fields_t operator[](int p)
  {
    return (*sub_)[p];
  }

  void put_as(struct psc_mfields *mflds_base, int mb, int me)
  {
    psc_mfields_put_as(mflds_, mflds_base, mb, me);
  }

  struct psc_mfields *mflds() { return mflds_; }
  
  sub_t* operator->() { return sub_; }
  sub_t* sub() { return sub_; }
  
private:
  struct psc_mfields *mflds_;
  sub_t *sub_;
};

using mfields_base_t = mfields_base<psc_mfields_base>;

#endif

