
#ifndef PSC_PUSH_FIELS_IMPL_HXX
#define PSC_PUSH_FIELS_IMPL_HXX

#include "psc.h" // FIXME, for foreach_3d macro

#include <type_traits>

template<bool IX = false, bool IY = false, bool IZ = false>
struct Invar
{
  using InvarX = std::integral_constant<bool, IX>;
  using InvarY = std::integral_constant<bool, IY>;
  using InvarZ = std::integral_constant<bool, IZ>;
};

using DIM_XYZ = Invar<false, false, false>;
using DIM_XY  = Invar<false, false, true >;
using DIM_XZ  = Invar<false, true , false>;
using DIM_YZ  = Invar<true , false, false>;

template<typename F>
struct fields_traits
{
};

template<typename F, typename D = DIM_XYZ>
class Fields3d
{
public:
  using fields_t = F;
  using real_t = typename fields_traits<F>::real_t;
  using DIM = D;

  Fields3d(const fields_t& f)
    : data_(f.data),
      n_comp_(f.nr_comp),
      first_comp_(f.first_comp)
  {
    for (int d = 0; d < 3; d++) {
      ib[d] = f.ib[d];
      im[d] = f.im[d];
    }
  }

  const real_t operator()(int m, int i, int j, int k) const
  {
    return data_[index(m, i, j, k)];
  }

  real_t& operator()(int m, int i, int j, int k)
  {
    return data_[index(m, i, j, k)];
  }

private:
  int index(int m, int i_, int j_, int k_)
  {
    int i = DIM::InvarX::value ? 0 : i_;
    int j = DIM::InvarY::value ? 0 : j_;
    int k = DIM::InvarZ::value ? 0 : k_;

#ifdef BOUNDS_CHECK
    assert(m >= first_comp_ && m < n_comp_);
    assert(i >= ib[0] && i < ib[0] + im[0]);
    assert(j >= ib[1] && j < ib[1] + im[1]);
    assert(k >= ib[2] && k < ib[2] + im[2]);
#endif
  
      return ((((((m) - first_comp_)
	       * im[2] + (k - ib[2]))
	      * im[1] + (j - ib[1]))
	     * im[0] + (i - ib[0])));
      }

private:
  real_t *data_;
  int ib[3], im[3];
  int n_comp_;
  int first_comp_;
};

#include "psc_fields_single.h"
#include "psc_fields_c.h"

template<>
struct fields_traits<fields_single_t>
{
  using real_t = fields_single_real_t;
  static constexpr const char* name = "single";
};

template<>
struct fields_traits<fields_c_t>
{
  using real_t = fields_c_real_t;
  static constexpr const char* name = "c";
};

// ----------------------------------------------------------------------
// Foreach_3d

template<class F>
static void Foreach_3d(F f, int l, int r)
{
  foreach_3d(ppsc, 0, i,j,k, l, r) {
    f.x(i,j,k);
    f.y(i,j,k);
    f.z(i,j,k);
  } foreach_3d_end;
}

// ----------------------------------------------------------------------

template<typename Fields>
class PushBase
{
public:
  using real_t = typename Fields::real_t;
  using fields_t = typename Fields::fields_t;
  using DIM = typename Fields::DIM;
  
  PushBase(struct psc* psc, double dt_fac)
  {
    for (int p = 1; p < psc->nr_patches; p++) {
      for (int d = 0; d < 3; d++) {
	assert(psc->patch[0].ldims[d] == psc->patch[p].ldims[d]);
	assert(psc->patch[0].dx[d] == psc->patch[p].dx[d]);
      }
    }
    
    dth = dt_fac * psc->dt;

    // FIXME, it'd be even better to not even calculate derivates
    // that will be multiplied by 0 
    cnx = DIM::InvarX::value ? 0 : dth / psc->patch[0].dx[0];
    cny = DIM::InvarY::value ? 0 : dth / psc->patch[0].dx[1];
    cnz = DIM::InvarZ::value ? 0 : dth / psc->patch[0].dx[2];
  }

protected:
  real_t dth;
  real_t cnx, cny, cnz;
};
  
template<typename Fields>
class PushE : PushBase<Fields>
{
public:
  using Base = PushBase<Fields>;
  using typename Base::real_t;
  using typename Base::fields_t;
  
  PushE(const fields_t& flds, struct psc* psc, double dt_fac)
    : Base(psc, dt_fac),
      F(flds)
  {
  }
  
  void x(int i, int j,int k)
  {
    F(EX, i,j,k) += (cny * (F(HZ, i,j,k) - F(HZ, i,j-1,k)) - cnz * (F(HY, i,j,k) - F(HY, i,j,k-1)) -
    		     dth * F(JXI, i,j,k));
  }

  void y(int i, int j, int k)
  {
    F(EY, i,j,k) += (cnz * (F(HX, i,j,k) - F(HX, i,j,k-1)) - cnx * (F(HZ, i,j,k) - F(HZ, i-1,j,k)) -
		     dth * F(JYI, i,j,k));
  }

  void z(int i, int j, int k)
  {
    F(EZ, i,j,k) += (cnx * (F(HY, i,j,k) - F(HY, i-1,j,k)) - cny * (F(HX, i,j,k) - F(HX, i,j-1,k)) -
		     dth * F(JZI, i,j,k));
  }

protected:
  Fields F;
  using Base::dth;
  using Base::cnx;
  using Base::cny;
  using Base::cnz;
};

template<typename Fields>
class PushH : PushBase<Fields>
{
public:
  using Base = PushBase<Fields>;
  using typename Base::real_t;
  using typename Base::fields_t;
  
  PushH(const fields_t& flds, struct psc* psc, double dt_fac)
    : Base(psc, dt_fac),
      F(flds)
  {
  }
  
  void x(int i, int j,int k)
  {
    F(HX, i,j,k) -= (cny * (F(EZ, i,j+1,k) - F(EZ, i,j,k)) - cnz * (F(EY, i,j,k+1) - F(EY, i,j,k)));
  }

  void y(int i, int j, int k)
  {
    F(HY, i,j,k) -= (cnz * (F(EX, i,j,k+1) - F(EX, i,j,k)) - cnx * (F(EZ, i+1,j,k) - F(EZ, i,j,k)));
  }

  void z(int i, int j, int k)
  {
    F(HZ, i,j,k) -= (cnx * (F(EY, i+1,j,k) - F(EY, i,j,k)) - cny * (F(EX, i,j+1,k) - F(EX, i,j,k)));
  }

protected:
  Fields F;
  using Base::dth;
  using Base::cnx;
  using Base::cny;
  using Base::cnz;
};


// ----------------------------------------------------------------------
// psc_push_fields_push_E

template<typename DIM, typename fields_t>
void psc_push_fields_push_E(struct psc_push_fields* push, fields_t flds,
			    struct psc *psc, double dt_fac)
{
  using Fields = Fields3d<fields_t, DIM>;
  PushE<Fields> push_E(flds, psc, dt_fac);

  Foreach_3d(push_E, 1, 2);
}

// ----------------------------------------------------------------------
// psc_push_fields_push_H

template<typename DIM, typename fields_t>
void psc_push_fields_push_H(struct psc_push_fields* push, fields_t flds,
			    struct psc *psc, double dt_fac)
{
  using Fields = Fields3d<fields_t, DIM>;
  PushH<Fields> push_H(flds, psc, dt_fac);

  Foreach_3d(push_H, 2, 1);
}

// ----------------------------------------------------------------------
// psc_push_fields_sub_push_mflds_E
//
// E-field propagation E^(n)    , H^(n), j^(n) 
//                  -> E^(n+0.5), H^(n), j^(n)
// Ex^{n}[-.5:+.5][-1:1][-1:1] -> Ex^{n+.5}[-.5:+.5][-1:1][-1:1]
// using Hx^{n}[-1:1][-1.5:1.5][-1.5:1.5]
//       jx^{n+1}[-.5:.5][-1:1][-1:1]


template<typename fields_t>
static void psc_push_fields_sub_push_mflds_E(struct psc_push_fields *push,
					     struct psc_mfields *mflds_base,
					     double dt_fac)
{
  const char* fields_type = fields_traits<fields_t>::name;

  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, fields_type, JXI, HX + 3);

  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    int *gdims = ppsc->domain.gdims;
    if (gdims[0] > 1 && gdims[1] > 1 && gdims[2] > 1) {
      psc_push_fields_push_E<DIM_XYZ>(push, flds, ppsc, dt_fac);
    } else if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
      psc_push_fields_push_E<DIM_YZ>(push, flds, ppsc, dt_fac);
    } else if (gdims[0] > 1 && gdims[1] == 1 && gdims[2] > 1) {
      psc_push_fields_push_E<DIM_XZ>(push, flds, ppsc, dt_fac);
    } else {
      assert(0);
    }
  }

  psc_mfields_put_as(mflds, mflds_base, EX, EX + 3);
}

// ----------------------------------------------------------------------
// psc_push_fields_sub_push_mflds_H
//
// B-field propagation E^(n+0.5), H^(n    ), j^(n), m^(n+0.5)
//                  -> E^(n+0.5), H^(n+0.5), j^(n), m^(n+0.5)
// Hx^{n}[:][-.5:.5][-.5:.5] -> Hx^{n+.5}[:][-.5:.5][-.5:.5]
// using Ex^{n+.5}[-.5:+.5][-1:1][-1:1]

template<typename fields_t>
static void psc_push_fields_sub_push_mflds_H(struct psc_push_fields *push,
					     struct psc_mfields *mflds_base,
					     double dt_fac)
{
  const char* fields_type = fields_traits<fields_t>::name;

  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, fields_type, EX, HX + 3);

  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    int *gdims = ppsc->domain.gdims;
    if (gdims[0] > 1 && gdims[1] > 1 && gdims[2] > 1) {
      psc_push_fields_push_H<DIM_XYZ>(push, flds, ppsc, dt_fac);
    } else if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
      psc_push_fields_push_H<DIM_YZ>(push, flds, ppsc, dt_fac);
    } else if (gdims[0] > 1 && gdims[1] == 1 && gdims[2] > 1) {
      psc_push_fields_push_H<DIM_XZ>(push, flds, ppsc, dt_fac);
    } else {
      assert(0);
    }
  }

  psc_mfields_put_as(mflds, mflds_base, HX, HX + 3);
}

#endif

