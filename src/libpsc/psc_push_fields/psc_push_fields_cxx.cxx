
#include "psc_push_fields_iface.h"

#include "psc.h"

#include <type_traits>

template<bool IX = false, bool IY = false, bool IZ = false>
struct Invar
{
  using InvarX = std::integral_constant<bool, IX>;
  using InvarY = std::integral_constant<bool, IY>;
  using InvarZ = std::integral_constant<bool, IZ>;
};

using DIM_XYZ = Invar<false, false, false>;
using DIM_XZ = Invar<false, true, false>;

template<typename R, typename F, typename D = DIM_XYZ>
class Fields3d
{
public:
  using real_t = R;
  using fields_t = F;
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
class PushE
{
public:
  using real_t = typename Fields::real_t;
  using fields_t = typename Fields::fields_t;
  using DIM = typename Fields::DIM;
  
  PushE(const fields_t& flds, struct psc* psc, double dt_fac)
    : F(flds)
  {
    for (int p = 1; p < psc->nr_patches; p++) {
      for (int d = 0; d < 3; d++) {
	assert(psc->patch[0].ldims[d] == psc->patch[p].ldims[d]);
	assert(psc->patch[0].dx[d] == psc->patch[p].dx[d]);
      }
    }
    
    dth = dt_fac * psc->dt;
    
    cnx = DIM::InvarX::value ? 0 : dth / psc->patch[0].dx[0];
    cny = DIM::InvarY::value ? 0 : dth / psc->patch[0].dx[1];
    cnz = DIM::InvarZ::value ? 0 : dth / psc->patch[0].dx[2];
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

private:
  Fields F;
  real_t dth;
  real_t cnx, cny, cnz;
};

// ----------------------------------------------------------------------
// psc_push_fields_push_E

template<typename Fields>
void psc_push_fields_push_E(struct psc_push_fields *push, typename Fields::fields_t flds,
			    struct psc *psc, double dt_fac)
{
  PushE<Fields> push_E(flds, psc, dt_fac);

  MHERE;
  Foreach_3d(push_E, 1, 2);
}

void
psc_push_fields_single_push_E_xz(struct psc_push_fields *push, fields_single_t flds,
				 struct psc *psc, double dt_fac)
{
  using Fields = Fields3d<fields_single_real_t, fields_single_t, DIM_XZ>;
  psc_push_fields_push_E<Fields>(push, flds, psc, dt_fac);
}

