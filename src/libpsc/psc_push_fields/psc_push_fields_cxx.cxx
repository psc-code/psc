
#include "psc_push_fields_iface.h"

#include "psc.h"
#include "psc_fields_as_single.h"

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

template<typename R, typename DIM = DIM_XYZ>
class Field3d
{
public:
  using real_t = R;

  Field3d(const fields_t& f)
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

using Field = Field3d<float, DIM_XZ>;

// ----------------------------------------------------------------------

struct params_push_fields {
  fields_real_t dth;
  fields_real_t cnx;
  fields_real_t cny;
  fields_real_t cnz;
  int ldims[3];
};

static struct params_push_fields prm;

// ----------------------------------------------------------------------
// params_push_fields_set

static void
params_push_fields_set(struct psc *psc, double dt_fac)
{
  for (int p = 1; p < psc->nr_patches; p++) {
    for (int d = 0; d < 3; d++) {
      assert(psc->patch[0].ldims[d] == psc->patch[p].ldims[d]);
      assert(psc->patch[0].dx[d] == psc->patch[p].dx[d]);
    }
  }

  prm.dth = dt_fac * psc->dt;

  prm.cnx = prm.dth / psc->patch[0].dx[0];
  prm.cny = prm.dth / psc->patch[0].dx[1];
  prm.cnz = prm.dth / psc->patch[0].dx[2];

  if (psc->domain.gdims[0] == 1) {
    prm.cnx = 0.;
  }
  if (psc->domain.gdims[1] == 1) {
    prm.cny = 0.;
  }
  if (psc->domain.gdims[2] == 1) {
    prm.cnz = 0.;
  }

  for (int d = 0; d < 3; d++) {
    prm.ldims[d] = psc->patch[0].ldims[d];
  }
}

// ----------------------------------------------------------------------
// psc_push_fields_single_push_E_xz

void
psc_push_fields_single_push_E_xz(struct psc_push_fields *push, fields_t flds,
				 struct psc *psc, double dt_fac)
{
  params_push_fields_set(psc, dt_fac);

  Field F(flds);
  MHERE;
  foreach_3d(ppsc, 0, i,j,k, 1, 2) {
    F(EX, i,j,k) +=
      prm.cny * (F(HZ, i,j,k) - F(HZ, i,j-1,k)) -
      prm.cnz * (F(HY, i,j,k) - F(HY, i,j,k-1)) -
      prm.dth * F(JXI, i,j,k);

    F(EY, i,j,k) +=
      prm.cnz * (F(HX, i,j,k) - F(HX, i,j,k-1)) -
      prm.cnx * (F(HZ, i,j,k) - F(HZ, i-1,j,k)) -
      prm.dth * F(JYI, i,j,k);

    F(EZ, i,j,k) +=
      prm.cnx * (F(HY, i,j,k) - F(HY, i-1,j,k)) -
      prm.cny * (F(HX, i,j,k) - F(HX, i,j-1,k)) -
      prm.dth * F(JZI, i,j,k);
  } foreach_3d_end;
}

