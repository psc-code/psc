
#include "psc_push_fields_iface.h"

#include "psc.h"
#include "psc_fields_as_single.h"

#define F3(flds, m, i,j,k) _F3(flds, m, i,0,k)

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
  MHERE;
  foreach_3d(ppsc, 0, i,j,k, 1, 2) {
    F3(flds, EX, i,j,k) +=
      prm.cny * (F3(flds, HZ, i,j,k) - F3(flds, HZ, i,j-1,k)) -
      prm.cnz * (F3(flds, HY, i,j,k) - F3(flds, HY, i,j,k-1)) -
      prm.dth * F3(flds, JXI, i,j,k);

    F3(flds, EY, i,j,k) +=
      prm.cnz * (F3(flds, HX, i,j,k) - F3(flds, HX, i,j,k-1)) -
      prm.cnx * (F3(flds, HZ, i,j,k) - F3(flds, HZ, i-1,j,k)) -
      prm.dth * F3(flds, JYI, i,j,k);

    F3(flds, EZ, i,j,k) +=
      prm.cnx * (F3(flds, HY, i,j,k) - F3(flds, HY, i-1,j,k)) -
      prm.cny * (F3(flds, HX, i,j,k) - F3(flds, HX, i,j-1,k)) -
      prm.dth * F3(flds, JZI, i,j,k);
  } foreach_3d_end;
}

