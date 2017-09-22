
#include "psc.h"

struct params_push_fields {
  fields_real_t dt;
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
params_push_fields_set(struct psc *psc)
{
  for (int p = 1; p < psc->nr_patches; p++) {
    for (int d = 0; d < 3; d++) {
      assert(psc->patch[0].ldims[d] == psc->patch[p].ldims[d]);
      assert(psc->patch[0].dx[d] == psc->patch[p].dx[d]);
    }
  }

  prm.dt = psc->dt;
  prm.dth = .5f * psc->dt;

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

static void
psc_push_fields_sub_push_E_gen(struct psc_push_fields *push, fields_t flds)
{
  psc_foreach_3d(ppsc, 0, ix, iy, iz, 2, 2) {
    _F3(flds, EX, ix,iy,iz) +=
      prm.cny * (_F3(flds, HZ, ix,iy,iz) - _F3(flds, HZ, ix,iy-1,iz)) -
      prm.cnz * (_F3(flds, HY, ix,iy,iz) - _F3(flds, HY, ix,iy,iz-1)) -
      prm.dth * _F3(flds, JXI, ix,iy,iz);
    
    _F3(flds, EY, ix,iy,iz) +=
      prm.cnz * (_F3(flds, HX, ix,iy,iz) - _F3(flds, HX, ix,iy,iz-1)) -
      prm.cnx * (_F3(flds, HZ, ix,iy,iz) - _F3(flds, HZ, ix-1,iy,iz)) -
      prm.dth * _F3(flds, JYI, ix,iy,iz);
    
    _F3(flds, EZ, ix,iy,iz) +=
      prm.cnx * (_F3(flds, HY, ix,iy,iz) - _F3(flds, HY, ix-1,iy,iz)) -
      prm.cny * (_F3(flds, HX, ix,iy,iz) - _F3(flds, HX, ix,iy-1,iz)) -
      prm.dth * _F3(flds, JZI, ix,iy,iz);
  } foreach_3d_end;
}

static void
psc_push_fields_sub_push_E_yz(struct psc_push_fields *push, fields_t flds)
{
  for (int iz = -1; iz < prm.ldims[2] + 2; iz++) {
    for (int iy = -1; iy < prm.ldims[1] + 2; iy++) {
      _F3(flds, EX, 0,iy,iz) +=
	prm.cny * (_F3(flds, HZ, 0,iy,iz) - _F3(flds, HZ, 0,iy-1,iz)) -
	prm.cnz * (_F3(flds, HY, 0,iy,iz) - _F3(flds, HY, 0,iy,iz-1)) -
	prm.dth * _F3(flds, JXI, 0,iy,iz);
    }

    for (int iy = -2; iy < prm.ldims[1] + 2; iy++) {
      _F3(flds, EY, 0,iy,iz) +=
	prm.cnz * (_F3(flds, HX, 0,iy,iz) - _F3(flds, HX, 0,iy,iz-1)) -
	0 -
	prm.dth * _F3(flds, JYI, 0,iy,iz);
    }
  }
      
  for (int iz = -2; iz < prm.ldims[2] + 2; iz++) {
    for (int iy = -1; iy < prm.ldims[1] + 2; iy++) {
      _F3(flds, EZ, 0,iy,iz) +=
	0 -
	prm.cny * (_F3(flds, HX, 0,iy,iz) - _F3(flds, HX, 0,iy-1,iz)) -
	prm.dth * _F3(flds, JZI, 0,iy,iz);
    }
  }
}

// E-field propagation E^(n)    , H^(n), j^(n) 
//                  -> E^(n+0.5), H^(n), j^(n)
// Ex^{n}[-.5:+.5][-1:1][-1:1] -> Ex^{n+.5}[-.5:+.5][-1:1][-1:1]
// using Hx^{n}[-1:1][-1.5:1.5][-1.5:1.5]
//       jx^{n+1}[-.5:.5][-1:1][-1:1]

static void
psc_push_fields_sub_push_mflds_E(struct psc_push_fields *push, struct psc_mfields *mflds_base)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, JXI, HX + 3);

  params_push_fields_set(ppsc);
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    int *gdims = ppsc->domain.gdims;
    if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
      psc_push_fields_sub_push_E_yz(push, flds);
    } else {
      psc_push_fields_sub_push_E_gen(push, flds);
    }
  }

  psc_mfields_put_as(mflds, mflds_base, EX, EX + 3);
}

static void
psc_push_fields_sub_push_H_gen(struct psc_push_fields *push, fields_t flds)
{
  psc_foreach_3d(ppsc, 0, ix, iy, iz, 1, 2) {
    _F3(flds, HX, ix,iy,iz) -=
      prm.cny * (_F3(flds, EZ, ix,iy+1,iz) - _F3(flds, EZ, ix,iy,iz)) -
      prm.cnz * (_F3(flds, EY, ix,iy,iz+1) - _F3(flds, EY, ix,iy,iz));
    
    _F3(flds, HY, ix,iy,iz) -=
      prm.cnz * (_F3(flds, EX, ix,iy,iz+1) - _F3(flds, EX, ix,iy,iz)) -
      prm.cnx * (_F3(flds, EZ, ix+1,iy,iz) - _F3(flds, EZ, ix,iy,iz));
    
    _F3(flds, HZ, ix,iy,iz) -=
      prm.cnx * (_F3(flds, EY, ix+1,iy,iz) - _F3(flds, EY, ix,iy,iz)) -
      prm.cny * (_F3(flds, EX, ix,iy+1,iz) - _F3(flds, EX, ix,iy,iz));
  } foreach_3d_end;
}

static void
psc_push_fields_sub_push_H_yz(struct psc_push_fields *push, fields_t flds)
{
  for (int iz = -1; iz < prm.ldims[2] + 1; iz++) {
    for (int iy = -1; iy < prm.ldims[1] + 1; iy++) {
      _F3(flds, HX, 0,iy,iz) -=
	prm.cny * (_F3(flds, EZ, 0,iy+1,iz) - _F3(flds, EZ, 0,iy,iz)) -
	prm.cnz * (_F3(flds, EY, 0,iy,iz+1) - _F3(flds, EY, 0,iy,iz));
    }

    for (int iy = -1; iy < prm.ldims[1] + 2; iy++) {
      _F3(flds, HY, 0,iy,iz) -=
	prm.cnz * (_F3(flds, EX, 0,iy,iz+1) - _F3(flds, EX, 0,iy,iz)) -
	0;
    }
  }
      
  for (int iz = -1; iz < prm.ldims[2] + 2; iz++) {
    for (int iy = -1; iy < prm.ldims[1] + 1; iy++) {
      _F3(flds, HZ, 0,iy,iz) -=
	0 -
	prm.cny * (_F3(flds, EX, 0,iy+1,iz) - _F3(flds, EX, 0,iy,iz));
    }
  }
}

// B-field propagation E^(n+0.5), H^(n    ), j^(n), m^(n+0.5)
//                  -> E^(n+0.5), H^(n+0.5), j^(n), m^(n+0.5)
// Hx^{n}[:][-.5:.5][-.5:.5] -> Hx^{n+.5}[:][-.5:.5][-.5:.5]
// using Ex^{n+.5}[-.5:+.5][-1:1][-1:1]

static void
psc_push_fields_sub_push_mflds_H(struct psc_push_fields *push, struct psc_mfields *mflds_base)
{
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, FIELDS_TYPE, EX, HX + 3);

  params_push_fields_set(ppsc);
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_t flds = fields_t_mflds(mflds, p);
    int *gdims = ppsc->domain.gdims;
    if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
      psc_push_fields_sub_push_H_yz(push, flds);
    } else {
      psc_push_fields_sub_push_H_gen(push, flds);
    }
  }

  psc_mfields_put_as(mflds, mflds_base, HX, HX + 3);
}

