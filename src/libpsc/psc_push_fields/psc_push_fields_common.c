
#include "psc.h"

static void
psc_push_fields_sub_push_E_gen(struct psc_push_fields *push, struct psc_fields *flds)
{
  struct psc_patch *patch = &ppsc->patch[flds->p];
  fields_real_t cnx = .5 * ppsc->dt / patch->dx[0];
  fields_real_t cny = .5 * ppsc->dt / patch->dx[1];
  fields_real_t cnz = .5 * ppsc->dt / patch->dx[2];

  if (ppsc->domain.gdims[0] == 1) {
    cnx = 0.;
  }
  if (ppsc->domain.gdims[1] == 1) {
    cny = 0.;
  }
  if (ppsc->domain.gdims[2] == 1) {
    cnz = 0.;
  }

  psc_foreach_3d(ppsc, flds->p, ix, iy, iz, 2, 2) {
    F3(flds, EX, ix,iy,iz) +=
      cny * (F3(flds, HZ, ix,iy,iz) - F3(flds, HZ, ix,iy-1,iz)) -
      cnz * (F3(flds, HY, ix,iy,iz) - F3(flds, HY, ix,iy,iz-1)) -
      .5 * ppsc->dt * F3(flds, JXI, ix,iy,iz);
    
    F3(flds, EY, ix,iy,iz) +=
      cnz * (F3(flds, HX, ix,iy,iz) - F3(flds, HX, ix,iy,iz-1)) -
      cnx * (F3(flds, HZ, ix,iy,iz) - F3(flds, HZ, ix-1,iy,iz)) -
      .5 * ppsc->dt * F3(flds, JYI, ix,iy,iz);
    
    F3(flds, EZ, ix,iy,iz) +=
      cnx * (F3(flds, HY, ix,iy,iz) - F3(flds, HY, ix-1,iy,iz)) -
      cny * (F3(flds, HX, ix,iy,iz) - F3(flds, HX, ix,iy-1,iz)) -
      .5 * ppsc->dt * F3(flds, JZI, ix,iy,iz);
  } foreach_3d_end;
}

static void
psc_push_fields_sub_push_E_yz(struct psc_push_fields *push, struct psc_fields *flds)
{
  struct psc_patch *patch = &ppsc->patch[flds->p];
  fields_real_t cny = .5 * ppsc->dt / patch->dx[1];
  fields_real_t cnz = .5 * ppsc->dt / patch->dx[2];

  int *ldims = patch->ldims;
  for (int iz = -1; iz < ldims[2] + 2; iz++) {
    for (int iy = -1; iy < ldims[1] + 2; iy++) {
      F3(flds, EX, 0,iy,iz) +=
	cny * (F3(flds, HZ, 0,iy,iz) - F3(flds, HZ, 0,iy-1,iz))
	- cnz * (F3(flds, HY, 0,iy,iz) - F3(flds, HY, 0,iy,iz-1))
	- .5 * ppsc->dt * F3(flds, JXI, 0,iy,iz);
    }

    for (int iy = -2; iy < ldims[1] + 2; iy++) {
      F3(flds, EY, 0,iy,iz) +=
	cnz * (F3(flds, HX, 0,iy,iz) - F3(flds, HX, 0,iy,iz-1))
	- .5 * ppsc->dt * F3(flds, JYI, 0,iy,iz);
    }
  }
      
  for (int iz = -2; iz < ldims[2] + 2; iz++) {
    for (int iy = -1; iy < ldims[1] + 2; iy++) {
      F3(flds, EZ, 0,iy,iz) +=
	- cny * (F3(flds, HX, 0,iy,iz) - F3(flds, HX, 0,iy-1,iz))
	- .5 * ppsc->dt * F3(flds, JZI, 0,iy,iz);
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

  for (int p = 0; p < mflds->nr_patches; p++) {
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
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
psc_push_fields_sub_push_H_gen(struct psc_push_fields *push, struct psc_fields *flds)
{
  struct psc_patch *patch = &ppsc->patch[flds->p];
  fields_real_t cnx = .5 * ppsc->dt / patch->dx[0];
  fields_real_t cny = .5 * ppsc->dt / patch->dx[1];
  fields_real_t cnz = .5 * ppsc->dt / patch->dx[2];

  if (ppsc->domain.gdims[0] == 1) {
    cnx = 0.;
  }
  if (ppsc->domain.gdims[1] == 1) {
    cny = 0.;
  }
  if (ppsc->domain.gdims[2] == 1) {
    cnz = 0.;
  }

  psc_foreach_3d(ppsc, flds->p, ix, iy, iz, 1, 2) {
    F3(flds, HX, ix,iy,iz) -=
      cny * (F3(flds, EZ, ix,iy+1,iz) - F3(flds, EZ, ix,iy,iz)) -
      cnz * (F3(flds, EY, ix,iy,iz+1) - F3(flds, EY, ix,iy,iz));
    
    F3(flds, HY, ix,iy,iz) -=
      cnz * (F3(flds, EX, ix,iy,iz+1) - F3(flds, EX, ix,iy,iz)) -
      cnx * (F3(flds, EZ, ix+1,iy,iz) - F3(flds, EZ, ix,iy,iz));
    
    F3(flds, HZ, ix,iy,iz) -=
      cnx * (F3(flds, EY, ix+1,iy,iz) - F3(flds, EY, ix,iy,iz)) -
	cny * (F3(flds, EX, ix,iy+1,iz) - F3(flds, EX, ix,iy,iz));
  } foreach_3d_end;
}

static void
psc_push_fields_sub_push_H_yz(struct psc_push_fields *push, struct psc_fields *flds)
{
  struct psc_patch *patch = &ppsc->patch[flds->p];
  fields_real_t cny = .5 * ppsc->dt / patch->dx[1];
  fields_real_t cnz = .5 * ppsc->dt / patch->dx[2];

  int *ldims = patch->ldims;
  for (int iz = -1; iz < ldims[2] + 1; iz++) {
    for (int iy = -1; iy < ldims[1] + 1; iy++) {
      F3(flds, HX, 0,iy,iz) -=
	  cny * (F3(flds, EZ, 0,iy+1,iz) - F3(flds, EZ, 0,iy,iz))
	- cnz * (F3(flds, EY, 0,iy,iz+1) - F3(flds, EY, 0,iy,iz));
    }

    for (int iy = -1; iy < ldims[1] + 2; iy++) {
      F3(flds, HY, 0,iy,iz) -=
	  cnz * (F3(flds, EX, 0,iy,iz+1) - F3(flds, EX, 0,iy,iz));
    }
  }
      
  for (int iz = -1; iz < ldims[2] + 2; iz++) {
    for (int iy = -1; iy < ldims[1] + 1; iy++) {
      F3(flds, HZ, 0,iy,iz) -=
	- cny * (F3(flds, EX, 0,iy+1,iz) - F3(flds, EX, 0,iy,iz));
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

  for (int p = 0; p < mflds->nr_patches; p++) {
    int *gdims = ppsc->domain.gdims;
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    if (gdims[0] == 1 && gdims[1] > 1 && gdims[2] > 1) {
      psc_push_fields_sub_push_H_yz(push, flds);
    } else {
      psc_push_fields_sub_push_H_gen(push, flds);
    }
  }

  psc_mfields_put_as(mflds, mflds_base, HX, HX + 3);
}

