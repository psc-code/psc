
#include "psc.h"

static void
psc_push_fields_c_push_a_E(struct psc_push_fields *push, struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, JXI, HX + 3);
  
  f_real cnx = .5 * ppsc->dt / ppsc->dx[0];
  f_real cny = .5 * ppsc->dt / ppsc->dx[1];
  f_real cnz = .5 * ppsc->dt / ppsc->dx[2];

  if (ppsc->domain.gdims[0] == 1) {
    cnx = 0.;
  }
  if (ppsc->domain.gdims[1] == 1) {
    cny = 0.;
  }
  if (ppsc->domain.gdims[2] == 1) {
    cnz = 0.;
  }

  // E-field propagation E^(n)    , H^(n), j^(n) 
  //                  -> E^(n+0.5), H^(n), j^(n)

  psc_foreach_3d(ppsc, flds->p, ix, iy, iz, 1, 1) {
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

  psc_fields_put_as(flds, flds_base, EX, EX + 3);
}

static void
psc_push_fields_c_push_a_H(struct psc_push_fields *push, struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, HX + 3);
  
  f_real cnx = .5 * ppsc->dt / ppsc->dx[0];
  f_real cny = .5 * ppsc->dt / ppsc->dx[1];
  f_real cnz = .5 * ppsc->dt / ppsc->dx[2];

  if (ppsc->domain.gdims[0] == 1) {
    cnx = 0.;
  }
  if (ppsc->domain.gdims[1] == 1) {
    cny = 0.;
  }
  if (ppsc->domain.gdims[2] == 1) {
    cnz = 0.;
  }

  // B-field propagation E^(n+0.5), H^(n    ), j^(n), m^(n+0.5)
  //                  -> E^(n+0.5), H^(n+0.5), j^(n), m^(n+0.5)

  psc_foreach_3d(ppsc, flds->p, ix, iy, iz, 1, 1) {
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

  psc_fields_put_as(flds, flds_base, HX, HX + 3);
}

static void
psc_push_fields_c_push_b_H(struct psc_push_fields *push, struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, EX, HX + 3);
  
  f_real cnx = .5 * ppsc->dt / ppsc->dx[0];
  f_real cny = .5 * ppsc->dt / ppsc->dx[1];
  f_real cnz = .5 * ppsc->dt / ppsc->dx[2];

  if (ppsc->domain.gdims[0] == 1) {
    cnx = 0.;
  }
  if (ppsc->domain.gdims[1] == 1) {
    cny = 0.;
  }
  if (ppsc->domain.gdims[2] == 1) {
    cnz = 0.;
  }

  // B-field propagation E^(n+0.5), B^(n+0.5), j^(n+1.0), m^(n+0.5)
  //                  -> E^(n+0.5), B^(n+1.0), j^(n+1.0), m^(n+0.5)

  psc_foreach_3d(ppsc, flds->p, ix, iy, iz, 1, 1) {
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

  psc_fields_put_as(flds, flds_base, HX, HX + 3);
}

static void
psc_push_fields_c_push_b_E(struct psc_push_fields *push, struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, FIELDS_TYPE, JXI, HX + 3);
  
  f_real cnx = .5 * ppsc->dt / ppsc->dx[0];
  f_real cny = .5 * ppsc->dt / ppsc->dx[1];
  f_real cnz = .5 * ppsc->dt / ppsc->dx[2];

  if (ppsc->domain.gdims[0] == 1) {
    cnx = 0.;
  }
  if (ppsc->domain.gdims[1] == 1) {
    cny = 0.;
  }
  if (ppsc->domain.gdims[2] == 1) {
    cnz = 0.;
  }

  // E-field propagation E^(n+0.5), B^(n+1.0), j^(n+1.0) 
  //                  -> E^(n+1.0), B^(n+1.0), j^(n+1.0)

  psc_foreach_3d(ppsc, flds->p, ix, iy, iz, 1, 1) {
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

  psc_fields_put_as(flds, flds_base, EX, EX + 3);
}

