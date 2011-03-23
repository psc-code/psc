
#include "psc_push_fields_private.h"

#include "psc.h"

static void
psc_push_fields_c_push_a_e(struct psc_push_fields *push, mfields_base_t *flds)
{
  assert(psc.domain.bnd_fld_lo[0] == BND_FLD_PERIODIC);
  assert(psc.domain.bnd_fld_lo[1] == BND_FLD_PERIODIC);
  assert(psc.domain.bnd_fld_lo[2] == BND_FLD_PERIODIC);
  assert(psc.domain.bnd_fld_hi[0] == BND_FLD_PERIODIC);
  assert(psc.domain.bnd_fld_hi[1] == BND_FLD_PERIODIC);
  assert(psc.domain.bnd_fld_hi[2] == BND_FLD_PERIODIC);

  f_real cnx = .5 * psc.dt / psc.dx[0];
  f_real cny = .5 * psc.dt / psc.dx[1];
  f_real cnz = .5 * psc.dt / psc.dx[2];

  if (psc.domain.gdims[0] == 1) {
    cnx = 0.;
  }
  if (psc.domain.gdims[1] == 1) {
    cny = 0.;
  }
  if (psc.domain.gdims[2] == 1) {
    cnz = 0.;
  }

  // E-field propagation E^(n)    , H^(n), j^(n) 
  //                  -> E^(n+0.5), H^(n), j^(n)

  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    foreach_3d(p, ix, iy, iz, 1, 1) {
      F3_BASE(pf, EX, ix,iy,iz) +=
	cny * (F3_BASE(pf, HZ, ix,iy,iz) - F3_BASE(pf, HZ, ix,iy-1,iz)) -
	cnz * (F3_BASE(pf, HY, ix,iy,iz) - F3_BASE(pf, HY, ix,iy,iz-1)) -
	.5 * psc.dt * F3_BASE(pf, JXI, ix,iy,iz);
      
      F3_BASE(pf, EY, ix,iy,iz) +=
	cnz * (F3_BASE(pf, HX, ix,iy,iz) - F3_BASE(pf, HX, ix,iy,iz-1)) -
	cnx * (F3_BASE(pf, HZ, ix,iy,iz) - F3_BASE(pf, HZ, ix-1,iy,iz)) -
	.5 * psc.dt * F3_BASE(pf, JYI, ix,iy,iz);
      
      F3_BASE(pf, EZ, ix,iy,iz) +=
	cnx * (F3_BASE(pf, HY, ix,iy,iz) - F3_BASE(pf, HY, ix-1,iy,iz)) -
	cny * (F3_BASE(pf, HX, ix,iy,iz) - F3_BASE(pf, HX, ix,iy-1,iz)) -
	.5 * psc.dt * F3_BASE(pf, JZI, ix,iy,iz);
    } foreach_3d_end;
  }
}

static void
psc_push_fields_c_push_a_h(struct psc_push_fields *push, mfields_base_t *flds)
{
  f_real cnx = .5 * psc.dt / psc.dx[0];
  f_real cny = .5 * psc.dt / psc.dx[1];
  f_real cnz = .5 * psc.dt / psc.dx[2];

  if (psc.domain.gdims[0] == 1) {
    cnx = 0.;
  }
  if (psc.domain.gdims[1] == 1) {
    cny = 0.;
  }
  if (psc.domain.gdims[2] == 1) {
    cnz = 0.;
  }

  // B-field propagation E^(n+0.5), H^(n    ), j^(n), m^(n+0.5)
  //                  -> E^(n+0.5), H^(n+0.5), j^(n), m^(n+0.5)

  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    foreach_3d(p, ix, iy, iz, 1, 1) {
      F3_BASE(pf, HX, ix,iy,iz) -=
	cny * (F3_BASE(pf, EZ, ix,iy+1,iz) - F3_BASE(pf, EZ, ix,iy,iz)) -
	cnz * (F3_BASE(pf, EY, ix,iy,iz+1) - F3_BASE(pf, EY, ix,iy,iz));
      
      F3_BASE(pf, HY, ix,iy,iz) -=
	cnz * (F3_BASE(pf, EX, ix,iy,iz+1) - F3_BASE(pf, EX, ix,iy,iz)) -
	cnx * (F3_BASE(pf, EZ, ix+1,iy,iz) - F3_BASE(pf, EZ, ix,iy,iz));
      
      F3_BASE(pf, HZ, ix,iy,iz) -=
	cnx * (F3_BASE(pf, EY, ix+1,iy,iz) - F3_BASE(pf, EY, ix,iy,iz)) -
	cny * (F3_BASE(pf, EX, ix,iy+1,iz) - F3_BASE(pf, EX, ix,iy,iz));
    } foreach_3d_end;
  }
}

static void
psc_push_fields_c_push_b_h(struct psc_push_fields *push, mfields_base_t *flds)
{
  f_real cnx = .5 * psc.dt / psc.dx[0];
  f_real cny = .5 * psc.dt / psc.dx[1];
  f_real cnz = .5 * psc.dt / psc.dx[2];

  if (psc.domain.gdims[0] == 1) {
    cnx = 0.;
  }
  if (psc.domain.gdims[1] == 1) {
    cny = 0.;
  }
  if (psc.domain.gdims[2] == 1) {
    cnz = 0.;
  }

  // B-field propagation E^(n+0.5), B^(n+0.5), j^(n+1.0), m^(n+0.5)
  //                  -> E^(n+0.5), B^(n+1.0), j^(n+1.0), m^(n+0.5)

  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    foreach_3d(p, ix, iy, iz, 1, 1) {
      F3_BASE(pf, HX, ix,iy,iz) -=
	cny * (F3_BASE(pf, EZ, ix,iy+1,iz) - F3_BASE(pf, EZ, ix,iy,iz)) -
	cnz * (F3_BASE(pf, EY, ix,iy,iz+1) - F3_BASE(pf, EY, ix,iy,iz));
      
      F3_BASE(pf, HY, ix,iy,iz) -=
	cnz * (F3_BASE(pf, EX, ix,iy,iz+1) - F3_BASE(pf, EX, ix,iy,iz)) -
	cnx * (F3_BASE(pf, EZ, ix+1,iy,iz) - F3_BASE(pf, EZ, ix,iy,iz));
      
      F3_BASE(pf, HZ, ix,iy,iz) -=
	cnx * (F3_BASE(pf, EY, ix+1,iy,iz) - F3_BASE(pf, EY, ix,iy,iz)) -
	cny * (F3_BASE(pf, EX, ix,iy+1,iz) - F3_BASE(pf, EX, ix,iy,iz));
    } foreach_3d_end;
  }
}

static void
psc_push_fields_c_push_b_e(struct psc_push_fields *push, mfields_base_t *flds)
{
  f_real cnx = .5 * psc.dt / psc.dx[0];
  f_real cny = .5 * psc.dt / psc.dx[1];
  f_real cnz = .5 * psc.dt / psc.dx[2];

  if (psc.domain.gdims[0] == 1) {
    cnx = 0.;
  }
  if (psc.domain.gdims[1] == 1) {
    cny = 0.;
  }
  if (psc.domain.gdims[2] == 1) {
    cnz = 0.;
  }

  // E-field propagation E^(n+0.5), B^(n+1.0), j^(n+1.0) 
  //                  -> E^(n+1.0), B^(n+1.0), j^(n+1.0)

  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    foreach_3d(p, ix, iy, iz, 1, 1) {
      F3_BASE(pf, EX, ix,iy,iz) +=
	cny * (F3_BASE(pf, HZ, ix,iy,iz) - F3_BASE(pf, HZ, ix,iy-1,iz)) -
	cnz * (F3_BASE(pf, HY, ix,iy,iz) - F3_BASE(pf, HY, ix,iy,iz-1)) -
	.5 * psc.dt * F3_BASE(pf, JXI, ix,iy,iz);
      
      F3_BASE(pf, EY, ix,iy,iz) +=
	cnz * (F3_BASE(pf, HX, ix,iy,iz) - F3_BASE(pf, HX, ix,iy,iz-1)) -
	cnx * (F3_BASE(pf, HZ, ix,iy,iz) - F3_BASE(pf, HZ, ix-1,iy,iz)) -
	.5 * psc.dt * F3_BASE(pf, JYI, ix,iy,iz);
      
      F3_BASE(pf, EZ, ix,iy,iz) +=
	cnx * (F3_BASE(pf, HY, ix,iy,iz) - F3_BASE(pf, HY, ix-1,iy,iz)) -
	cny * (F3_BASE(pf, HX, ix,iy,iz) - F3_BASE(pf, HX, ix,iy-1,iz)) -
	.5 * psc.dt * F3_BASE(pf, JZI, ix,iy,iz);
    } foreach_3d_end;
  }
}

// ======================================================================
// psc_push_fields: subclass "c"

struct psc_push_fields_ops psc_push_fields_c_ops = {
  .name                  = "c",
  .push_a_e              = psc_push_fields_c_push_a_e,
  .push_a_h              = psc_push_fields_c_push_a_h,
  .push_b_h              = psc_push_fields_c_push_b_h,
  .push_b_e              = psc_push_fields_c_push_b_e,
};
