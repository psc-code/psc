
#include "psc.h"
#include <mrc_profile.h>

static void
c_push_field_a_nopml()
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_field_a", 1., 0, 0);
  }
  prof_start(pr);

  fields_base_t *pf = &psc.pf;

  assert(psc.domain.bnd_fld_lo[0] == BND_FLD_PERIODIC);
  assert(psc.domain.bnd_fld_lo[1] == BND_FLD_PERIODIC);
  assert(psc.domain.bnd_fld_lo[2] == BND_FLD_PERIODIC);
  assert(psc.domain.bnd_fld_hi[0] == BND_FLD_PERIODIC);
  assert(psc.domain.bnd_fld_hi[1] == BND_FLD_PERIODIC);
  assert(psc.domain.bnd_fld_hi[2] == BND_FLD_PERIODIC);

  f_real lx = psc.dt / psc.dx[0];
  f_real ly = psc.dt / psc.dx[1];
  f_real lz = psc.dt / psc.dx[2];

  f_real cnx = .5 * lx;
  f_real cny = .5 * ly;
  f_real cnz = .5 * lz;

  if (psc.domain.ihi[0] - psc.domain.ilo[0] == 1) {
    cnx = 0.;
  }
  if (psc.domain.ihi[1] - psc.domain.ilo[1] == 1) {
    cny = 0.;
  }
  if (psc.domain.ihi[2] - psc.domain.ilo[2] == 1) {
    cnz = 0.;
  }

  // E-field propagation E^(n)    , H^(n), j^(n) 
  //                  -> E^(n+0.5), H^(n), j^(n)

  foreach_3d(ix, iy, iz, 1, 1) {
    XF3_BASE(pf, EX, ix,iy,iz) +=
      cny * (XF3_BASE(pf, HZ, ix,iy,iz) - XF3_BASE(pf, HZ, ix,iy-1,iz)) -
      cnz * (XF3_BASE(pf, HY, ix,iy,iz) - XF3_BASE(pf, HY, ix,iy,iz-1)) -
      .5 * psc.dt * XF3_BASE(pf, JXI, ix,iy,iz);
    
    XF3_BASE(pf, EY, ix,iy,iz) +=
      cnz * (XF3_BASE(pf, HX, ix,iy,iz) - XF3_BASE(pf, HX, ix,iy,iz-1)) -
      cnx * (XF3_BASE(pf, HZ, ix,iy,iz) - XF3_BASE(pf, HZ, ix-1,iy,iz)) -
      .5 * psc.dt * XF3_BASE(pf, JYI, ix,iy,iz);
    
    XF3_BASE(pf, EZ, ix,iy,iz) +=
      cnx * (XF3_BASE(pf, HY, ix,iy,iz) - XF3_BASE(pf, HY, ix-1,iy,iz)) -
      cny * (XF3_BASE(pf, HX, ix,iy,iz) - XF3_BASE(pf, HX, ix,iy-1,iz)) -
      .5 * psc.dt * XF3_BASE(pf, JZI, ix,iy,iz);
  } foreach_3d_end;

  psc_fill_ghosts(&psc.pf, EX, EX + 3);

  // B-field propagation E^(n+0.5), H^(n    ), j^(n), m^(n+0.5)
  //                  -> E^(n+0.5), H^(n+0.5), j^(n), m^(n+0.5)

  foreach_3d(ix, iy, iz, 1, 1) {
    XF3_BASE(pf, HX, ix,iy,iz) -=
      cny * (XF3_BASE(pf, EZ, ix,iy+1,iz) - XF3_BASE(pf, EZ, ix,iy,iz)) -
      cnz * (XF3_BASE(pf, EY, ix,iy,iz+1) - XF3_BASE(pf, EY, ix,iy,iz));
    
    XF3_BASE(pf, HY, ix,iy,iz) -=
      cnz * (XF3_BASE(pf, EX, ix,iy,iz+1) - XF3_BASE(pf, EX, ix,iy,iz)) -
      cnx * (XF3_BASE(pf, EZ, ix+1,iy,iz) - XF3_BASE(pf, EZ, ix,iy,iz));
    
    XF3_BASE(pf, HZ, ix,iy,iz) -=
      cnx * (XF3_BASE(pf, EY, ix+1,iy,iz) - XF3_BASE(pf, EY, ix,iy,iz)) -
      cny * (XF3_BASE(pf, EX, ix,iy+1,iz) - XF3_BASE(pf, EX, ix,iy,iz));
  } foreach_3d_end;

  psc_fill_ghosts(&psc.pf, HX, HX + 3);

  prof_stop(pr);
}

static void
c_push_field_b_nopml()
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_field_b", 1., 0, 0);
  }
  prof_start(pr);

  fields_base_t *pf = &psc.pf;

  f_real lx = psc.dt / psc.dx[0];
  f_real ly = psc.dt / psc.dx[1];
  f_real lz = psc.dt / psc.dx[2];

  f_real cnx = .5 * lx;
  f_real cny = .5 * ly;
  f_real cnz = .5 * lz;

  if (psc.domain.ihi[0] - psc.domain.ilo[0] == 1) {
    cnx = 0.;
  }
  if (psc.domain.ihi[1] - psc.domain.ilo[1] == 1) {
    cny = 0.;
  }
  if (psc.domain.ihi[2] - psc.domain.ilo[2] == 1) {
    cnz = 0.;
  }

  // B-field propagation E^(n+0.5), B^(n+0.5), j^(n+1.0), m^(n+0.5)
  //                  -> E^(n+0.5), B^(n+1.0), j^(n+1.0), m^(n+0.5)

  foreach_3d(ix, iy, iz, 1, 1) {
    XF3_BASE(pf, HX, ix,iy,iz) -=
      cny * (XF3_BASE(pf, EZ, ix,iy+1,iz) - XF3_BASE(pf, EZ, ix,iy,iz)) -
      cnz * (XF3_BASE(pf, EY, ix,iy,iz+1) - XF3_BASE(pf, EY, ix,iy,iz));
    
    XF3_BASE(pf, HY, ix,iy,iz) -=
      cnz * (XF3_BASE(pf, EX, ix,iy,iz+1) - XF3_BASE(pf, EX, ix,iy,iz)) -
      cnx * (XF3_BASE(pf, EZ, ix+1,iy,iz) - XF3_BASE(pf, EZ, ix,iy,iz));
    
    XF3_BASE(pf, HZ, ix,iy,iz) -=
      cnx * (XF3_BASE(pf, EY, ix+1,iy,iz) - XF3_BASE(pf, EY, ix,iy,iz)) -
      cny * (XF3_BASE(pf, EX, ix,iy+1,iz) - XF3_BASE(pf, EX, ix,iy,iz));
  } foreach_3d_end;

  psc_fill_ghosts(&psc.pf, HX, HX + 3);

  // E-field propagation E^(n+0.5), B^(n+1.0), j^(n+1.0) 
  //                  -> E^(n+1.0), B^(n+1.0), j^(n+1.0)

  foreach_3d(ix, iy, iz, 1, 1) {
    XF3_BASE(pf, EX, ix,iy,iz) +=
      cny * (XF3_BASE(pf, HZ, ix,iy,iz) - XF3_BASE(pf, HZ, ix,iy-1,iz)) -
      cnz * (XF3_BASE(pf, HY, ix,iy,iz) - XF3_BASE(pf, HY, ix,iy,iz-1)) -
      .5 * psc.dt * XF3_BASE(pf, JXI, ix,iy,iz);
    
    XF3_BASE(pf, EY, ix,iy,iz) +=
      cnz * (XF3_BASE(pf, HX, ix,iy,iz) - XF3_BASE(pf, HX, ix,iy,iz-1)) -
      cnx * (XF3_BASE(pf, HZ, ix,iy,iz) - XF3_BASE(pf, HZ, ix-1,iy,iz)) -
      .5 * psc.dt * XF3_BASE(pf, JYI, ix,iy,iz);
    
    XF3_BASE(pf, EZ, ix,iy,iz) +=
      cnx * (XF3_BASE(pf, HY, ix,iy,iz) - XF3_BASE(pf, HY, ix-1,iy,iz)) -
      cny * (XF3_BASE(pf, HX, ix,iy,iz) - XF3_BASE(pf, HX, ix,iy-1,iz)) -
      .5 * psc.dt * XF3_BASE(pf, JZI, ix,iy,iz);
  } foreach_3d_end;

  psc_fill_ghosts(&psc.pf, EX, EX + 3);

  prof_stop(pr);
}

static void
c_push_field_a(void)
{
  if (psc.domain.use_pml) {
    assert(0);
  } else {
    c_push_field_a_nopml();
  }
}

static void
c_push_field_b(void)
{
  if (psc.domain.use_pml) {
    assert(0);
  } else {
    c_push_field_b_nopml();
  }
}

struct psc_push_field_ops psc_push_field_ops_c = {
  .name         = "c",
  .push_field_a = c_push_field_a,
  .push_field_b = c_push_field_b,
};

