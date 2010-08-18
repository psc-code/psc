
#include "psc.h"
#include "util/profile.h"

static void
c_push_field_a_nopml()
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_field_a", 1., 0, 0);
  }
  prof_start(pr);

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

  for (int iz = psc.ilo[2] - 1; iz < psc.ihi[2] + 1; iz++) {
    for (int iy = psc.ilo[1] - 1; iy < psc.ihi[1] + 1; iy++) {
      for (int ix = psc.ilo[0] - 1; ix < psc.ihi[0] + 1; ix++) {
	F3_BASE(EX, ix,iy,iz) +=
	  cny * (F3_BASE(HZ, ix,iy,iz) - F3_BASE(HZ, ix,iy-1,iz)) -
	  cnz * (F3_BASE(HY, ix,iy,iz) - F3_BASE(HY, ix,iy,iz-1)) -
	  .5 * psc.dt * F3_BASE(JXI, ix,iy,iz);

	F3_BASE(EY, ix,iy,iz) +=
	  cnz * (F3_BASE(HX, ix,iy,iz) - F3_BASE(HX, ix,iy,iz-1)) -
	  cnx * (F3_BASE(HZ, ix,iy,iz) - F3_BASE(HZ, ix-1,iy,iz)) -
	  .5 * psc.dt * F3_BASE(JYI, ix,iy,iz);

	F3_BASE(EZ, ix,iy,iz) +=
	  cnx * (F3_BASE(HY, ix,iy,iz) - F3_BASE(HY, ix-1,iy,iz)) -
	  cny * (F3_BASE(HX, ix,iy,iz) - F3_BASE(HX, ix,iy-1,iz)) -
	  .5 * psc.dt * F3_BASE(JZI, ix,iy,iz);
      }
    }
  }

  psc_fill_ghosts(EX, EX + 3);

  // B-field propagation E^(n+0.5), H^(n    ), j^(n), m^(n+0.5)
  //                  -> E^(n+0.5), H^(n+0.5), j^(n), m^(n+0.5)

  for (int iz = psc.ilo[2] - 1; iz < psc.ihi[2] + 1; iz++) {
    for (int iy = psc.ilo[1] - 1; iy < psc.ihi[1] + 1; iy++) {
      for (int ix = psc.ilo[0] - 1; ix < psc.ihi[0] + 1; ix++) {
	F3_BASE(HX, ix,iy,iz) -=
	  cny * (F3_BASE(EZ, ix,iy+1,iz) - F3_BASE(EZ, ix,iy,iz)) -
	  cnz * (F3_BASE(EY, ix,iy,iz+1) - F3_BASE(EY, ix,iy,iz));

	F3_BASE(HY, ix,iy,iz) -=
	  cnz * (F3_BASE(EX, ix,iy,iz+1) - F3_BASE(EX, ix,iy,iz)) -
	  cnx * (F3_BASE(EZ, ix+1,iy,iz) - F3_BASE(EZ, ix,iy,iz));

	F3_BASE(HZ, ix,iy,iz) -=
	  cnx * (F3_BASE(EY, ix+1,iy,iz) - F3_BASE(EY, ix,iy,iz)) -
	  cny * (F3_BASE(EX, ix,iy+1,iz) - F3_BASE(EX, ix,iy,iz));
      }
    }
  }

  psc_fill_ghosts(HX, HX + 3);

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

  for (int iz = psc.ilo[2] - 1; iz < psc.ihi[2] + 1; iz++) {
    for (int iy = psc.ilo[1] - 1; iy < psc.ihi[1] + 1; iy++) {
      for (int ix = psc.ilo[0] - 1; ix < psc.ihi[0] + 1; ix++) {
	F3_BASE(HX, ix,iy,iz) -=
	  cny * (F3_BASE(EZ, ix,iy+1,iz) - F3_BASE(EZ, ix,iy,iz)) -
	  cnz * (F3_BASE(EY, ix,iy,iz+1) - F3_BASE(EY, ix,iy,iz));

	F3_BASE(HY, ix,iy,iz) -=
	  cnz * (F3_BASE(EX, ix,iy,iz+1) - F3_BASE(EX, ix,iy,iz)) -
	  cnx * (F3_BASE(EZ, ix+1,iy,iz) - F3_BASE(EZ, ix,iy,iz));

	F3_BASE(HZ, ix,iy,iz) -=
	  cnx * (F3_BASE(EY, ix+1,iy,iz) - F3_BASE(EY, ix,iy,iz)) -
	  cny * (F3_BASE(EX, ix,iy+1,iz) - F3_BASE(EX, ix,iy,iz));
      }
    }
  }

  psc_fill_ghosts(HX, HX + 3);

  // E-field propagation E^(n+0.5), B^(n+1.0), j^(n+1.0) 
  //                  -> E^(n+1.0), B^(n+1.0), j^(n+1.0)

  for (int iz = psc.ilo[2] - 1; iz < psc.ihi[2] + 1; iz++) {
    for (int iy = psc.ilo[1] - 1; iy < psc.ihi[1] + 1; iy++) {
      for (int ix = psc.ilo[0] - 1; ix < psc.ihi[0] + 1; ix++) {
	F3_BASE(EX, ix,iy,iz) +=
	  cny * (F3_BASE(HZ, ix,iy,iz) - F3_BASE(HZ, ix,iy-1,iz)) -
	  cnz * (F3_BASE(HY, ix,iy,iz) - F3_BASE(HY, ix,iy,iz-1)) -
	  .5 * psc.dt * F3_BASE(JXI, ix,iy,iz);

	F3_BASE(EY, ix,iy,iz) +=
	  cnz * (F3_BASE(HX, ix,iy,iz) - F3_BASE(HX, ix,iy,iz-1)) -
	  cnx * (F3_BASE(HZ, ix,iy,iz) - F3_BASE(HZ, ix-1,iy,iz)) -
	  .5 * psc.dt * F3_BASE(JYI, ix,iy,iz);

	F3_BASE(EZ, ix,iy,iz) +=
	  cnx * (F3_BASE(HY, ix,iy,iz) - F3_BASE(HY, ix-1,iy,iz)) -
	  cny * (F3_BASE(HX, ix,iy,iz) - F3_BASE(HX, ix,iy-1,iz)) -
	  .5 * psc.dt * F3_BASE(JZI, ix,iy,iz);
      }
    }
  }

  psc_fill_ghosts(EX, EX + 3);

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

