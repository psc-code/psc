
#include "psc.h"
#include "util/profile.h"

static void
c_push_field_a()
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_push_field_a", 1., 0, 0);
  }
  prof_start(pr);

  f_real lx = psc.dt / psc.dx[0];
  f_real ly = psc.dt / psc.dx[1];
  f_real lz = psc.dt / psc.dx[2];

  f_real cnx = .5 * lx;
  f_real cny = .5 * ly;
  f_real cnz = .5 * lz;

  // E-field propagation E^(n)    , B^(n), j^(n) 
  //                  -> E^(n+0.5), B^(n), j^(n)

  for (int iz = psc.ilo[2] - 1; iz < psc.ihi[2] + 1; iz++) {
    for (int iy = psc.ilo[1] - 1; iy < psc.ihi[1] + 1; iy++) {
      for (int ix = psc.ilo[0] - 1; ix < psc.ihi[0] + 1; ix++) {
	F3_BASE(EX, ix,iy,iz) +=
	  cny * (F3_BASE(BZ, ix,iy,iz) - F3_BASE(BZ, ix,iy-1,iz)) -
	  cnz * (F3_BASE(BY, ix,iy,iz) - F3_BASE(BY, ix,iy,iz-1)) -
	  .5 * psc.dt * F3_BASE(JXI, ix,iy,iz);

	F3_BASE(EY, ix,iy,iz) +=
	  cnz * (F3_BASE(BX, ix,iy,iz) - F3_BASE(BX, ix,iy,iz-1)) -
	  cnx * (F3_BASE(BZ, ix,iy,iz) - F3_BASE(BZ, ix-1,iy,iz)) -
	  .5 * psc.dt * F3_BASE(JYI, ix,iy,iz);

	F3_BASE(EZ, ix,iy,iz) +=
	  cnx * (F3_BASE(BY, ix,iy,iz) - F3_BASE(BY, ix-1,iy,iz)) -
	  cny * (F3_BASE(BX, ix,iy,iz) - F3_BASE(BX, ix,iy-1,iz)) -
	  .5 * psc.dt * F3_BASE(JZI, ix,iy,iz);
      }
    }
  }

  psc_fill_ghosts(EX, EX + 3);

  // B-field propagation E^(n+0.5), B^(n    ), j^(n), m^(n+0.5)
  //                  -> E^(n+0.5), B^(n+0.5), j^(n), m^(n+0.5)

  for (int iz = psc.ilo[2] - 1; iz < psc.ihi[2] + 1; iz++) {
    for (int iy = psc.ilo[1] - 1; iy < psc.ihi[1] + 1; iy++) {
      for (int ix = psc.ilo[0] - 1; ix < psc.ihi[0] + 1; ix++) {
	F3_BASE(BX, ix,iy,iz) -=
	  cny * (F3_BASE(EZ, ix,iy+1,iz) - F3_BASE(EZ, ix,iy,iz)) -
	  cnz * (F3_BASE(EY, ix,iy,iz+1) - F3_BASE(EY, ix,iy,iz));

	F3_BASE(BY, ix,iy,iz) -=
	  cnz * (F3_BASE(EX, ix,iy,iz+1) - F3_BASE(EX, ix,iy,iz)) -
	  cnx * (F3_BASE(EZ, ix+1,iy,iz) - F3_BASE(EZ, ix,iy,iz));

	F3_BASE(BZ, ix,iy,iz) -=
	  cnx * (F3_BASE(EY, ix+1,iy,iz) - F3_BASE(EY, ix,iy,iz)) -
	  cny * (F3_BASE(EX, ix,iy+1,iz) - F3_BASE(EX, ix,iy,iz));
      }
    }
  }

  psc_fill_ghosts(BX, BX + 3);

  prof_stop(pr);
}

static void
c_push_field_b()
{
  static int pr;
  if (!pr) {
    pr = prof_register("c_push_field_b", 1., 0, 0);
  }
  prof_start(pr);

  f_real lx = psc.dt / psc.dx[0];
  f_real ly = psc.dt / psc.dx[1];
  f_real lz = psc.dt / psc.dx[2];

  f_real cnx = .5 * lx;
  f_real cny = .5 * ly;
  f_real cnz = .5 * lz;

  // B-field propagation E^(n+0.5), B^(n+0.5), j^(n+1.0), m^(n+0.5)
  //                  -> E^(n+0.5), B^(n+1.0), j^(n+1.0), m^(n+0.5)

  for (int iz = psc.ilo[2] - 1; iz < psc.ihi[2] + 1; iz++) {
    for (int iy = psc.ilo[1] - 1; iy < psc.ihi[1] + 1; iy++) {
      for (int ix = psc.ilo[0] - 1; ix < psc.ihi[0] + 1; ix++) {
	F3_BASE(BX, ix,iy,iz) -=
	  cny * (F3_BASE(EZ, ix,iy+1,iz) - F3_BASE(EZ, ix,iy,iz)) -
	  cnz * (F3_BASE(EY, ix,iy,iz+1) - F3_BASE(EY, ix,iy,iz));

	F3_BASE(BY, ix,iy,iz) -=
	  cnz * (F3_BASE(EX, ix,iy,iz+1) - F3_BASE(EX, ix,iy,iz)) -
	  cnx * (F3_BASE(EZ, ix+1,iy,iz) - F3_BASE(EZ, ix,iy,iz));

	F3_BASE(BZ, ix,iy,iz) -=
	  cnx * (F3_BASE(EY, ix+1,iy,iz) - F3_BASE(EY, ix,iy,iz)) -
	  cny * (F3_BASE(EX, ix,iy+1,iz) - F3_BASE(EX, ix,iy,iz));
      }
    }
  }

  psc_fill_ghosts(BX, BX + 3);

  // E-field propagation E^(n+0.5), B^(n+1.0), j^(n+1.0) 
  //                  -> E^(n+1.0), B^(n+1.0), j^(n+1.0)

  for (int iz = psc.ilo[2] - 1; iz < psc.ihi[2] + 1; iz++) {
    for (int iy = psc.ilo[1] - 1; iy < psc.ihi[1] + 1; iy++) {
      for (int ix = psc.ilo[0] - 1; ix < psc.ihi[0] + 1; ix++) {
	F3_BASE(EX, ix,iy,iz) +=
	  cny * (F3_BASE(BZ, ix,iy,iz) - F3_BASE(BZ, ix,iy-1,iz)) -
	  cnz * (F3_BASE(BY, ix,iy,iz) - F3_BASE(BY, ix,iy,iz-1)) -
	  .5 * psc.dt * F3_BASE(JXI, ix,iy,iz);

	F3_BASE(EY, ix,iy,iz) +=
	  cnz * (F3_BASE(BX, ix,iy,iz) - F3_BASE(BX, ix,iy,iz-1)) -
	  cnx * (F3_BASE(BZ, ix,iy,iz) - F3_BASE(BZ, ix-1,iy,iz)) -
	  .5 * psc.dt * F3_BASE(JYI, ix,iy,iz);

	F3_BASE(EZ, ix,iy,iz) +=
	  cnx * (F3_BASE(BY, ix,iy,iz) - F3_BASE(BY, ix-1,iy,iz)) -
	  cny * (F3_BASE(BX, ix,iy,iz) - F3_BASE(BX, ix,iy-1,iz)) -
	  .5 * psc.dt * F3_BASE(JZI, ix,iy,iz);
      }
    }
  }

  psc_fill_ghosts(EX, EX + 3);

  prof_stop(pr);
}

struct psc_push_field_ops psc_push_field_ops_c = {
  .name         = "c",
  .push_field_a = c_push_field_a,
  .push_field_b = c_push_field_b,
};

