
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
	FF3(EX, ix,iy,iz) +=
	  cny * (FF3(BZ, ix,iy,iz) - FF3(BZ, ix,iy-1,iz)) -
	  cnz * (FF3(BY, ix,iy,iz) - FF3(BY, ix,iy,iz-1)) -
	  .5 * psc.dt * FF3(JXI, ix,iy,iz);

	FF3(EY, ix,iy,iz) +=
	  cnz * (FF3(BX, ix,iy,iz) - FF3(BX, ix,iy,iz-1)) -
	  cnx * (FF3(BZ, ix,iy,iz) - FF3(BZ, ix-1,iy,iz)) -
	  .5 * psc.dt * FF3(JYI, ix,iy,iz);

	FF3(EZ, ix,iy,iz) +=
	  cnx * (FF3(BY, ix,iy,iz) - FF3(BY, ix-1,iy,iz)) -
	  cny * (FF3(BX, ix,iy,iz) - FF3(BX, ix,iy-1,iz)) -
	  .5 * psc.dt * FF3(JZI, ix,iy,iz);
      }
    }
  }

  psc_fill_ghosts(EX);
  psc_fill_ghosts(EY);
  psc_fill_ghosts(EZ);

  // B-field propagation E^(n+0.5), B^(n    ), j^(n), m^(n+0.5)
  //                  -> E^(n+0.5), B^(n+0.5), j^(n), m^(n+0.5)

  for (int iz = psc.ilo[2] - 1; iz < psc.ihi[2] + 1; iz++) {
    for (int iy = psc.ilo[1] - 1; iy < psc.ihi[1] + 1; iy++) {
      for (int ix = psc.ilo[0] - 1; ix < psc.ihi[0] + 1; ix++) {
	FF3(BX, ix,iy,iz) -=
	  cny * (FF3(EZ, ix,iy+1,iz) - FF3(EZ, ix,iy,iz)) -
	  cnz * (FF3(EY, ix,iy,iz+1) - FF3(EY, ix,iy,iz));

	FF3(BY, ix,iy,iz) -=
	  cnz * (FF3(EX, ix,iy,iz+1) - FF3(EX, ix,iy,iz)) -
	  cnx * (FF3(EZ, ix+1,iy,iz) - FF3(EZ, ix,iy,iz));

	FF3(BZ, ix,iy,iz) -=
	  cnx * (FF3(EY, ix+1,iy,iz) - FF3(EY, ix,iy,iz)) -
	  cny * (FF3(EX, ix,iy+1,iz) - FF3(EX, ix,iy,iz));
      }
    }
  }

  psc_fill_ghosts(BX);
  psc_fill_ghosts(BY);
  psc_fill_ghosts(BZ);

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
	FF3(BX, ix,iy,iz) -=
	  cny * (FF3(EZ, ix,iy+1,iz) - FF3(EZ, ix,iy,iz)) -
	  cnz * (FF3(EY, ix,iy,iz+1) - FF3(EY, ix,iy,iz));

	FF3(BY, ix,iy,iz) -=
	  cnz * (FF3(EX, ix,iy,iz+1) - FF3(EX, ix,iy,iz)) -
	  cnx * (FF3(EZ, ix+1,iy,iz) - FF3(EZ, ix,iy,iz));

	FF3(BZ, ix,iy,iz) -=
	  cnx * (FF3(EY, ix+1,iy,iz) - FF3(EY, ix,iy,iz)) -
	  cny * (FF3(EX, ix,iy+1,iz) - FF3(EX, ix,iy,iz));
      }
    }
  }

  psc_fill_ghosts(BX);
  psc_fill_ghosts(BY);
  psc_fill_ghosts(BZ);

  // E-field propagation E^(n+0.5), B^(n+1.0), j^(n+1.0) 
  //                  -> E^(n+1.0), B^(n+1.0), j^(n+1.0)

  for (int iz = psc.ilo[2] - 1; iz < psc.ihi[2] + 1; iz++) {
    for (int iy = psc.ilo[1] - 1; iy < psc.ihi[1] + 1; iy++) {
      for (int ix = psc.ilo[0] - 1; ix < psc.ihi[0] + 1; ix++) {
	FF3(EX, ix,iy,iz) +=
	  cny * (FF3(BZ, ix,iy,iz) - FF3(BZ, ix,iy-1,iz)) -
	  cnz * (FF3(BY, ix,iy,iz) - FF3(BY, ix,iy,iz-1)) -
	  .5 * psc.dt * FF3(JXI, ix,iy,iz);

	FF3(EY, ix,iy,iz) +=
	  cnz * (FF3(BX, ix,iy,iz) - FF3(BX, ix,iy,iz-1)) -
	  cnx * (FF3(BZ, ix,iy,iz) - FF3(BZ, ix-1,iy,iz)) -
	  .5 * psc.dt * FF3(JYI, ix,iy,iz);

	FF3(EZ, ix,iy,iz) +=
	  cnx * (FF3(BY, ix,iy,iz) - FF3(BY, ix-1,iy,iz)) -
	  cny * (FF3(BX, ix,iy,iz) - FF3(BX, ix,iy-1,iz)) -
	  .5 * psc.dt * FF3(JZI, ix,iy,iz);
      }
    }
  }

  psc_fill_ghosts(EX);
  psc_fill_ghosts(EY);
  psc_fill_ghosts(EZ);

  prof_stop(pr);
}

struct psc_push_field_ops psc_push_field_ops_c = {
  .name         = "c",
  .push_field_a = c_push_field_a,
  .push_field_b = c_push_field_b,
};

