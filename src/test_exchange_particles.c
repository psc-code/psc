
#include "psc.h"
#include "util/profile.h"
#include "util/params.h"

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

void
setup_particles(void)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // FIXME, realloc
  // only particles on proc 1, but some are out of bounds.
  // particles are right on nodes of the grid, but as far as
  // particle -> cell mapping goes, that's right in the center of the
  // particle ordering cell, for reasons that should be laid out more clearly.
  // for the old ordering, nodes aren't good because they're indeterminate
  // (could go either way), so let's shift them a bit so we get a unique answer
  // we can check.
  if (rank == 0) {
    int i = 0;
    for (int iz = psc.ilo[2]-1; iz < psc.ihi[2]+1; iz++) {
      for (int iy = psc.ilo[1]-1; iy < psc.ihi[1]+1; iy++) {
	for (int ix = psc.ilo[0]-1; ix < psc.ihi[0]+1; ix++) {
	  struct f_particle *p;
	  p = &psc.f_part[i++];
	  memset(p, 0, sizeof(*p));
	  p->xi = (ix + .01) * psc.dx[0];
	  p->yi = (iy + .01) * psc.dx[1];
	  p->zi = (iz + .01) * psc.dx[2];

	  p = &psc.f_part[i++];
	  memset(p, 0, sizeof(*p));
	  p->xi = (ix - .01) * psc.dx[0];
	  p->yi = (iy - .01) * psc.dx[1];
	  p->zi = (iz - .01) * psc.dx[2];
	}
      }
    }
    psc.n_part = i;
  } else {
    psc.n_part = 0;
  }
}

static void
check_particles_old_xz(void)
{
  f_real xb[3], xe[3];
  
  // These boundaries are, I believe, what the Fortran code guarantees.
  // The rest of the code works under these conditions, but I don't like
  // not having one strict owner for each location of a particle in the domain
  // might be.
  // OTOH, there might be a point to have valid regions overlap so that particles
  // gyrating near a boundary don't keep getting bounced between processors --
  // but that's the opposite of what's happening here.
  for (int d = 0; d < 3; d++) {
    xb[d] = (psc.ilo[d]-2) * psc.dx[d];
    xe[d] = (psc.ihi[d]+1) * psc.dx[d];
  }

  int fail_cnt = 0;
  for (int i = 0; i < psc.n_part; i++) {
    struct f_particle *p = &psc.f_part[i];
    if (p->xi < xb[0] || p->xi > xe[0] ||
	p->zi < xb[2] || p->zi > xe[2]) {
      if (fail_cnt++ < 10) {
	printf("FAIL old: xi %g [%g:%g]\n", p->xi, xb[0], xe[0]);
	printf("          zi %g [%g:%g]\n", p->zi, xb[2], xe[2]);
      }
    }
  }
  assert(fail_cnt == 0);
}

static void
dump_particles(const char *s)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  char fname[200];
  sprintf(fname, "%s-p%d.asc", s, rank);
  psc_dump_particles(fname);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  params_init(argc, argv);

  // test psc_exchange_particles()

  struct psc_mod_config conf_fortran = {
    .mod_bnd = "fortran",
  };

  psc_create_test_xz(&conf_fortran);
  setup_particles();
  dump_particles("part-0");
  psc_exchange_particles();
  dump_particles("part-1");
  check_particles_old_xz();
  psc_destroy();

  prof_print();

  MPI_Finalize();
}
