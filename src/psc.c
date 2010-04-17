
#include "psc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

struct psc psc;

// ----------------------------------------------------------------------
// assert_equal
//
// make sure that the two values are almost equal.

void
assert_equal(double x, double y)
{
  double max = fmax(fabs(x), fabs(y)) + 1e-10;
  double eps = fabs((x - y) / max);
  if (eps > 1e-5) {
    fprintf(stderr, "assert_equal: fail x = %g y = %g rel err = %g\n",
	    x, y, eps);
    abort();
  }
}

// ----------------------------------------------------------------------

static struct psc_ops *psc_ops_list[] = {
  &psc_ops_generic_c,
  &psc_ops_fortran,
  &psc_ops_cuda,
  NULL,
};

static struct psc_ops *
psc_find_ops(const char *ops_name)
{
  for (int i = 0; psc_ops_list[i]; i++) {
    if (strcmp(psc_ops_list[i]->name, ops_name) == 0)
      return psc_ops_list[i];
  }
  fprintf(stderr, "ERROR: psc_ops '%s' not available.\n", ops_name);
  abort();
}

void
psc_create(const char *ops_name)
{
  memset(&psc, 0, sizeof(psc));
  psc.ops = psc_find_ops(ops_name);
  if (psc.ops->create) {
    psc.ops->create();
  }
}

// ----------------------------------------------------------------------

void
psc_alloc(int ilo[3], int ihi[3], int ibn[3], int n_part)
{
  for (int d = 0; d < 3; d++) {
    psc.ilo[d] = ilo[d];
    psc.ilg[d] = ilo[d] - ibn[d];
    psc.ihi[d] = ihi[d];
    psc.ihg[d] = ihi[d] + ibn[d];
    psc.img[d] = ihi[d] - ilo[d] + 2 * ibn[d];
    psc.ibn[d] = ibn[d];
  }
  psc.fld_size = psc.img[0] * psc.img[1] * psc.img[2];
  for (int m = 0; m < NR_FIELDS; m++) {
    psc.f_fields[m] = calloc(psc.fld_size, sizeof(f_real));
  }

  psc.n_part = n_part;
  psc.f_part = calloc(n_part, sizeof(*psc.f_part));

  psc.allocated = true;
}

void
psc_destroy()
{
  if (psc.ops->destroy) {
    psc.ops->destroy();
  }

  if (psc.allocated) {
    for (int m = 0; m < NR_FIELDS; m++) {
      free(psc.f_fields[m]);
    }

    free(psc.f_part);
  }
}

void
psc_setup_parameters()
{
  psc.prm.cori = 2.;
  psc.prm.eta = 3.;
  psc.prm.alpha = 5.;
  psc.dt = 1.;
  psc.dx[0] = 1.;
  psc.dx[1] = 1.;
  psc.dx[2] = 1.;
}

void
psc_setup_fields_zero()
{
  for (int m = 0; m < NR_FIELDS; m++) {
    memset(psc.f_fields[m], 0, psc.fld_size * sizeof(f_real));
  }
}

void
psc_setup_fields_1()
{
  psc_setup_fields_zero();
  for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
    for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
      for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	FF3(EX, jx,jy,jz) = .1 * sin(.5 * jx) + .2 * sin(.4 * jy) + .3 * sin(.3 * jz);
	FF3(EY, jx,jy,jz) = .2 * sin(.4 * jx) + .3 * sin(.3 * jy) + .4 * sin(.2 * jz);
	FF3(EZ, jx,jy,jz) = .3 * sin(.3 * jx) + .4 * sin(.2 * jy) + .5 * sin(.1 * jz);
	FF3(BX, jx,jy,jz) = .1 * cos(.5 * jx) + .2 * cos(.4 * jy) + .3 * cos(.3 * jz);
	FF3(BY, jx,jy,jz) = .2 * cos(.4 * jx) + .3 * cos(.3 * jy) + .4 * cos(.2 * jz);
	FF3(BZ, jx,jy,jz) = .3 * cos(.3 * jx) + .4 * cos(.2 * jy) + .5 * cos(.1 * jz);
      }
    }
  }
}

void
psc_setup_particles_1()
{
  int n = 0;
  int n_per_cell = psc.n_part / ((psc.ihi[1]-psc.ilo[1])*(psc.ihi[2]-psc.ilo[2]));
  for (int iz = psc.ilo[2]; iz < psc.ihi[2]; iz++) {
    for (int iy = psc.ilo[1]; iy < psc.ihi[1]; iy++) {
      for (int cnt = 0; cnt < n_per_cell; cnt++) {
	psc.f_part[n].xi = .2;
	psc.f_part[n].yi = iy;
	psc.f_part[n].zi = iz;
	psc.f_part[n].pxi = 0.;
	psc.f_part[n].pyi = .02;
	psc.f_part[n].pzi = .01;
	psc.f_part[n].qni = -1.;
	psc.f_part[n].mni = 1.;
	psc.f_part[n].lni = 0.;
	psc.f_part[n].wni = 1.;
	n++;
      }
    }
  }
}

void
psc_dump_particles(const char *fname)
{
  printf("psc_dump_particles %s\n", fname);

  FILE *file = fopen(fname, "w");
  fprintf(file, "i\txi\tyi\tzi\tpxi\tpyi\tpzi\tqni\tmni\twni\n");
  for (int i = 0; i < psc.n_part; i++) {
    struct f_particle *p = &psc.f_part[i];
    fprintf(file, "%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	    i, p->xi, p->yi, p->zi,
	    p->pxi, p->pyi, p->pzi, p->qni, p->mni, p->wni);
  }
  fclose(file);
}

// ----------------------------------------------------------------------
// psc_push_part_yz

void
psc_push_part_yz()
{
  assert(psc.ops->push_part_yz);
  psc.ops->particles_from_fortran();
  psc.ops->push_part_yz();
  psc.ops->particles_to_fortran();
}

// ----------------------------------------------------------------------
// psc_push_part_z

void
psc_push_part_z()
{
  assert(psc.ops->push_part_z);
  psc.ops->particles_from_fortran();
  psc.ops->push_part_z();
  psc.ops->particles_to_fortran();
}

// ----------------------------------------------------------------------
// psc_push_part_yz_a

void
psc_push_part_yz_a()
{
  assert(psc.ops->push_part_yz_a);
  psc.ops->particles_from_fortran();
  psc.ops->push_part_yz_a();
  psc.ops->particles_to_fortran();
}

// ----------------------------------------------------------------------
// psc_push_part_yz_b

void
psc_push_part_yz_b()
{
  assert(psc.ops->push_part_yz_b);
  psc.ops->particles_from_fortran();
  psc.ops->fields_from_fortran();
  psc.ops->push_part_yz_b();
  psc.ops->fields_to_fortran();
  psc.ops->particles_to_fortran();
}


static struct f_particle *particle_ref;

// ----------------------------------------------------------------------
// psc_save_particles_ref
//
// save current particle data as reference solution

void
psc_save_particles_ref()
{
  if (!particle_ref) {
    particle_ref = calloc(psc.n_part, sizeof(*particle_ref));
  }
  for (int i = 0; i < psc.n_part; i++) {
    particle_ref[i] = psc.f_part[i];
  }
}

// ----------------------------------------------------------------------
// psc_check_particles_ref
//
// check current particle data agains previously saved reference solution

void
psc_check_particles_ref()
{
  assert(particle_ref);
  for (int i = 0; i < psc.n_part; i++) {
    assert_equal(psc.f_part[i].xi , particle_ref[i].xi);
    assert_equal(psc.f_part[i].yi , particle_ref[i].yi);
    assert_equal(psc.f_part[i].zi , particle_ref[i].zi);
    assert_equal(psc.f_part[i].pxi, particle_ref[i].pxi);
    assert_equal(psc.f_part[i].pyi, particle_ref[i].pyi);
    assert_equal(psc.f_part[i].pzi, particle_ref[i].pzi);
  }
}

// ----------------------------------------------------------------------
// psc_create_test_1
//
// set up test case 1

void
psc_create_test_1(const char *ops_name)
{
  int ilo[3] = { 0,  0,  0 };
  int ihi[3] = { 1, 16, 16 };
  int ibn[3] = { 2,  2,  2 }; // FIXME?

  int n_part = 1e4 * (ihi[2]-ilo[2]) * (ihi[1] - ilo[1]);

  psc_create(ops_name);
  psc_alloc(ilo, ihi, ibn, n_part);
  psc_setup_parameters();
  psc_setup_fields_1();
  psc_setup_particles_1();
}

