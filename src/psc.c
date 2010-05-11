
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
  double max = fmax(fabs(x), fabs(y)) + 1e-7;
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
  &psc_ops_sse2,
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

// Fills the fields with sequentially increasing numbers
// I use it to check macros and indexing

void
psc_setup_fields_seq_yz()
{
  for(int m = EX; m < NR_FIELDS; m++){
    for(int i = 0; i < psc.img[0]; i++){
      for(int j = 0; j < psc.img[1]; j++){
	for(int k = 0; k < psc.img[2]; k++){
	  psc.f_fields[m][i + j*(psc.img[0]) + k*(psc.img[0]*psc.img[1])] = j*k;
	}
      }
    }
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
  for (int n = 0; n < psc.n_part; n++) {
    psc.f_part[n].xi = .5;
    psc.f_part[n].yi = .5;
    psc.f_part[n].zi = .5 + n / (real) psc.n_part;
    psc.f_part[n].pxi = 0.;
    psc.f_part[n].pyi = .02;
    psc.f_part[n].pzi = .01;
    psc.f_part[n].qni = -1.;
    psc.f_part[n].mni = 1.;
    psc.f_part[n].lni = 0.;
    psc.f_part[n].wni = 1.;
  }
}

void
psc_setup_particles_random_yz()
{
  srand48(42.0);
  for (int n = 0; n < psc.n_part; n++) {
    psc.f_part[n].xi = .5;
    psc.f_part[n].yi = psc.ilo[1] + (psc.ihi[1] - psc.ilo[1])*drand48();
    psc.f_part[n].zi = psc.ilo[2] + (psc.ihi[2] - psc.ilo[2])*drand48();;
    psc.f_part[n].pxi = 0.;
    psc.f_part[n].pyi = .02;
    psc.f_part[n].pzi = .01;
    psc.f_part[n].qni = -1.;
    psc.f_part[n].mni = 1.;
    psc.f_part[n].lni = 0.;
    psc.f_part[n].wni = 1.;
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
  psc.ops->fields_from_fortran();
  psc.ops->push_part_yz();
  psc.ops->particles_to_fortran();
  psc.ops->fields_to_fortran();
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
  psc.ops->particles_to_fortran();
  psc.ops->fields_to_fortran();
}


static struct f_particle *particle_ref;
static f_real *field_ref[NR_FIELDS];
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
// psc_save_fields_ref
//
// save current field data as reference solution

void
psc_save_fields_ref()
{
  if (!field_ref[EX]) { //FIXME this is bad mojo
    for (int m = 0; m < NR_FIELDS; m++) {
      field_ref[m] = calloc(psc.fld_size, sizeof(f_real));
    }
  }
  for (int m = 0; m < NR_FIELDS; m++) {
    for (int n = 0; n < psc.fld_size; n++){
      field_ref[m][n] = psc.f_fields[m][n];
    }
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
// psc_check_fields_ref
//
// check current field data agains previously saved reference solution

void
psc_check_fields_ref()
{
  assert(field_ref[NE]); //FIXME: this is bad
  for (int m = 0; m < NR_FIELDS; m++){
    for (int n = 0; n < psc.fld_size; n++){
      assert_equal(psc.f_fields[m][n], field_ref[m][n]);
      //  fprintf(stderr, "%d\n", n);
    }
  }
}

// ----------------------------------------------------------------------
// psc_check_currents_ref
//
// check current current density data agains previously saved reference solution

void
psc_check_currents_ref()
{
  assert(field_ref[JXI]); //FIXME: this is bad
  for (int m = JXI; m <= JZI; m++){
    for (int n = 0; n < psc.fld_size; n++){
      assert_equal(psc.f_fields[m][n], field_ref[m][n]);
      //  fprintf(stderr, "%d\n", n);
    }
    HERE;
  }
}

// ----------------------------------------------------------------------
// psc_create_test_1
//
// set up test case 1

void
psc_create_test_1(const char *ops_name)
{
  int ilo[3] = { 0, 0, 0 };
  int ihi[3] = { 1, 20, 20 };
  int ibn[3] = { 2, 2, 2 }; // FIXME?

  int n_part = 1e6;

  psc_create(ops_name);
  psc_alloc(ilo, ihi, ibn, n_part);
  psc_setup_parameters();
  psc_setup_fields_zero();
  psc_setup_particles_1();
}

// ----------------------------------------------------------------------
// psc_create_test_2
//
// set up test case 2 with some non-zero fields

void
psc_create_test_2(const char *ops_name)
{
  int ilo[3] = { 0, 0, 0 };
  int ihi[3] = { 1, 8, 8 };
  int ibn[3] = { 2, 2, 2 }; // FIXME?

  int n_part = 1e6;

  psc_create(ops_name);
  psc_alloc(ilo, ihi, ibn, n_part);
  psc_setup_parameters();
  psc_setup_fields_1();
  psc_setup_particles_1();
}

// ----------------------------------------------------------------------
// psc_create_test_3
//
// A test case with a slightly more realistic setup

void
psc_create_test_3(const char *ops_name)
{
  int ilo[3] = { 0, 0, 0 };
  int ihi[3] = { 1, 200, 200 };
  int ibn[3] = { 2, 2, 2 }; // FIXME?

  int n_part = 4e6;

  psc_create(ops_name);
  psc_alloc(ilo, ihi, ibn, n_part);
  psc_setup_parameters();
  psc_setup_fields_1();
  psc_setup_particles_random_yz();
}
