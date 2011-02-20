
#include "psc_testing.h"

#include <math.h>
#include <limits.h>

// ----------------------------------------------------------------------
// assert_equal
//
// make sure that the two values are almost equal.

void
__assert_equal(double x, double y, const char *xs, const char *ys, double thres)
{
  double max = fmax(fabs(x), fabs(y));
  if (max < thres)
    return;

  double eps = fabs((x - y) / max);
  if (eps > thres) {
    fprintf(stderr, "assert_equal: fail %s = %g, %s = %g rel err = %g thres = %g\n",
	    xs, x, ys, y, eps, thres);
    abort();
  }
}

static particle_base_t *particle_ref;
static fields_base_t *field_ref;

// ----------------------------------------------------------------------
// psc_save_particles_ref
//
// save current particle data as reference solution

void
psc_save_particles_ref()
{
  if (!particle_ref) {
    particle_ref = calloc(psc.pp.n_part, sizeof(*particle_ref));
  }
  for (int i = 0; i < psc.pp.n_part; i++) {
    particle_ref[i] = *particles_base_get_one(&psc.pp, i);
  }
}

// ----------------------------------------------------------------------
// psc_save_fields_ref
//
// save current field data as reference solution

void
psc_save_fields_ref()
{
  if (!field_ref) {
    field_ref = calloc(1, sizeof(*field_ref));
    struct psc_patch *patch = &psc.patch[0];
    int ilg[3] = { -psc.ibn[0], -psc.ibn[1], -psc.ibn[2] };
    int ihg[3] = { patch->ldims[0] + psc.ibn[0],
		   patch->ldims[1] + psc.ibn[1],
		   patch->ldims[2] + psc.ibn[2] };
    fields_base_alloc(&psc.pf, ilg, ihg, NR_FIELDS);
  }
  int me = psc.domain.use_pml ? NR_FIELDS : HZ + 1;
  foreach_patch(patch) {
    for (int m = 0; m < me; m++) {
      foreach_3d_g(patch, ix, iy, iz) {
	F3_BASE(field_ref,m, ix,iy,iz) = F3_BASE(&psc.pf, m, ix,iy,iz);
      } foreach_3d_g_end;
    }
  }
} 

// ----------------------------------------------------------------------
// psc_check_particles_ref
//
// check current particle data agains previously saved reference solution

void
psc_check_particles_ref(double thres, const char *test_str)
{
  assert(particle_ref);
  particle_base_real_t xi = 0., yi = 0., zi = 0., pxi = 0., pyi = 0., pzi = 0.;
  for (int i = 0; i < psc.pp.n_part; i++) {
    particle_base_t *part = particles_base_get_one(&psc.pp, i);
    //    printf("i = %d\n", i);
    xi  = fmax(xi , fabs(part->xi  - particle_ref[i].xi));
    yi  = fmax(yi , fabs(part->yi  - particle_ref[i].yi));
    zi  = fmax(zi , fabs(part->zi  - particle_ref[i].zi));
    pxi = fmax(pxi, fabs(part->pxi - particle_ref[i].pxi));
    pyi = fmax(pyi, fabs(part->pyi - particle_ref[i].pyi));
    pzi = fmax(pzi, fabs(part->pzi - particle_ref[i].pzi));
    assert_equal(part->xi , particle_ref[i].xi, thres);
    assert_equal(part->yi , particle_ref[i].yi, thres);
    assert_equal(part->zi , particle_ref[i].zi, thres);
    assert_equal(part->pxi, particle_ref[i].pxi, thres);
    assert_equal(part->pyi, particle_ref[i].pyi, thres);
    assert_equal(part->pzi, particle_ref[i].pzi, thres);
  }
  printf("max delta: (%s)\n", test_str);
  printf("    xi ,yi ,zi  %g\t%g\t%g\n    pxi,pyi,pzi %g\t%g\t%g\n",
	 xi, yi, zi, pxi, pyi, pzi);
}


// ----------------------------------------------------------------------
// psc_check_fields_ref
//
// check field data against previously saved reference solution

void
psc_check_fields_ref(int *flds, double thres)
{
  assert(field_ref);
  for (int i = 0; flds[i] >= 0; i++) {
    int m = flds[i];
    foreach_patch(patch) {
      foreach_3d(patch, ix, iy, iz, 0, 0) {
	//	  printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	assert_equal(F3_BASE(&psc.pf, m, ix,iy,iz), F3_BASE(field_ref, m, ix,iy,iz), thres);
      } foreach_3d_end;
    }
  }
}

// ----------------------------------------------------------------------
// psc_check_currents_ref
//
// check current current density data agains previously saved reference solution

void
psc_check_currents_ref(double thres)
{
  assert(field_ref);
  foreach_patch(patch) {
    for (int m = JXI; m <= JZI; m++){
      foreach_3d_g(patch, ix, iy, iz) {
	double val = F3_BASE(&psc.pf, m, ix,iy,iz);
	if (fabs(val) > 0.) {
	  printf("cur %s: [%d,%d,%d] = %g\n", fldname[m],
		 ix, iy, iz, val);
	} foreach_3d_g_end;
      }
    }
  }
  foreach_patch(patch) {
    for (int m = JXI; m <= JZI; m++){
      double max_delta = 0.;
      foreach_3d_g(patch, ix, iy, iz) {
	//	  printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	assert_equal(F3_BASE(&psc.pf, m, ix,iy,iz), F3_BASE(field_ref,m, ix,iy,iz), thres);
	max_delta = fmax(max_delta, 
			 fabs(F3_BASE(&psc.pf, m, ix,iy,iz) - F3_BASE(field_ref, m, ix,iy,iz)));
      } foreach_3d_g_end;
      printf("max_delta (%s) %g\n", fldname[m], max_delta);
    }
  }
}

void
psc_check_currents_ref_noghost(double thres)
{
  assert(field_ref);
  for (int m = JXI; m <= JZI; m++){
    foreach_patch(patch) {
      foreach_3d(patch, ix, iy, iz, 0, 0) {
	//	  printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	assert_equal(F3_BASE(&psc.pf, m, ix,iy,iz), F3_BASE(field_ref, m, ix,iy,iz), thres);
      } foreach_3d_end;
    }
  }
}

// ----------------------------------------------------------------------
// psc_check_particles_sorted
//
// checks particles are sorted by cell index

void
psc_check_particles_sorted()
{
#if PARTICLES_BASE == PARTICLES_FORTRAN
  int last = INT_MIN;

  for (int i = 0; i < psc.pp.n_part; i++) {
    assert(psc.pp.particles[i].cni >= last);
    last = psc.pp.particles[i].cni;
  }
#else
  assert(0);
#endif
}

// ----------------------------------------------------------------------
// psc_create_test_xz

void
psc_create_test_xz(struct psc_mod_config *conf)
{
  // make sure if we call it again, we really get the same i.c.
  srandom(0);

  psc_create(conf);
  psc_init("test_xz");
}

// ----------------------------------------------------------------------
// psc_create_test_yz

void
psc_create_test_yz(struct psc_mod_config *conf)
{
  // make sure if we call it again, we really get the same i.c.
  srandom(0);

  psc_create(conf);
  psc_init("test_yz");
  psc.pp.n_part = 1;
  psc_sort();
}

