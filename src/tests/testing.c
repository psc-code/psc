
#include "psc_testing.h"
#include "psc_sort.h"
#include "psc_case.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"

#include <math.h>
#include <limits.h>
#include <stdlib.h>

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

static mparticles_t *particles_ref;
static struct psc_mfields *flds_ref;

// ----------------------------------------------------------------------
// psc_save_particles_ref
//
// save current particle data as reference solution

void
psc_save_particles_ref(struct psc *psc, mparticles_base_t *particles_base)
{
  mparticles_t *particles = psc_mparticles_get_cf(particles_base);

  if (!particles_ref) {
    int nr_particles_by_patch[psc->nr_patches];
    psc_foreach_patch(psc, p) {
      nr_particles_by_patch[p] = psc_mparticles_get_patch(particles, p)->n_part;
    }
    particles_ref = psc_mparticles_create(MPI_COMM_WORLD);
    psc_mparticles_set_domain_nr_particles(particles_ref, psc->mrc_domain,
					   nr_particles_by_patch);
    psc_mparticles_setup(particles_ref);
  }
  psc_foreach_patch(psc, p) {
    particles_t *pp = psc_mparticles_get_patch(particles, p);
    particles_t *pp_ref = psc_mparticles_get_patch(particles_ref, p);
    for (int i = 0; i < pp->n_part; i++) {
      *particles_get_one(pp_ref, i) = *particles_get_one(pp, i);
    }
  }

  psc_mparticles_put_cf(particles, particles_base);
}

// ----------------------------------------------------------------------
// psc_save_fields_ref
//
// save current field data as reference solution

void
psc_save_fields_ref(struct psc *psc, mfields_base_t *flds_base)
{
  int me = psc->domain.use_pml ? NR_FIELDS : HZ + 1;
  mfields_t *flds = psc_mfields_get_cf(flds_base, 0, me);

  if (!flds_ref) {
    flds_ref = psc_mfields_create(psc_comm(psc));
    psc_mfields_set_domain(flds_ref, psc->mrc_domain);
    psc_mfields_set_param_int(flds_ref, "nr_fields", NR_FIELDS);
    psc_mfields_set_param_int3(flds_ref, "ibn", psc->ibn);
    psc_mfields_setup(flds_ref);
  }

  psc_foreach_patch(psc, p) {
    fields_t *pf = psc_mfields_get_patch(flds, p);
    fields_t *pf_ref = psc_mfields_get_patch(flds_ref, p);
    for (int m = 0; m < me; m++) {
      psc_foreach_3d_g(psc, p, ix, iy, iz) {
	F3(pf_ref, m, ix,iy,iz) = F3(pf, m, ix,iy,iz);
      } psc_foreach_3d_g_end;
    }
  }
  psc_mfields_put_cf(flds, flds_base, 0, me);
} 

// ----------------------------------------------------------------------
// psc_check_particles_ref
//
// check current particle data agains previously saved reference solution

void
psc_check_particles_ref(struct psc *psc, mparticles_base_t *particles_base,
			double thres, const char *test_str)
{
  mparticles_t *particles = psc_mparticles_get_cf(particles_base);

  assert(particles_ref);
  particle_real_t xi = 0., yi = 0., zi = 0., pxi = 0., pyi = 0., pzi = 0.;
  psc_foreach_patch(psc, p) {
    particles_t *pp = psc_mparticles_get_patch(particles, p);
    particles_t *pp_ref = psc_mparticles_get_patch(particles_ref, p);
    
    for (int i = 0; i < pp->n_part; i++) {
      particle_t *part = particles_get_one(pp, i);
      particle_t *part_ref = particles_get_one(pp_ref, i);
      //    printf("i = %d\n", i);
      xi  = fmax(xi , fabs(part->xi  - part_ref->xi));
      yi  = fmax(yi , fabs(part->yi  - part_ref->yi));
      zi  = fmax(zi , fabs(part->zi  - part_ref->zi));
      pxi = fmax(pxi, fabs(part->pxi - part_ref->pxi));
      pyi = fmax(pyi, fabs(part->pyi - part_ref->pyi));
      pzi = fmax(pzi, fabs(part->pzi - part_ref->pzi));
      assert_equal(part->xi , part_ref->xi, thres);
      assert_equal(part->yi , part_ref->yi, thres);
      assert_equal(part->zi , part_ref->zi, thres);
      assert_equal(part->pxi, part_ref->pxi, thres);
      assert_equal(part->pyi, part_ref->pyi, thres);
      assert_equal(part->pzi, part_ref->pzi, thres);
    }
  }
  printf("max delta: (%s)\n", test_str);
  printf("    xi ,yi ,zi  %g\t%g\t%g\n    pxi,pyi,pzi %g\t%g\t%g\n",
	 xi, yi, zi, pxi, pyi, pzi);

  psc_mparticles_put_cf(particles, particles_base); // FIXME, no copy-back needed
}


// ----------------------------------------------------------------------
// psc_check_fields_ref
//
// check field data against previously saved reference solution

void
psc_check_fields_ref(struct psc *psc, mfields_base_t *flds_base, int *m_flds, double thres)
{
  mfields_t *flds = psc_mfields_get_cf(flds_base, 0, 12);

  psc_foreach_patch(psc, p) {
    fields_t *pf = psc_mfields_get_patch(flds, p);
    fields_t *pf_ref = psc_mfields_get_patch(flds_ref, p);
    for (int i = 0; m_flds[i] >= 0; i++) {
      int m = m_flds[i];
      psc_foreach_3d(psc, p, ix, iy, iz, 0, 0) {
	//	  printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	assert_equal(F3(pf, m, ix,iy,iz), F3(pf_ref, m, ix,iy,iz), thres);
      } psc_foreach_3d_end;
    }
  }
  psc_mfields_put_cf(flds, flds_base, 0, 0);
}

// ----------------------------------------------------------------------
// psc_check_currents_ref
//
// check current current density data agains previously saved reference solution

void
psc_check_currents_ref(struct psc *psc, mfields_base_t *flds_base, double thres)
{
  mfields_t *flds = psc_mfields_get_cf(flds_base, JXI, JXI + 3);

#if 0
  foreach_patch(p) {
    fields_base_t *pf = psc_mfields_get_patch(flds, p);
    for (int m = JXI; m <= JZI; m++){
      foreach_3d_g(p, ix, iy, iz) {
	double val = F3(pf, m, ix,iy,iz);
	if (fabs(val) > 0.) {
	  printf("cur %s: [%d,%d,%d] = %g\n", fldname[m],
		 ix, iy, iz, val);
	} foreach_3d_g_end;
      }
    }
  }
#endif
  mfields_t *diff = psc_mfields_create(psc_comm(psc));
  psc_mfields_set_domain(diff, psc->mrc_domain);
  psc_mfields_set_param_int3(diff, "ibn", psc->ibn);
  psc_mfields_setup(diff);
  // FIXME, make funcs for this (waxpy, norm)

  for (int m = JXI; m <= JZI; m++){
    fields_real_t max_delta = 0.;
    psc_foreach_patch(psc, p) {
      fields_t *pf = psc_mfields_get_patch(flds, p);
      fields_t *pf_ref = psc_mfields_get_patch(flds_ref, p);
      fields_t *pf_diff = psc_mfields_get_patch(diff, p);
      psc_foreach_3d_g(psc, p, ix, iy, iz) {
	F3(pf_diff, 0, ix,iy,iz) =
	  F3(pf, m, ix,iy,iz) - F3(pf_ref, m, ix,iy,iz);
	max_delta = fmax(max_delta, fabs(F3(pf_diff, 0, ix,iy,iz)));
      } psc_foreach_3d_g_end;
    }
    printf("max_delta (%s) %g / thres %g\n", fldname[m], max_delta, thres);
    if (max_delta > thres) {
      //      psc_dump_field(diff, 0, "diff");
      abort();
    }
  }
  psc_mfields_destroy(diff);

  psc_mfields_put_cf(flds, flds_base, 0, 0);
}

void
psc_check_currents_ref_noghost(struct psc *psc, mfields_base_t *flds_base, double thres)
{
  mfields_t *flds = psc_mfields_get_cf(flds_base, JXI, JXI + 3);

  psc_foreach_patch(psc, p) {
    fields_t *pf = psc_mfields_get_patch(flds, p);
    fields_t *pf_ref = psc_mfields_get_patch(flds_ref, p);
    for (int m = JXI; m <= JZI; m++){
      psc_foreach_3d(psc, p, ix, iy, iz, 0, 0) {
	//	  printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	assert_equal(F3(pf, m, ix,iy,iz), F3(pf_ref, m, ix,iy,iz), thres);
      } psc_foreach_3d_end;
    }
  }

  psc_mfields_put_cf(flds, flds_base, 0, 0);
}

// ----------------------------------------------------------------------
// psc_check_particles_sorted
//
// checks particles are sorted by cell index

void
psc_check_particles_sorted(struct psc *psc, mparticles_base_t *particles_base)
{
  mparticles_t *particles = psc_mparticles_get_cf(particles_base);

  int last = INT_MIN;

  particle_real_t dxi = 1.f / psc->dx[0];
  particle_real_t dyi = 1.f / psc->dx[1];
  particle_real_t dzi = 1.f / psc->dx[2];

  int *ibn = psc->ibn;
  psc_foreach_patch(psc, p) {
    particles_t *pp = psc_mparticles_get_patch(particles, p);
    struct psc_patch *patch = &psc->patch[p];
    int *ldims = patch->ldims;
    for (int n = 0; n < pp->n_part; n++) {
      // FIXME, duplicated
      particle_t *part = particles_get_one(pp, n);

      particle_real_t u = (part->xi - patch->xb[0]) * dxi;
      particle_real_t v = (part->yi - patch->xb[1]) * dyi;
      particle_real_t w = (part->zi - patch->xb[2]) * dzi;
      int j0 = particle_real_nint(u) + ibn[0];
      int j1 = particle_real_nint(v) + ibn[1];
      int j2 = particle_real_nint(w) + ibn[2];
      
      assert(j0 >= 0 && j0 < ldims[0] + 2*ibn[0]);
      assert(j1 >= 0 && j1 < ldims[1] + 2*ibn[1]);
      assert(j2 >= 0 && j2 < ldims[2] + 2*ibn[2]);
      int cni = ((j2) * (ldims[1] + 2*ibn[1]) + j1) * (ldims[0] + 2*ibn[0]) + j0;

      assert(cni >= last);
      last = cni;
    }
  }
  psc_mparticles_put_cf(particles, particles_base);
}

// ----------------------------------------------------------------------
// psc_create_test_xy

struct psc_case *
psc_create_test_xy()
{
  // make sure if we call it again, we really get the same i.c.
  srandom(0);

  struct psc_case *_case = psc_case_create(MPI_COMM_WORLD);
  psc_case_set_type(_case, "test_xy");
  psc_case_set_from_options(_case);
  return _case;
}


// ----------------------------------------------------------------------
// psc_create_test_xz

struct psc_case *
psc_create_test_xz()
{
  // make sure if we call it again, we really get the same i.c.
  srandom(0);

  struct psc_case *_case = psc_case_create(MPI_COMM_WORLD);
  psc_case_set_type(_case, "test_xz");
  psc_case_set_from_options(_case);
  return _case;
}

// ----------------------------------------------------------------------
// psc_create_test_yz

struct psc_case *
psc_create_test_yz(void)
{
  // make sure if we call it again, we really get the same i.c.
  srandom(0);

  struct psc_case *_case = psc_case_create(MPI_COMM_WORLD);
  psc_case_set_type(_case, "test_yz");
  psc_case_set_from_options(_case);
  return _case;
}

// ----------------------------------------------------------------------
// psc_create_test_z

struct psc_case *
psc_create_test_z(void)
{
  // make sure if we call it again, we really get the same i.c.
  srandom(0);

  struct psc_case *_case = psc_case_create(MPI_COMM_WORLD);
  psc_case_set_type(_case, "test_z");
  psc_case_set_from_options(_case);
  return _case;
}

