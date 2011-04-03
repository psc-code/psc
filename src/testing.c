
#include "psc_testing.h"
#include "psc_sort.h"
#include "psc_case.h"

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

static mparticles_base_t particles_ref;
static mfields_base_t flds_ref;

// ----------------------------------------------------------------------
// psc_save_particles_ref
//
// save current particle data as reference solution

void
psc_save_particles_ref(struct psc *psc, mparticles_base_t *particles)
{
  if (!particles_ref.p) {
    particles_ref.p = calloc(psc->nr_patches, sizeof(*particles_ref.p));
    psc_foreach_patch(psc, p) {
      particles_base_t *pp = &particles->p[p];
      particles_base_alloc(&particles_ref.p[p], pp->n_part);
    }
  }
  psc_foreach_patch(psc, p) {
    particles_base_t *pp = &particles->p[p];
    particles_base_t *pp_ref = &particles_ref.p[p];
    for (int i = 0; i < pp->n_part; i++) {
      *particles_base_get_one(pp_ref, i) = *particles_base_get_one(pp, i);
    }
  }
}

// ----------------------------------------------------------------------
// psc_save_fields_ref
//
// save current field data as reference solution

void
psc_save_fields_ref(struct psc *psc, mfields_base_t *flds)
{
  if (!flds_ref.f) {
    mfields_base_alloc(psc->mrc_domain, &flds_ref, NR_FIELDS, psc->ibn);
  }
  int me = psc->domain.use_pml ? NR_FIELDS : HZ + 1;
  psc_foreach_patch(psc, p) {
    fields_base_t *pf = &flds->f[p];
    fields_base_t *pf_ref = &flds_ref.f[p];
    for (int m = 0; m < me; m++) {
      psc_foreach_3d_g(psc, p, ix, iy, iz) {
	F3_BASE(pf_ref, m, ix,iy,iz) = F3_BASE(pf, m, ix,iy,iz);
      } psc_foreach_3d_g_end;
    }
  }
} 

// ----------------------------------------------------------------------
// psc_check_particles_ref
//
// check current particle data agains previously saved reference solution

void
psc_check_particles_ref(struct psc *psc, mparticles_base_t *particles,
			double thres, const char *test_str)
{
  assert(particles_ref.p);
  particle_base_real_t xi = 0., yi = 0., zi = 0., pxi = 0., pyi = 0., pzi = 0.;
  psc_foreach_patch(psc, p) {
    particles_base_t *pp = &particles->p[p];
    particles_base_t *pp_ref = &particles_ref.p[p];
    
    for (int i = 0; i < pp->n_part; i++) {
      particle_base_t *part = particles_base_get_one(pp, i);
      particle_base_t *part_ref = particles_base_get_one(pp_ref, i);
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
}


// ----------------------------------------------------------------------
// psc_check_fields_ref
//
// check field data against previously saved reference solution

void
psc_check_fields_ref(struct psc *psc, mfields_base_t *flds, int *m_flds, double thres)
{
  psc_foreach_patch(psc, p) {
    fields_base_t *pf = &flds->f[p];
    fields_base_t *pf_ref = &flds_ref.f[p];
    for (int i = 0; m_flds[i] >= 0; i++) {
      int m = m_flds[i];
      psc_foreach_3d(psc, p, ix, iy, iz, 0, 0) {
	//	  printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	assert_equal(F3_BASE(pf, m, ix,iy,iz), F3_BASE(pf_ref, m, ix,iy,iz), thres);
      } psc_foreach_3d_end;
    }
  }
}

// ----------------------------------------------------------------------
// psc_check_currents_ref
//
// check current current density data agains previously saved reference solution

void
psc_check_currents_ref(struct psc *psc, mfields_base_t *flds, double thres)
{
#if 0
  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    for (int m = JXI; m <= JZI; m++){
      foreach_3d_g(p, ix, iy, iz) {
	double val = F3_BASE(pf, m, ix,iy,iz);
	if (fabs(val) > 0.) {
	  printf("cur %s: [%d,%d,%d] = %g\n", fldname[m],
		 ix, iy, iz, val);
	} foreach_3d_g_end;
      }
    }
  }
#endif
  psc_foreach_patch(psc, p) {
    fields_base_t *pf = &flds->f[p];
    fields_base_t *pf_ref = &flds_ref.f[p];
    for (int m = JXI; m <= JZI; m++){
      double max_delta = 0.;
      psc_foreach_3d_g(psc, p, ix, iy, iz) {
	//	  printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	assert_equal(F3_BASE(pf, m, ix,iy,iz), F3_BASE(pf_ref,m, ix,iy,iz), thres);
	max_delta = fmax(max_delta, 
			 fabs(F3_BASE(pf, m, ix,iy,iz) - F3_BASE(pf_ref, m, ix,iy,iz)));
      } psc_foreach_3d_g_end;
      printf("max_delta (%s) %g\n", fldname[m], max_delta);
    }
  }
}

void
psc_check_currents_ref_noghost(struct psc *psc, mfields_base_t *flds, double thres)
{
  psc_foreach_patch(psc, p) {
    fields_base_t *pf = &flds->f[p];
    fields_base_t *pf_ref = &flds_ref.f[p];
    for (int m = JXI; m <= JZI; m++){
      psc_foreach_3d(psc, p, ix, iy, iz, 0, 0) {
	//	  printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	assert_equal(F3_BASE(pf, m, ix,iy,iz), F3_BASE(pf_ref, m, ix,iy,iz), thres);
      } psc_foreach_3d_end;
    }
  }
}

// ----------------------------------------------------------------------
// psc_check_particles_sorted
//
// checks particles are sorted by cell index

void
psc_check_particles_sorted(struct psc *psc, mparticles_base_t *particles)
{
#if PARTICLES_BASE == PARTICLES_FORTRAN
  int last = INT_MIN;

  psc_foreach_patch(psc, p) {
    particles_base_t *pp = &particles->p[p];
    for (int i = 0; i < pp->n_part; i++) {
      assert(pp->particles[i].cni >= last);
      last = pp->particles[i].cni;
    }
  }
#else
  int last = INT_MIN;

  particle_base_real_t dxi = 1.f / psc->dx[0];
  particle_base_real_t dyi = 1.f / psc->dx[1];
  particle_base_real_t dzi = 1.f / psc->dx[2];

  int *ibn = psc->ibn;
  psc_foreach_patch(psc, p) {
    particles_base_t *pp = &particles->p[p];
    struct psc_patch *patch = &psc->patch[p];
    int *ldims = patch->ldims;
    for (int n = 0; n < pp->n_part; n++) {
      // FIXME, duplicated
      particle_base_t *part = particles_base_get_one(pp, n);

      particle_base_real_t u = (part->xi - patch->xb[0]) * dxi;
      particle_base_real_t v = (part->yi - patch->xb[1]) * dyi;
      particle_base_real_t w = (part->zi - patch->xb[2]) * dzi;
      int j0 = particle_base_real_nint(u) + ibn[0];
      int j1 = particle_base_real_nint(v) + ibn[1];
      int j2 = particle_base_real_nint(w) + ibn[2];
      
      assert(j0 >= 0 && j0 < ldims[0] + 2*ibn[0]);
      assert(j1 >= 0 && j1 < ldims[1] + 2*ibn[1]);
      assert(j2 >= 0 && j2 < ldims[2] + 2*ibn[2]);
      int cni = ((j2) * (ldims[1] + 2*ibn[1]) + j1) * (ldims[0] + 2*ibn[0]) + j0;

      assert(cni >= last);
      last = cni;
    }
  }
#endif
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

void
psc_create_test_yz(struct psc_mod_config *conf)
{
  // make sure if we call it again, we really get the same i.c.
  srandom(0);

  struct psc_case *_case = psc_case_create(MPI_COMM_WORLD);
  psc_case_set_type(_case, "test_yz");
  struct psc *psc = psc_case_get_psc(_case);
  psc_set_conf(psc, conf);
  psc_case_set_from_options(_case);
  psc_case_setup(_case);
  psc_case_view(_case);
  psc->particles.p[0].n_part = 1;
  psc_sort_run(psc->sort, &psc->particles);
}

