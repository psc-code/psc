
#include "psc_testing.h"
#include "psc_sort.h"
#include "psc_particles_as_c.h"
#include "psc_fields_as_c.h"
#include "psc_push_particles.h"
#include "psc_bnd.h"
#include "psc_bnd_particles.h"
#include "psc_randomize.h"
#include "psc_output_fields_item.h"

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <limits.h>
#include <stdlib.h>

bool opt_checks_verbose = false;
bool opt_testing_dump = false;
bool opt_testing_check_currents = true;
bool opt_testing_check_particles = true;
bool opt_testing_check_densities = true;

// ----------------------------------------------------------------------
// FIXME, should be consolidated?

static struct psc_mfields *
fld_create(struct psc *psc, int nr_fields)
{
  struct psc_mfields *fld = psc_mfields_create(psc_comm(psc));
  psc_mfields_set_type(fld, "c");
  psc_mfields_set_domain(fld, psc->mrc_domain);
  psc_mfields_set_param_int3(fld, "ibn", psc->ibn);
  psc_mfields_set_param_int(fld, "nr_fields", nr_fields);
  psc_mfields_setup(fld);

  return fld;
}

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

static struct psc_mparticles *mprts_ref;
static struct psc_mfields *mflds_ref;

// ----------------------------------------------------------------------
// psc_save_particles_ref
//
// save current particle data as reference solution

void
psc_save_particles_ref(struct psc *psc, struct psc_mparticles *mprts_base)
{
  if (!mprts_ref) {
    int n_prts_by_patch[psc->nr_patches];
    psc_mparticles_get_size_all(mprts_base, n_prts_by_patch);

    mprts_ref = psc_mparticles_create(MPI_COMM_WORLD);
    psc_mparticles_set_param_int(mprts_ref, "nr_patches", psc->nr_patches);
    psc_mparticles_setup(mprts_ref);
    psc_mparticles_reserve_all(mprts_ref, n_prts_by_patch);
  }

  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, "c", 0);
  psc_foreach_patch(psc, p) {
    int n_prts = mparticles_get_n_prts(mprts, p);
    for (int n = 0; n < n_prts; n++) {
      mparticles_patch_push_back(mprts_ref, p, *mparticles_get_one(mprts, p, n));
    }
  }
  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);
}

// ----------------------------------------------------------------------
// psc_save_fields_ref
//
// save current field data as reference solution

void
psc_save_fields_ref(struct psc *psc, struct psc_mfields *mflds_base)
{
  int me = HZ + 1;

  if (!mflds_ref) {
    mflds_ref = psc_mfields_create(psc_comm(psc));
    psc_mfields_set_domain(mflds_ref, psc->mrc_domain);
    psc_mfields_set_param_int(mflds_ref, "nr_fields", NR_FIELDS);
    psc_mfields_set_param_int3(mflds_ref, "ibn", psc->ibn);
    psc_mfields_setup(mflds_ref);
  }

  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "c", 0, me);
  psc_foreach_patch(psc, p) {
    struct psc_fields *pf = psc_mfields_get_patch(mflds, p);
    struct psc_fields *pf_ref = psc_mfields_get_patch(mflds_ref, p);
    for (int m = 0; m < me; m++) {
      psc_foreach_3d_g(psc, p, ix, iy, iz) {
	F3(pf_ref, m, ix,iy,iz) = F3(pf, m, ix,iy,iz);
      } psc_foreach_3d_g_end;
    }
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
} 

// ----------------------------------------------------------------------
// psc_check_particles_ref
//
// check current particle data agains previously saved reference solution

void
psc_check_particles_ref(struct psc *psc, struct psc_mparticles *particles_base,
			double thres, const char *test_str)
{
  if (!opt_testing_check_particles) {
    return;
  }

  assert(mprts_ref);
  struct psc_mparticles *mprts = psc_mparticles_get_as(particles_base, "c", 0);
  particle_real_t xi = 0., yi = 0., zi = 0., pxi = 0., pyi = 0., pzi = 0.;
  psc_foreach_patch(psc, p) {
    particle_range_t prts = particle_range_mprts(mprts, p);
    particle_range_t prts_ref = particle_range_mprts(mprts_ref, p);
  
    assert(particle_range_size(prts) == particle_range_size(prts_ref));
    for (particle_iter_t prt_iter = prts.begin, prt_ref_iter = prts_ref.end;
	 !particle_iter_equal(prt_iter, prts.end);
	 prt_iter = particle_iter_next(prt_iter), prt_ref_iter = particle_iter_next(prt_ref_iter)) {
      *particle_iter_deref(prt_ref_iter) = *particle_iter_deref(prt_iter);
      particle_t *part = particle_iter_deref(prt_iter);
      particle_t *part_ref = particle_iter_deref(prt_ref_iter);
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
  psc_mparticles_put_as(mprts, particles_base, MP_DONT_COPY);

  if (opt_checks_verbose) {
    mprintf("max delta: (%s)\n", test_str);
    mprintf("    xi ,yi ,zi  %g\t%g\t%g\n    pxi,pyi,pzi %g\t%g\t%g\n",
	    xi, yi, zi, pxi, pyi, pzi);
  }
}


// ----------------------------------------------------------------------
// psc_check_fields_ref
//
// check field data against previously saved reference solution

void
psc_check_fields_ref(struct psc *psc, struct psc_mfields *mflds_base, int *m_flds, double thres)
{
   struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "c", 0, 12);
   psc_foreach_patch(psc, p) {
    struct psc_fields *pf = psc_mfields_get_patch(mflds, p);
    struct psc_fields *pf_ref = psc_mfields_get_patch(mflds_ref, p);
    for (int i = 0; m_flds[i] >= 0; i++) {
      int m = m_flds[i];
      psc_foreach_3d(psc, p, ix, iy, iz, 0, 0) {
	//	  printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	assert_equal(F3(pf, m, ix,iy,iz), F3(pf_ref, m, ix,iy,iz), thres);
      } psc_foreach_3d_end;
    }
  }
   psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

// ----------------------------------------------------------------------
// psc_check_currents_ref
//
// check current current density data agains previously saved reference solution

void
psc_check_currents_ref(struct psc *psc, struct psc_mfields *mflds_base, double thres, int sw)
{
#if 0
  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "c", JXI, JXI + 3);
  foreach_patch(p) {
    struct psc_fields *pf = psc_mfields_get_patch(mflds, p);
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
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
#endif
  struct psc_mfields *diff = psc_mfields_create(psc_comm(psc));
  psc_mfields_set_domain(diff, psc->mrc_domain);
  psc_mfields_set_param_int3(diff, "ibn", psc->ibn);
  psc_mfields_setup(diff);
  // FIXME, make funcs for this (waxpy, norm)

  struct psc_mfields *mflds = psc_mfields_get_as(mflds_base, "c", JXI, JXI + 3);
  for (int m = JXI; m <= JZI; m++){
    fields_real_t max_delta = 0.;
    psc_foreach_patch(psc, p) {
      struct psc_fields *pf = psc_mfields_get_patch(mflds, p);
      struct psc_fields *pf_ref = psc_mfields_get_patch(mflds_ref, p);
      struct psc_fields *pf_diff = psc_mfields_get_patch(diff, p);
      psc_foreach_3d(psc, p, ix, iy, iz, sw, sw) {
	F3(pf_diff, 0, ix,iy,iz) =
	  F3(pf, m, ix,iy,iz) - F3(pf_ref, m, ix,iy,iz);
	max_delta = fmax(max_delta, fabs(F3(pf_diff, 0, ix,iy,iz)));
      } psc_foreach_3d_g_end;
    }
    if (opt_checks_verbose || max_delta > thres) {
      mprintf("max_delta (%s) %g / thres %g\n", fldname[m], max_delta, thres);
    }
    if (max_delta > thres) {
      //      psc_dump_field(diff, 0, "diff");
      assert(0);
    }
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
  psc_mfields_destroy(diff);
}

// ----------------------------------------------------------------------
// psc_testing_check_densities_ref

void
psc_testing_check_densities_ref(struct psc *psc, struct psc_mparticles *particles,
				double eps)
{
  struct psc_mfields *dens_ref = fld_create(psc, 3);
  struct psc_mfields *dens = fld_create(psc, 3);

  struct psc_output_fields_item *item = psc_output_fields_item_create(psc_comm(psc));
  psc_output_fields_item_set_type(item, "n");
  psc_output_fields_item_set_psc_bnd(item, psc->bnd);
  psc_output_fields_item_setup(item);
  psc_output_fields_item_run(item, mflds_ref, mprts_ref, dens_ref);
  psc_output_fields_item_run(item, psc->flds, particles, dens);
  psc_output_fields_item_destroy(item);

  // dens -= dens_ref
  psc_mfields_axpy(dens, -1., dens_ref);

  // FIXME, do this generically
  for (int m = 0; m < 2; m++) {
    double max_err = 0.;
    psc_foreach_patch(psc, p) {
      struct psc_fields *p_dens = psc_mfields_get_patch(dens, p);
      psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
	fields_c_real_t val = F3_C(p_dens,m, jx,jy,jz);
	max_err = fmax(max_err, fabs(val));
	if (fabs(val) > eps) {
	  printf("(%d,%d,%d): diff %g\n", jx, jy, jz, val);
	}
      } psc_foreach_3d_end;
    }
    if (opt_checks_verbose) {
      mprintf("densities: m = %d max_err = %g (thres %g)\n", m, max_err, eps);
    }
    assert(max_err <= eps);
  }
}

// ----------------------------------------------------------------------
// psc_check_particles_sorted
//
// checks particles are sorted by cell index

void
psc_check_particles_sorted(struct psc *psc, struct psc_mparticles *particles_base)
{
  int last = INT_MIN;

  struct psc_mparticles *mprts = psc_mparticles_get_as(particles_base, "c", 0);
  int *ibn = psc->ibn;
  psc_foreach_patch(psc, p) {
    struct psc_patch *patch = &psc->patch[p];
    particle_range_t prts = particle_range_mprts(mprts, p);
    particle_real_t dxi = 1.f / patch->dx[0];
    particle_real_t dyi = 1.f / patch->dx[1];
    particle_real_t dzi = 1.f / patch->dx[2];

    int *ldims = patch->ldims;

    PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
      particle_t *part = particle_iter_deref(prt_iter);
      // FIXME, duplicated

      particle_real_t u = part->xi * dxi;
      particle_real_t v = part->yi * dyi;
      particle_real_t w = part->zi * dzi;
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
  psc_mparticles_put_as(mprts, particles_base, 0);
}

// ======================================================================

struct psc *
psc_testing_create_test_yz(const char *s_push_particles, unsigned int mask)
{
  struct psc *psc = psc_create(MPI_COMM_WORLD);
  psc->domain.gdims[0] = 1; // make yz
  psc_push_particles_set_type(psc->push_particles, s_push_particles);
  psc_sort_set_type(psc->sort, "countsort2");
  psc_sort_set_param_int(ppsc->sort, "mask", mask);
  psc_randomize_set_type(psc->randomize, "c");
  psc_set_from_options(psc);

  return psc;
}

struct psc *
psc_testing_create_test_xz()
{
  struct psc *psc = psc_create(MPI_COMM_WORLD);
  psc->domain.gdims[1] = 1; // make xz
  psc_set_from_options(psc);

  return psc;
}

void
psc_testing_dump(struct psc *psc, const char *basename)
{
  if (!opt_testing_dump)
    return;

  static int cnt = 0;

  char s[200];
  sprintf(s, "part_%s_%d", basename, cnt);
  psc_dump_particles(psc->particles, s);
  sprintf(s, "jx_%s_%d", basename, cnt);
  psc_dump_field(psc->flds, JXI, s);
  sprintf(s, "jy_%s_%d", basename, cnt);
  psc_dump_field(psc->flds, JYI, s);
  sprintf(s, "jz_%s_%d", basename, cnt);
  psc_dump_field(psc->flds, JZI, s);

  cnt ++;
}

void
psc_testing_push_particles(struct psc *psc, const char *s_push_particles)
{
  if (opt_checks_verbose) {
    mprintf("=== testing push_part() %s\n", s_push_particles);
  }

  psc_randomize_run(psc->randomize, psc->particles);
  psc_bnd_particles_exchange(psc->bnd_particles, psc->particles);
  psc_sort_run(psc->sort, psc->particles);

  psc_testing_dump(psc, s_push_particles);

  psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);
  psc_bnd_particles_exchange(psc->bnd_particles, psc->particles);
  psc_sort_run(psc->sort, psc->particles);
  psc_bnd_add_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);

  psc_testing_dump(psc, s_push_particles);
}

void
psc_testing_save_ref(struct psc *psc)
{
  psc_save_particles_ref(psc, psc->particles);
  psc_save_fields_ref(psc, psc->flds);
}

void
psc_testing_push_particles_check(struct psc *psc, double eps_particles, double eps_fields)
{
  //psc_check_continuity(psc, psc->particles, psc->flds, eps_fields);
  psc_check_particles_ref(psc, psc->particles, eps_particles, "push_particles");
  if (opt_testing_check_densities) {
    psc_testing_check_densities_ref(psc, psc->particles, eps_particles);
  }
  if (opt_testing_check_currents) {
    psc_check_currents_ref(psc, psc->flds, eps_fields, 0);
  }
}

// ======================================================================

void
psc_testing_init(int *argc, char ***argv)
{
  MPI_Init(argc, argv);
  libmrc_params_init(*argc, *argv);
  mrc_set_flags(MRC_FLAG_SUPPRESS_UNPREFIXED_OPTION_WARNING);

  mrc_params_get_option_bool("verbose", &opt_checks_verbose);
  mrc_params_get_option_bool("dump", &opt_testing_dump);
  mrc_params_get_option_bool("check_currents", &opt_testing_check_currents);
  mrc_params_get_option_bool("check_particles", &opt_testing_check_particles);
  mrc_params_get_option_bool("check_densities", &opt_testing_check_densities);
}

void
psc_testing_finalize()
{
  if (opt_checks_verbose) {
    prof_print();
  }

  MPI_Finalize();
}

