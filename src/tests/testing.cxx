
#include "psc_testing.h"
#include "psc_sort.h"
#include "psc_particles_as_double.h"
#include "psc_fields_as_c.h"
#include "psc_push_particles.h"
#include "psc_bnd.h"
#include "psc_bnd_particles.h"
#include "psc_randomize.h"
#include "psc_output_fields_item.h"

#include "fields.hxx"

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <limits.h>
#include <stdlib.h>

using Fields = Fields3d<fields_t>;
using real_t = Fields::real_t;

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
  psc_mfields_set_param_obj(fld, "domain", psc->mrc_domain);
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

  mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  mparticles_t mp_ref = mparticles_t(mprts_ref);
  psc_foreach_patch(psc, p) {
    int n_prts = mprts[p].size();
    for (int n = 0; n < n_prts; n++) {
      mp_ref[p].push_back(mprts[p][n]);
    }
  }
  mprts.put_as(mprts_base, MP_DONT_COPY);
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
    psc_mfields_set_param_obj(mflds_ref, "domain", psc->mrc_domain);
    psc_mfields_set_param_int(mflds_ref, "nr_fields", NR_FIELDS);
    psc_mfields_set_param_int3(mflds_ref, "ibn", psc->ibn);
    psc_mfields_setup(mflds_ref);
  }

  mfields_t mf = mflds_base->get_as<mfields_t>(0, me);
  mfields_t mf_ref(mflds_ref);
  psc_foreach_patch(psc, p) {
    Fields F = mf[p];
    Fields R = mf_ref[p];
    for (int m = 0; m < me; m++) {
      psc_foreach_3d_g(psc, p, ix, iy, iz) {
	R(m, ix,iy,iz) = F(m, ix,iy,iz);
      } psc_foreach_3d_g_end;
    }
  }
  mf.put_as(mflds_base, 0, 0);
} 

// ----------------------------------------------------------------------
// psc_check_particles_ref
//
// check current particle data agains previously saved reference solution

void
psc_check_particles_ref(struct psc *psc, struct psc_mparticles *mprts_base,
			double thres, const char *test_str)
{
  if (!opt_testing_check_particles) {
    return;
  }

  assert(mprts_ref);
  mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  real_t xi = 0., yi = 0., zi = 0., pxi = 0., pyi = 0., pzi = 0.;
  psc_foreach_patch(psc, p) {
    auto& prts = mprts[p];
    auto& prts_ref = mparticles_t(mprts_ref)[p];
  
    assert(prts.size() == prts_ref.size());
    for (auto prt_iter = prts.begin(), prt_ref_iter = prts_ref.begin();
	 prt_iter != prts.end(); ++prt_iter, ++prt_ref_iter) {
      *prt_ref_iter = *prt_iter;
      particle_t *part = &*prt_iter;
      particle_t *part_ref = &*prt_ref_iter;
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
  mprts.put_as(mprts_base, 0);

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
  mfields_t mf = mflds_base->get_as<mfields_t>(0, 12);
  mfields_t mf_ref(mflds_ref);
  psc_foreach_patch(psc, p) {
    Fields F(mf[p]);
    Fields R(mf_ref[p]);
    for (int i = 0; m_flds[i] >= 0; i++) {
      int m = m_flds[i];
      psc_foreach_3d(psc, p, ix, iy, iz, 0, 0) {
	//	  printf("m %d %d,%d,%d\n", m, ix,iy,iz);
	assert_equal(F(m, ix,iy,iz), R(m, ix,iy,iz), thres);
      } psc_foreach_3d_end;
    }
  }
  mf.put_as(mflds_base, 0, 0);
}

// ----------------------------------------------------------------------
// psc_check_currents_ref
//
// check current current density data agains previously saved reference solution

void
psc_check_currents_ref(struct psc *psc, struct psc_mfields *mflds_base, double thres, int sw)
{
#if 0
  mfields_t mf = mflds_base->get_as<mfields_t>(JXI, JXI+3);
  foreach_patch(p) {
    Fields F(mf[p]);
    for (int m = JXI; m <= JZI; m++){
      foreach_3d_g(p, ix, iy, iz) {
	double val = F(m, ix,iy,iz);
	if (fabs(val) > 0.) {
	  printf("cur %s: [%d,%d,%d] = %g\n", fldname[m],
		 ix, iy, iz, val);
	} foreach_3d_g_end;
      }
    }
  }
  psc_mfields_put_as(mflds, mflds_base, 0, 0);
  struct psc_mfields *diff = psc_mfields_create(psc_comm(psc));
  psc_mfields_set_param_obj(diff, "domain", psc->mrc_domain);
  psc_mfields_set_param_int3(diff, "ibn", psc->ibn);
  psc_mfields_setup(diff);
  // FIXME, make funcs for this (waxpy, norm)

  mfields_t mf_ref(mflds_ref), mf_diff(diff);
  mfields_t mf = mflds_base->get_as<mfields_t>(JXI, JXI + 3);
  for (int m = JXI; m <= JZI; m++){
    fields_t::real_t max_delta = 0.;
    psc_foreach_patch(psc, p) {
      Fields F(mf[p]);
      Fields R(mf_ref[p]);
      Fields D(mf_diff[p]);
      psc_foreach_3d(psc, p, ix, iy, iz, sw, sw) {
	D(0, ix,iy,iz) = F(m, ix,iy,iz) - R(m, ix,iy,iz);
	max_delta = fmax(max_delta, fabs(D(0, ix,iy,iz)));
      } psc_foreach_3d_g_end;
    }
    if (opt_checks_verbose || max_delta > thres) {
      mprintf("max_delta (%s) %g / thres %g\n",
	      psc_mfields_comp_name(mflds, m), max_delta, thres);
    }
    if (max_delta > thres) {
      assert(0);
    }
  }
  mf.put_as(mflds_base, 0, 0);
  psc_mfields_destroy(diff);
#endif
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
  dens->axpy(-1., dens_ref);
  mfields_t mf_dens(dens);

  // FIXME, do this generically
  for (int m = 0; m < 2; m++) {
    double max_err = 0.;
    psc_foreach_patch(psc, p) {
      Fields D(mf_dens[p]);
      psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
	fields_t::real_t val = D(m, jx,jy,jz);
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
psc_check_particles_sorted(struct psc *psc, struct psc_mparticles *mprts_base)
{
  int last = INT_MIN;

  mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  int *ibn = psc->ibn;
  psc_foreach_patch(psc, p) {
    auto& prts = mprts[p];
    real_t dxi = 1.f / psc->grid.dx[0];
    real_t dyi = 1.f / psc->grid.dx[1];
    real_t dzi = 1.f / psc->grid.dx[2];

    int *ldims = psc->grid.ldims;

    PARTICLE_ITER_LOOP(prt_iter, prts.begin(), prts.end()) {
      particle_t *part = &*prt_iter;
      // FIXME, duplicated

      real_t u = part->xi * dxi;
      real_t v = part->yi * dyi;
      real_t w = part->zi * dzi;
      int j0 = nint(u) + ibn[0];
      int j1 = nint(v) + ibn[1];
      int j2 = nint(w) + ibn[2];
      
      assert(j0 >= 0 && j0 < ldims[0] + 2*ibn[0]);
      assert(j1 >= 0 && j1 < ldims[1] + 2*ibn[1]);
      assert(j2 >= 0 && j2 < ldims[2] + 2*ibn[2]);
      int cni = ((j2) * (ldims[1] + 2*ibn[1]) + j1) * (ldims[0] + 2*ibn[0]) + j0;

      assert(cni >= last);
      last = cni;
    }
  }
  mprts.put_as(mprts_base, 0);
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
psc_testing_push_particles(struct psc *psc, const char *s_push_particles)
{
  if (opt_checks_verbose) {
    mprintf("=== testing push_part() %s\n", s_push_particles);
  }

  psc_randomize_run(psc->randomize, psc->particles);
  psc_bnd_particles_exchange(psc->bnd_particles, psc->particles);
  psc_sort_run(psc->sort, psc->particles);

  //  psc_testing_dump(psc, s_push_particles);

  psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);
  psc_bnd_particles_exchange(psc->bnd_particles, psc->particles);
  psc_sort_run(psc->sort, psc->particles);
  psc_bnd_add_ghosts(psc->bnd, psc->flds, JXI, JXI + 3);

  //  psc_testing_dump(psc, s_push_particles);
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

