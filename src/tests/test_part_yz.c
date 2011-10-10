
#include "psc_testing.h"
#include "psc_push_particles.h"
#include "psc_sort.h"
#include "psc_bnd.h"
#include "psc_moments.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <string.h>
#include <mpi.h>

static bool do_dump = false;
static bool check_currents = true;

// ======================================================================
// psc_calc_rho

typedef particle_c_real_t creal;

static inline creal
creal_sqrt(creal x)
{
  return sqrt(x);
}

static void
do_shift_particle_positions(particles_c_t *pp, double dt)
{
  for (int n = 0; n < pp->n_part; n++) {
    particle_c_t *part = particles_c_get_one(pp, n);
      
    creal root = 1.f / creal_sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
    creal vxi = part->pxi * root;
    creal vyi = part->pyi * root;
    creal vzi = part->pzi * root;
    part->xi += vxi * dt;
    part->yi += vyi * dt;
    part->zi += vzi * dt;
  }
}

static void
psc_shift_particle_positions(struct psc *psc, mparticles_base_t *particles_base,
			     double dt)
{
  mparticles_c_t *particles = psc_mparticles_get_c(particles_base, 0);

  psc_foreach_patch(psc, p) {
    do_shift_particle_positions(psc_mparticles_get_patch_c(particles, p), dt);
  }

  psc_mparticles_put_c(particles, particles_base);
}

static void
dump(const char *basename, int cnt)
{
  if (!do_dump)
    return;

  char s[200];
  sprintf(s, "part_%s_%d", basename, cnt);
  psc_dump_particles(ppsc->particles, s);
  sprintf(s, "jx_%s_%d", basename, cnt);
  psc_dump_field(ppsc->flds, JXI, s);
  sprintf(s, "jy_%s_%d", basename, cnt);
  psc_dump_field(ppsc->flds, JYI, s);
  sprintf(s, "jz_%s_%d", basename, cnt);
  psc_dump_field(ppsc->flds, JZI, s);
}

static struct psc_case *
create_test(const char *s_push_particles)
{
  struct psc_case *_case = psc_create_test_yz();
  psc_push_particles_set_type(ppsc->push_particles, s_push_particles);
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_sort_set_param_int3(ppsc->sort, "blocksize", (int [3]) { 1, 8, 8 }); // FIXME
  psc_case_setup(_case);
  psc_bnd_exchange_particles(ppsc->bnd, ppsc->particles);
  psc_sort_run(ppsc->sort, ppsc->particles);
  return _case;
}

// ======================================================================
// psc_calc_div_j

static void
do_calc_div_j(struct psc *psc, int p, fields_c_t *flds, fields_c_t *div_j)
{
  creal h[3];
  for (int d = 0; d < 3; d++) {
    if (psc->domain.gdims[d] == 1) {
      h[d] = 0.;
    } else {
      h[d] = 1. / psc->dx[d];
    }
  }

  psc_foreach_3d_g(psc, p, jx, jy, jz) {
    F3_C(div_j,0, jx,jy,jz) =
      (F3_C(flds,JXI, jx,jy,jz) - F3_C(flds,JXI, jx-1,jy,jz)) * h[0] +
      (F3_C(flds,JYI, jx,jy,jz) - F3_C(flds,JYI, jx,jy-1,jz)) * h[1] +
      (F3_C(flds,JZI, jx,jy,jz) - F3_C(flds,JZI, jx,jy,jz-1)) * h[2];
  } psc_foreach_3d_g_end;
}

static void
psc_calc_div_j(struct psc *psc, mfields_base_t *flds_base, mfields_c_t *div_j)
{
  mfields_c_t *flds = psc_mfields_get_c(flds_base, JXI, JXI + 3);

  psc_foreach_patch(psc, p) {
    do_calc_div_j(psc, p, psc_mfields_get_patch_c(flds, p),
		  psc_mfields_get_patch_c(div_j, p));
  }

  psc_mfields_put_c(flds, flds_base, 0, 0);

  psc_bnd_add_ghosts(psc->bnd, div_j, 0, 1);
}

// ======================================================================
// psc_check_continuity

static mfields_c_t *
fld_create(struct psc *psc, int nr_fields)
{
  mfields_c_t *fld = psc_mfields_create(psc_comm(psc));
  psc_mfields_set_type(fld, "c");
  psc_mfields_set_domain(fld, psc->mrc_domain);
  psc_mfields_set_param_int3(fld, "ibn", psc->ibn);
  psc_mfields_set_param_int(fld, "nr_fields", nr_fields);
  psc_mfields_setup(fld);

  return fld;
}

static void
psc_calc_rho(struct psc *psc, mparticles_base_t *particles, mfields_c_t *rho)
{
  mfields_c_t *dens = fld_create(psc, 3);

  psc_moments_calc_densities(psc->moments, NULL, particles, dens);
  // rho = NE + NI
  psc_mfields_copy_comp(rho, 0, dens, 0); // FIXME, waxpy would be nicer
  psc_mfields_axpy_comp(rho, 0, 1., dens, 1);

  psc_mfields_destroy(dens);
}

static void
psc_check_continuity(struct psc *psc, mparticles_base_t *particles,
		     mfields_base_t *flds, double eps)
{
  mfields_c_t *rho_m = fld_create(psc, 1);
  mfields_c_t *rho_p = fld_create(psc, 1);
  mfields_c_t *div_j = fld_create(psc, 1);

  psc_shift_particle_positions(psc, particles, -.5 * psc->dt);
  psc_calc_rho(psc, particles, rho_m);
  psc_shift_particle_positions(psc, particles,  1. * psc->dt);
  psc_calc_rho(psc, particles, rho_p);
  psc_shift_particle_positions(psc, particles, -.5 * psc->dt);

  // rho_p = (rho_p - rho_m) / dt
  psc_mfields_axpy(rho_p, -1., rho_m);
  psc_mfields_scale(rho_p, 1. / psc->dt);

  psc_calc_div_j(psc, flds, div_j);

  //  psc_dump_field(div_j, 0, "div_j");
  //  psc_dump_field(rho_p, 0, "dt_rho");

  double max_err = 0.;
  psc_foreach_patch(psc, p) {
    fields_c_t *p_rho_p = psc_mfields_get_patch_c(rho_p, p);
    fields_c_t *p_div_j = psc_mfields_get_patch_c(div_j, p);
    psc_foreach_3d(psc, p, jx, jy, jz, 0, 0) {
      creal dt_rho = F3_C(p_rho_p,0, jx,jy,jz);
      creal div_j = F3_C(p_div_j,0, jx,jy,jz);
      max_err = fmax(max_err, fabs(dt_rho + div_j));
      if (fabs(dt_rho + div_j) > eps) {
	printf("(%d,%d,%d): %g -- %g diff %g\n", jx, jy, jz,
	       dt_rho, -div_j, dt_rho + div_j);
      }
    } psc_foreach_3d_end;
  }
  printf("continuity: max_err = %g (thres %g)\n", max_err, eps);

  //  psc_mfields_axpy(rho_p, +1., div_j);
  //  psc_dump_field(rho_p, 0, "cont_diff");

  assert(max_err <= eps);

  psc_mfields_destroy(rho_m);
  psc_mfields_destroy(rho_p);
  psc_mfields_destroy(div_j);
}

static void
run_test(bool is_ref, const char *s_push_particles, double eps_particles, double eps_fields,
	 struct psc_case *(*create_test)(const char *), const char *push)
{
  printf("=== testing push_part_yz%s() %s %s\n", push, s_push_particles,
	 is_ref ? "(ref)" : "");

  struct psc_case *_case = create_test(s_push_particles);
  dump(s_push_particles, 0);
  if (strlen(push) == 0) {
    psc_push_particles_run(ppsc->push_particles, ppsc->particles, ppsc->flds);
  } else if (strcmp(push, "_a") == 0) {
    psc_push_particles_push_yz_a(ppsc->push_particles, ppsc->particles, ppsc->flds);
  } else if (strcmp(push, "_b") == 0) {
    psc_push_particles_push_yz_b(ppsc->push_particles, ppsc->particles, ppsc->flds);
  }
  psc_bnd_exchange_particles(ppsc->bnd, ppsc->particles);
  psc_sort_run(ppsc->sort, ppsc->particles);
  dump(s_push_particles, 1);
  if (strlen(push) == 0) { // only check continuity for full pusher
    psc_check_continuity(ppsc, ppsc->particles, ppsc->flds, 1e-14);
  }
  if (is_ref) {
    psc_save_particles_ref(ppsc, ppsc->particles);
    psc_save_fields_ref(ppsc, ppsc->flds);
  } else {
    psc_check_particles_ref(ppsc, ppsc->particles, eps_particles, "push_part_yz()");
    if (check_currents && strlen(push) == 0) { // only check currents for full pusher
      psc_check_currents_ref(ppsc, ppsc->flds, eps_fields);
    }
  }
  psc_case_destroy(_case);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  mrc_params_get_option_bool("dump", &do_dump);
  mrc_params_get_option_bool("check_currents", &check_currents);

  // ----------------------------------------------------------------------
  // push_yz_a

  run_test(true, "fortran", 0., 0., create_test, "_a");
  run_test(false, "generic_c", 1e-7, 1e-7, create_test, "_a");
#ifdef xUSE_CUDA
  run_test(false, "cuda", 1e-3, 1e-2, create_test, "_a");
#endif
#ifdef USE_SSE2
  run_test(false, "sse2", 1e-7, 2e-6, create_test, "_a");
#endif

  // ----------------------------------------------------------------------
  // push_yz_b

  run_test(true, "fortran", 0., 0., create_test, "_b");
  run_test(false, "generic_c", 1e-7, 1e-7, create_test, "_b");
#ifdef xUSE_CUDA
  run_test(false, "cuda", 2e-3, 1e-2, create_test, "_b");
#endif
#ifdef USE_SSE2
  run_test(false, "sse2", 1e-7, 2e-6, create_test, "_b");
#endif

  // ----------------------------------------------------------------------
  // push_yz

  run_test(true, "fortran", 0., 0., create_test, "");
  run_test(false, "generic_c", 1e-7, 1e-7, create_test, "");
#ifdef xUSE_CUDA
  run_test(false, "cuda", 2e-3, 1e-3, create_test, "");
#endif
#ifdef USE_SSE2
  run_test(false, "sse2", 1e-7, 2e-6, create_test, "");
#endif

  prof_print();

  MPI_Finalize();
}
