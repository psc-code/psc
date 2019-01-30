
#include "psc_testing.h"
#include "psc_push_particles.h"
#include "psc_bnd.h"
#include "psc_sort.h"
#include "psc_particles_as_c.h"

#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <mpi.h>

#if 0
static bool do_dump = false;

static void
add_particle(double xi, double yi, double zi, double pxi, double pyi, double pzi,
	     double qni, double mni)
{
  struct psc_mparticles *mprts = psc_mparticles_get_as(psc->particles, "c", 0);
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, 0);

  int n = prts->n_part++;
  particle_t *part = mprts[0][];
  part->xi = xi;
  part->yi = yi;
  part->zi = zi;
  part->pxi = pxi;
  part->pyi = pyi;
  part->pzi = pzi;
  part->qni = qni;
  part->mni = mni;
  part->wni = ppsc->prm.nicell; // FIXME, better set nicell to 1 or get rid of it altogether

  psc_mparticles_put_as(mprts, psc->particles, 0);
}

static struct psc_case *
create_test_base(const char *s_push_particles)
{
  struct psc_case *_case = psc_create_test_z();
  psc_push_particles_set_type(ppsc->push_particles, s_push_particles);
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);

  ppsc->dt = ppsc->dx[2];

  struct psc_mparticles *mprts = psc_mparticles_get_as(psc->particles, "c", 0);
  struct psc_particles *prts = psc_mparticles_get_patch(mprts, 0);
  prts->n_part = 0;
  psc_mparticles_put_as(mprts, psc->particles, 0);
  return _case;
}

static struct psc_case *
create_test_1(const char *s_push_particles)
{
  struct psc_case *_case = psc_create_test_z();
  psc_push_particles_set_type(ppsc->push_particles, s_push_particles);
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  psc_sort_run(ppsc->sort, ppsc->particles);
  return _case;
}

static struct psc_case *
create_test_2(const char *s_push_particles)
{
  struct psc_case *_case = create_test_base(s_push_particles);

  double *dx = ppsc->dx;
  add_particle(0.,0.,2.*dx[2], 0.,0.,1000., 1., 1.);

  psc_sort_run(ppsc->sort, ppsc->particles);
  return _case;
}

static struct psc_case *
create_test_3(const char *s_push_particles)
{
  struct psc_case *_case = create_test_base(s_push_particles);

  double *dx = ppsc->dx;
  add_particle(0.,0.,2.*dx[2], 0.,0.,1000., 1., 1.);
  add_particle(0.,0.,3.*dx[2], 0.,0.,1000., 1., 1.);

  psc_sort_run(ppsc->sort, ppsc->particles);
  return _case;
}

static struct psc_case *
create_test_4(const char *s_push_particles)
{
  struct psc_case *_case = create_test_base(s_push_particles);

  double *dx = ppsc->dx;
  add_particle(0.,0.,.4999*dx[2], 0.,0.,1000., 1., 1.);

  psc_sort_run(ppsc->sort, ppsc->particles);
  return _case;
}

static struct psc_case *
create_test_5(const char *s_push_particles)
{
  struct psc_case *_case = create_test_base(s_push_particles);

  add_particle(0.,0.,0., .0,.0, .1, 1., 1.);
  add_particle(0.,0.,0., .0,.0, .1, 1., 1.);
  add_particle(0.,0.,0., .0,.0, .1, 1., 1.);
  add_particle(0.,0.,0., .0,.0, .1, 1., 1.);

  psc_sort_run(ppsc->sort, ppsc->particles);
  return _case;
}

static struct psc_case *
create_test_6(const char *s_push_particles)
{
  struct psc_case *_case = create_test_base(s_push_particles);

  add_particle(0.,0.,1e-6, -.43, .15,  .21, -1., 1.);
  //add_particle(0.,0.,1e-6, 0.03,-.22, -.23, -1., 1.);

  psc_sort_run(ppsc->sort, ppsc->particles);
  return _case;
}

static void
run_test(bool is_ref, const char *s_push_particles, double eps_particles, double eps_fields,
	 struct psc_case *(*create_test)(const char *s_push_particles))
{
  printf("=== testing push_part_z() %s %s\n", s_push_particles, is_ref ? "(ref)" : "");

  struct psc_case *_case = create_test(s_push_particles);
  dump(s_push_particles, 0);
  psc_push_particles_run(ppsc->push_particles, ppsc->particles, ppsc->flds);
  psc_bnd_exchange_particles(ppsc->bnd, ppsc->particles);
  psc_sort_run(ppsc->sort, ppsc->particles);
  dump(s_push_particles, 1);
  if (is_ref) {
    psc_save_particles_ref(ppsc, ppsc->particles);
    psc_save_fields_ref(ppsc, ppsc->flds);
  } else {
    psc_check_particles_ref(ppsc, ppsc->particles, eps_particles, "push_part_z()");
    psc_check_currents_ref(ppsc, ppsc->flds, eps_fields, 3);
  }
  psc_case_destroy(_case);
}
#endif

int
main(int argc, char **argv)
{
#if 0
  psc_testing_init(&argc, &argv);

  int testcase = 1;
  mrc_params_get_option_int("case", &testcase);
  mrc_params_get_option_bool("dump", &do_dump);

  struct psc_case *(*create_test)(const char *);
  switch (testcase) {
  case 1: create_test = create_test_1; break;
  case 2: create_test = create_test_2; break;
  case 3: create_test = create_test_3; break;
  case 4: create_test = create_test_4; break;
  case 5: create_test = create_test_5; break;
  case 6: create_test = create_test_6; break;
  default: assert(0);
  }

  run_test(true, "fortran", 0., 0., create_test);
  run_test(false, "generic_c", 1e-7, 1e-7, create_test);

#ifdef xUSE_CUDA // FIXME
  run_test(false, "cuda", 1e-4, 1e-4, create_test);
#endif

#ifdef USE_SSE2
  run_test(false, "sse2", 1e-7, 2e-6, create_test);
#endif

  psc_testing_finalize();
#endif
}
