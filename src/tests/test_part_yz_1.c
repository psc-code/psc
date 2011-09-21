
#include "psc_testing.h"
#include "psc_push_particles.h"
#include "psc_sort.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

static bool do_dump = false;
static bool check_currents = true;

struct psc_test {
};

#define to_psc_test(psc) mrc_to_subobj(psc, struct psc_test)

// ----------------------------------------------------------------------
// psc_test_create

static void
psc_test_create(struct psc *psc)
{
  // new defaults (dimensionless) for this case
  psc->prm.qq = 1.;
  psc->prm.mm = 1.;
  psc->prm.tt = 1.;
  psc->prm.cc = 1.;
  psc->prm.eps0 = 1.;

  psc->prm.lw = 2.*M_PI;
  psc->prm.i0 = 0.;
  psc->prm.n0 = 1.;
  psc->prm.e0 = 1.;

  psc->prm.nicell = 100;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 16;
  psc->domain.gdims[2] = 16;

  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_part[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part[2] = BND_PART_PERIODIC;

  psc_sort_set_type(psc->sort, "countsort2");
}

// ----------------------------------------------------------------------
// psc_test_set_from_options

static void
psc_test_set_from_options(struct psc *psc)
{
  psc->domain.length[0] = 1.;
  psc->domain.length[1] = 1.;
  psc->domain.length[2] = 1.;
}

// ----------------------------------------------------------------------
// psc_test_init_field

static double
psc_test_init_field(struct psc *psc, double x[3], int m)
{
  switch (m) {
  case EX: return x[1] + x[2];
  case EY: return x[1] + x[2];
  case EZ: return x[1] + x[2];
  case HX: return x[1] + x[2];
  case HY: return x[1] + x[2];
  case HZ: return x[1] + x[2];
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_test_init_npt

static void
psc_test_init_npt(struct psc *psc, int kind, double x[3],
		    struct psc_particle_npt *npt)
{
  npt->n = 1.;
  switch (kind) {
  case 0: // electrons
    npt->q = -1.;
    npt->m = 1.;
    break;
  case 1: // ions
    npt->q = 1.;
    npt->m = 100;
    break;
  default:
    assert(0);
  }
}

// ======================================================================
// psc_test_ops

struct psc_ops psc_test_ops = {
  .name             = "test",
  .size             = sizeof(struct psc_test),
  .create           = psc_test_create,
  .set_from_options = psc_test_set_from_options,
  .init_field       = psc_test_init_field,
  .init_npt         = psc_test_init_npt,
};

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

static struct psc *
create_test(const char *s_push_particles)
{
  struct psc *psc = psc_create(MPI_COMM_WORLD);
  psc_set_from_options(psc);
  psc_push_particles_set_type(psc->push_particles, s_push_particles);
  psc_sort_set_type(psc->sort, "countsort2");
  psc_setup(psc);
  psc_sort_run(psc->sort, psc->particles);
  return psc;
}

static void
run_test(bool is_ref, const char *s_push_particles, double eps_particles, double eps_fields,
	 struct psc *(*create_test)(const char *), const char *push)
{
  printf("=== testing push_part_yz%s() %s %s\n", push, s_push_particles,
	 is_ref ? "(ref)" : "");

  struct psc *psc = create_test(s_push_particles);
  dump(s_push_particles, 0);
  if (strlen(push) == 0) {
    psc_push_particles_run(psc->push_particles, psc->particles, psc->flds);
  } else if (strcmp(push, "_a") == 0) {
    psc_push_particles_push_yz_a(psc->push_particles, psc->particles, psc->flds);
  } else if (strcmp(push, "_b") == 0) {
    psc_push_particles_push_yz_b(psc->push_particles, psc->particles, psc->flds);
  }
  psc_sort_run(psc->sort, psc->particles);
  dump(s_push_particles, 1);
  if (is_ref) {
    psc_save_particles_ref(psc, psc->particles);
    psc_save_fields_ref(psc, psc->flds);
  } else {
    psc_check_particles_ref(psc, psc->particles, eps_particles, "push_part_yz()");
    if (check_currents && strlen(push) == 0) { // only check currents for full pusher
      psc_check_currents_ref(psc, psc->flds, eps_fields);
    }
  }
  psc_destroy(psc);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  mrc_class_register_subclass(&mrc_class_psc, &psc_test_ops);

  mrc_params_get_option_bool("dump", &do_dump);
  mrc_params_get_option_bool("check_currents", &check_currents);

  // ----------------------------------------------------------------------
  // push_yz 1st order

  run_test(true, "1st", 0., 0., create_test, "");
  // since the fields are linear functions of position, 1st order / 2nd order
  // field interpolation should give the same result
  run_test(false, "generic_c", 1e-7, 1e-0, create_test, "");
#ifdef USE_CUDA
  run_test(false, "cuda_1st", 1e-6, 1e-0, create_test, "");
#endif

  prof_print();

  MPI_Finalize();
}
