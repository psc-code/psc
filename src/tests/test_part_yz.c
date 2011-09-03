
#include "psc_testing.h"
#include "psc_push_particles.h"
#include "psc_sort.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <mpi.h>

static bool do_dump = false;
static bool check_currents = true;

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
  psc_case_setup(_case);
  psc_sort_run(ppsc->sort, ppsc->particles);
  return _case;
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
  psc_sort_run(ppsc->sort, ppsc->particles);
  dump(s_push_particles, 1);
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
#ifdef USE_CUDA
  run_test(false, "cuda", 1e-3, 1e-2, create_test, "_a");
#endif
#ifdef USE_SSE2
  run_test(false, "sse2", 1e-7, 2e-6, create_test, "_a");
#endif

  // ----------------------------------------------------------------------
  // push_yz_b

  run_test(true, "fortran", 0., 0., create_test, "_b");
  run_test(false, "generic_c", 1e-7, 1e-7, create_test, "_b");
#ifdef USE_CUDA
  run_test(false, "cuda", 2e-3, 1e-2, create_test, "_b");
#endif
#ifdef USE_SSE2
  run_test(false, "sse2", 1e-7, 2e-6, create_test, "_b");
#endif

  // ----------------------------------------------------------------------
  // push_yz

  run_test(true, "fortran", 0., 0., create_test, "");
  run_test(false, "generic_c", 1e-7, 1e-7, create_test, "");
#ifdef USE_CUDA
  run_test(false, "cuda", 2e-3, 1e-3, create_test, "");
#endif
#ifdef USE_SSE2
  run_test(false, "sse2", 1e-7, 2e-6, create_test, "");
#endif

  prof_print();

  MPI_Finalize();
}
