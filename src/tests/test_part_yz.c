
#include "psc_testing.h"
#include "psc_push_particles.h"
#include "psc_sort.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <stdio.h>
#include <mpi.h>

static bool do_dump = false;

static void
dump(const char *basename, int cnt)
{
  if (!do_dump)
    return;

  char s[200];
  sprintf(s, "part_%s_%d", basename, cnt);
  psc_dump_particles(ppsc->particles, s);
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
	 struct psc_case *(*create_test)(const char *s_push_particles))
{
  printf("=== testing push_part_yz() %s %s\n", s_push_particles, is_ref ? "(ref)" : "");

  struct psc_case *_case = create_test(s_push_particles);
  dump(s_push_particles, 0);
  psc_push_particles_run(ppsc->push_particles, ppsc->particles, ppsc->flds);
  psc_sort_run(ppsc->sort, ppsc->particles);
  dump(s_push_particles, 1);
  if (is_ref) {
    psc_save_particles_ref(ppsc, ppsc->particles);
    psc_save_fields_ref(ppsc, ppsc->flds);
  } else {
    psc_check_particles_ref(ppsc, ppsc->particles, eps_particles, "push_part_z -- generic_c");
    psc_check_currents_ref(ppsc, ppsc->flds, eps_fields);
  }
  psc_case_destroy(_case);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  mrc_params_get_option_bool("dump", &do_dump);

  printf("=== testing push_part_yz_a() ref: fortran\n");

  struct psc_case *_case = psc_create_test_yz();
  psc_push_particles_set_type(ppsc->push_particles, "fortran");
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  mparticles_base_t *particles = ppsc->particles;
  psc_sort_run(ppsc->sort, particles);
  psc_push_particles_push_yz_a(ppsc->push_particles, particles, ppsc->flds);
  psc_save_particles_ref(ppsc, particles);
  psc_case_destroy(_case);

  printf("=== testing push_part_yz_a() generic_c\n");
  _case = psc_create_test_yz();
  psc_push_particles_set_type(ppsc->push_particles, "generic_c");
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  particles = ppsc->particles;
  psc_sort_run(ppsc->sort, particles);
  psc_push_particles_push_yz_a(ppsc->push_particles, particles, ppsc->flds);
  psc_check_particles_ref(ppsc, particles, 1e-6, "push_part_yz_a -- generic_c");
  psc_case_destroy(_case);

#ifdef USE_CUDA
  printf("=== testing push_part_yz_a() cuda\n");
  _case = psc_create_test_yz();
  psc_push_particles_set_type(ppsc->push_particles, "cuda");
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  particles = ppsc->particles;
  psc_sort_run(ppsc->sort, particles);
  psc_push_particles_push_yz_a(ppsc->push_particles, particles, ppsc->flds);
  psc_check_particles_ref(ppsc, particles, 1e-6, "push_part_yz_a -- cuda");
  psc_case_destroy(_case);
#endif 

#ifdef USE_SSE2
  printf("=== testing push_part_yz_a() sse2\n");
  _case = psc_create_test_yz();
  psc_push_particles_set_type(ppsc->push_particles, "sse2");
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  particles = &ppsc->particles;
  psc_sort_run(ppsc->sort, particles);
  psc_push_particles_push_yz_a(ppsc->push_particles, particles, ppsc->flds);
  psc_check_particles_ref(ppsc, particles, 1e-6, "push_part_yz_a -- sse2");
  psc_case_destroy(_case);
#endif

  // ----------------------------------------------------------------------
  printf("=== testing push_part_yz_b() ref: fortran\n");

  _case = psc_create_test_yz();
  psc_push_particles_set_type(ppsc->push_particles, "fortran");
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  particles = ppsc->particles;
  psc_sort_run(ppsc->sort, particles);
  psc_push_particles_push_yz_b(ppsc->push_particles, particles, ppsc->flds);
  psc_save_particles_ref(ppsc, particles);
  psc_case_destroy(_case);

  printf("=== testing push_part_yz_b() generic_c\n");
  _case = psc_create_test_yz();
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_push_particles_set_type(ppsc->push_particles, "generic_c");
  psc_case_setup(_case);
  particles = ppsc->particles;
  psc_sort_run(ppsc->sort, particles);
  psc_push_particles_push_yz_b(ppsc->push_particles, particles, ppsc->flds);
  psc_check_particles_ref(ppsc, particles, 1e-6, "push_part_yz_b -- generic_c");
  psc_case_destroy(_case);

#ifdef USE_CUDA
  printf("=== testing push_part_yz_b() cuda\n");
  _case = psc_create_test_yz();
  psc_push_particles_set_type(ppsc->push_particles, "cuda");
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  particles = ppsc->particles;
  psc_sort_run(ppsc->sort, particles);
  psc_push_particles_push_yz_b(ppsc->push_particles, particles, ppsc->flds);
  psc_check_particles_ref(ppsc, particles, 1e-3, "push_part_yz_b -- cuda");
  psc_case_destroy(_case);
#endif

#ifdef USE_SSE2
  printf("=== testing push_part_yz_b() sse2\n");
  _case = psc_create_test_yz();
  psc_push_particles_set_type(ppsc->push_particles, "sse2");
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  particles = &ppsc->particles;
  psc_sort_run(ppsc->sort, particles);
  psc_push_particles_push_yz_b(ppsc->push_particles, particles, ppsc->flds);
  psc_check_particles_ref(ppsc, particles, 1e-6, "push_part_yz_b -- sse2");
  psc_case_destroy(_case);
#endif

  // ----------------------------------------------------------------------

  run_test(true, "fortran", 0., 0., create_test);
  run_test(false, "generic_c", 1e-7, 1e-7, create_test);

#ifdef USE_CUDA
  run_test(false, "cuda", 1e-3, 1e-2, create_test);
#endif

#ifdef USE_SSE2
  run_test(false, "sse2", 1e-7, 2e-6, create_test);
#endif

  prof_print();

  MPI_Finalize();
}
