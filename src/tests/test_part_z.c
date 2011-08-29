
#include "psc_testing.h"
#include "psc_push_particles.h"
#include "psc_sort.h"
#include <mrc_profile.h>

#include <stdio.h>
#include <mpi.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  // ----------------------------------------------------------------------
  printf("=== testing push_part_z() ref: fortran\n");

  struct psc_case *_case = psc_create_test_z();
  psc_push_particles_set_type(ppsc->push_particles, "fortran");
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  mparticles_base_t *particles = ppsc->particles;
  psc_sort_run(ppsc->sort, particles);
  //  psc_dump_particles("part-0");
  psc_push_particles_run(ppsc->push_particles, particles, ppsc->flds);
  //  psc_dump_particles("part-1");
  psc_save_particles_ref(ppsc, particles);
  psc_save_fields_ref(ppsc, ppsc->flds);
  psc_case_destroy(_case);

  printf("=== testing push_part_z() generic_c\n");
  _case = psc_create_test_z();
  psc_push_particles_set_type(ppsc->push_particles, "generic_c");
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  particles = ppsc->particles;
  psc_sort_run(ppsc->sort, particles);
  psc_push_particles_run(ppsc->push_particles, particles, ppsc->flds);
  psc_check_particles_ref(ppsc, particles, 1e-7, "push_part_z -- generic_c");
  psc_check_currents_ref(ppsc, ppsc->flds, 1e-7);
  psc_case_destroy(_case);

#ifdef USE_CUDA
  printf("=== testing push_part_z() cuda\n");
  _case = psc_create_test_z();
  psc_push_particles_set_type(ppsc->push_particles, "cuda");
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  particles = ppsc->particles;
  psc_sort_run(ppsc->sort, particles);
  psc_push_particles_run(ppsc->push_particles, particles, ppsc->flds);
  psc_check_particles_ref(ppsc, particles, 1e-3, "push_part_z -- cuda");
  psc_check_currents_ref(ppsc, ppsc->flds, 1e-2); 
  psc_case_destroy(_case);
#endif

#ifdef USE_SSE2
  printf("=== testing push_part_z() sse2\n");
  _case = psc_create_test_z();
  psc_push_particles_set_type(ppsc->push_particles, "sse2");
  psc_sort_set_type(ppsc->sort, "countsort2");
  psc_case_setup(_case);
  particles = ppsc->particles;
  psc_sort_run(ppsc->sort, particles);
  psc_push_particles_run(ppsc->push_particles, particles, ppsc->flds);
  psc_check_particles_ref(ppsc, particles, 1e-8, "push_part_z -- sse2");
  psc_check_currents_ref(ppsc, ppsc->flds, 2e-6); 
  psc_case_destroy(_case);
#endif

  prof_print();

  MPI_Finalize();
}
