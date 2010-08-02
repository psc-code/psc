
#include "psc.h"
#include "util/profile.h"

#include <stdio.h>
#include <mpi.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  struct psc_mod_config conf_fortran = {
    .mod_particle = "fortran",
  };
  struct psc_mod_config conf_generic_c = {
    .mod_particle = "generic_c",
  };

  printf("=== testing push_part_yz_a()\n");

  psc_create_test_yz(&conf_fortran);
  //  psc_dump_particles("part-0");
  psc_push_part_yz_a();
  //  psc_dump_particles("part-1");
  psc_save_particles_ref();
  psc_destroy();

  psc_create_test_yz(&conf_generic_c);
  psc_push_part_yz_a();
  //  psc_dump_particles("part-2");
  psc_check_particles_ref(1e-4);
  psc_destroy();

#ifdef USE_CUDA
  struct psc_mod_config conf_cuda = {
    .mod_particle = "cuda",
  };
  psc_create_test_yz(&conf_cuda);
  psc_push_part_yz_a();
  psc_check_particles_ref();
  psc_destroy();
#endif 

#ifdef USE_SSE2
  struct psc_mod_config conf_sse2 = {
    .mod_particle = "sse2",
  };
  psc_create_test_yz(&conf_sse2);
  psc_push_part_yz_a();
  //  psc_dump_particles("part-3");
  psc_check_particles_ref(5e-4);
  psc_destroy();
#endif

#if 1
  printf("=== testing push_part_yz_b()\n");

  psc_create_test_yz(&conf_fortran);
  //  psc_dump_particles("part-0");
  psc_push_part_yz_b();
  //  psc_dump_particles("part-1");
  psc_save_particles_ref();
  psc_destroy();

  psc_create_test_yz(&conf_generic_c);
  psc_push_part_yz_b();
  //  psc_dump_particles("part-2");
  psc_check_particles_ref(1e-2);
  psc_destroy();

#ifdef USE_CUDA
  psc_create_test_yz(&conf_cuda);
  psc_push_part_yz_b();
  //  psc_dump_particles("part-3");
  psc_check_particles_ref();
  psc_destroy();
#endif

#ifdef USE_SSE2
  psc_create_test_yz(&conf_sse2);
  psc_push_part_yz_b();
  //  psc_dump_particles("part-3");
  psc_check_particles_ref(1e-2);
  psc_destroy();
#endif
#endif

  printf("=== testing push_part_yz()\n");

  psc_create_test_yz(&conf_fortran);
  //  psc_dump_particles("part-0");
  psc_push_part_yz();
  //  psc_dump_particles("part-1");
  psc_save_particles_ref();
  psc_save_fields_ref();
  psc_destroy();

#ifdef USE_SSE2
  psc_create_test_yz(&conf_sse2);
  psc_push_part_yz();
  //  psc_dump_particles("part-3");
  psc_check_particles_ref(1e-4);
  psc_check_currents_ref(1e-4); 
  psc_destroy();
#endif


  prof_print();

  MPI_Finalize();
}
