
#include "psc_testing.h"
#include "psc_push_particles.h"
#include <mrc_profile.h>

#include <stdio.h>
#include <mpi.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  struct psc_mod_config conf_fortran = {
    .mod_particle = "fortran",
    .mod_sort = "qsort",
  };
  struct psc_mod_config conf_generic_c = {
    .mod_particle = "generic_c",
    .mod_sort = "qsort",
  };

#ifdef USE_SSE2
  struct psc_mod_config conf_sse2 = {
    .mod_particle = "sse2",
    .mod_sort = "qsort",
  };
#endif

  printf("=== testing push_part_yz_a()\n");

  psc_create_test_yz(&conf_fortran);
  mparticles_base_t *particles = &psc.particles;
  mfields_base_t *flds = &psc.flds;
  //  psc_dump_particles("part-0");
  psc_push_particles_push_yz_a(psc.push_particles, particles, flds);
  //  psc_dump_particles("part-1");
  psc_save_particles_ref(&psc, particles);
  psc_destroy(&psc);

  psc_create_test_yz(&conf_generic_c);
  psc_push_particles_push_yz_a(psc.push_particles, particles, flds);
  //  psc_dump_particles("part-2");
  psc_check_particles_ref(&psc, particles, 1e-6, "push_part_yz_a -- generic_c");
  psc_destroy(&psc);


#ifdef USE_CUDA
  struct psc_mod_config conf_cuda = {
    .mod_particle = "cuda",
    .mod_sort = "qsort",
  };
  psc_create_test_yz(&conf_cuda);
  psc_push_particles_push_yz_a(psc.push_particles, particles, flds);
  psc_check_particles_ref(1e-6, "push_part_yz_a -- cuda");
  psc_destroy();
#endif 

#ifdef USE_SSE2
  psc_create_test_yz(&conf_sse2);
  psc_push_particles_push_yz_a(psc.push_particles, particles, flds);
  //  psc_dump_particles("part-3");
  psc_check_particles_ref(&psc, particles, 1e-6, "push_part_yz_a -- sse2");
  psc_destroy(&psc);
#endif

  printf("=== testing push_part_yz_b()\n");

  psc_create_test_yz(&conf_fortran);
  //  psc_dump_particles("part-0");
  psc_push_particles_push_yz_b(psc.push_particles, particles, flds);
  //  psc_dump_particles("part-1");
  psc_save_particles_ref(&psc, particles);
  psc_destroy(&psc);

  psc_create_test_yz(&conf_generic_c);
  psc_push_particles_push_yz_b(psc.push_particles, particles, flds);
  //  psc_dump_particles("part-2");
  psc_check_particles_ref(&psc, particles, 1e-6, "push_part_yz_b -- generic_c");
  psc_destroy(&psc);

#ifdef USE_CUDA
  psc_create_test_yz(&conf_cuda);
  psc_push_particles_push_yz_b(psc.push_particles, particles, flds);
  //  psc_dump_particles("part-3");
  psc_check_particles_ref(1e-3, "push_part_yz_b -- cuda");
  psc_destroy();
#endif

#ifdef USE_SSE2
  psc_create_test_yz(&conf_sse2);
  psc_push_particles_push_yz_b(psc.push_particles, particles, flds);
  //  psc_dump_particles("part-3");
  psc_check_particles_ref(&psc, particles, 1e-6, "push_part_yz_b -- sse2");
  psc_destroy(&psc);
#endif

  printf("=== testing push_part_yz()\n");

  psc_create_test_yz(&conf_fortran);
  //  psc_dump_particles("part-0");
  psc_push_particles_run(psc.push_particles, particles, flds);
  //  psc_dump_particles("part-1");
  psc_save_particles_ref(&psc, particles);
  psc_save_fields_ref(&psc, flds);
  psc_destroy(&psc);

  psc_create_test_yz(&conf_generic_c);
  //  psc_dump_particles("part-0");
  psc_push_particles_run(psc.push_particles, particles, flds);
  //  psc_dump_particles("part-1");
  psc_check_particles_ref(&psc, particles, 1e-7, "push_part_yz -- generic_c");
  psc_check_currents_ref(&psc, flds, 1e-7);
  psc_destroy(&psc);

#ifdef USE_SSE2
  psc_create_test_yz(&conf_sse2);
  psc_push_particles_run(psc.push_particles, particles, flds);
  //  psc_dump_particles("part-3");
  psc_check_particles_ref(&psc, particles, 1e-8, "push_part_yz -- sse2");
  psc_check_currents_ref(&psc, flds, 2e-6); 
  psc_destroy(&psc);
#endif


  prof_print();

  MPI_Finalize();
}
