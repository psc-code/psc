#include "psc_testing.h"
#include <mrc_profile.h>

#include <stdio.h>
#include <mpi.h>

int
main(int argc, char **argv)
{
 
  MPI_Init(&argc, &argv);
  
  struct psc_mod_config conf_cbe = {
    .mod_particle = "cbe",
    .mod_sort = "countsort2",
  };
  psc_create(&conf_cbe);
  //PIC_find_cell_indices();
  //psc_sort();
  psc_push_part_yz();
  //  psc_dump_particles("part-3");
  //psc_check_particles_ref(1e-7);
  psc_destroy();
  //  prof_print();

  MPI_Finalize();

}
