
#include "psc.h"
#include "util/profile.h"
#include "util/params.h"

#include <stdio.h>
#include <mpi.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  params_init(argc, argv);

  struct psc_mod_config conf_fortran = {
    .mod_particle = "fortran",
  };
  struct psc_mod_config conf_generic_c = {
    .mod_particle = "generic_c",
  };

  psc_create_test_xz(&conf_fortran);
  //  psc_dump_particles("part-0.asc");
  psc_push_particles();
  //  psc_dump_particles("part-1.asc");
  psc_save_particles_ref();
  psc_save_fields_ref();
  psc_destroy();

  psc_create_test_xz(&conf_generic_c);
  psc_push_particles();
  //  psc_dump_particles("part-2.asc");
  psc_check_particles_ref(1e-7);
  psc_check_currents_ref(1e-7);
  psc_destroy();

  prof_print();

  MPI_Finalize();
}
