
#include <stdio.h>

#include "psc.h"
#include "profile/profile.h"

int
main()
{
  psc_create_test_3("fortran");
  //  psc_dump_particles("part-0.asc");
  psc_push_part_yz();
  //  psc_dump_particles("part-1.asc");
  psc_save_particles_ref();
  psc_save_fields_ref();
  psc_destroy();

  //  psc_create_test_1("generic_c");
  //  psc_push_part_yz_a();
  //  psc_dump_particles("part-2.asc");
  // psc_check_particles_ref();
  //  psc_destroy();
  
  psc_create_test_3("sse2");
  psc_push_part_yz();
  //  psc_dump_particles("part-3.asc");
  //These checks will fail for single precision
  //psc_check_particles_ref();
  //psc_check_currents_ref();
  psc_destroy();

  prof_print();
}
