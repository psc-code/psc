
#include <stdio.h>

#include "psc.h"

int
main()
{
  psc_create_test_1("fortran");
  //  psc_dump_particles("part-0.asc");
  PIC_push_part_yz_a();
  //  psc_dump_particles("part-1.asc");
  psc_save_particles_ref();
  psc_free();

  psc_create_test_1("generic_c");
  psc_push_part_yz_a();
  //  psc_dump_particles("part-2.asc");
  psc_check_particles_ref();
  psc_free();
}
