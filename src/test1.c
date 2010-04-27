
#include "psc.h"
#include "util/profile.h"

int
main()
{
  psc_create_test_1("fortran");
  
  // psc_dump_particles("part-0.asc");
  psc_push_part_yz();
  // psc_dump_particles("part-1.asc");

  psc_destroy();

  prof_print();
}
