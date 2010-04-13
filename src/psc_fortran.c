
#include "psc.h"

static void
fortran_particles_from_fortran()
{
}

static void
fortran_particles_to_fortran()
{
}

static void
fortran_push_part_yz()
{
  PIC_push_part_yz();
}

static void
fortran_push_part_z()
{
  PIC_push_part_z();
}

static void
fortran_push_part_yz_a()
{
  PIC_push_part_yz_a();
}

struct psc_ops psc_ops_fortran = {
  .name = "fortran",
  .particles_from_fortran = fortran_particles_from_fortran,
  .particles_to_fortran   = fortran_particles_to_fortran,
  .push_part_yz           = fortran_push_part_yz,
  .push_part_z            = fortran_push_part_z,
  .push_part_yz_a         = fortran_push_part_yz_a,
};
