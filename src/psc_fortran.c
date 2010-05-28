
#include "psc.h"
#include "util/profile.h"

static void
fortran_particles_from_fortran()
{
}

static void
fortran_particles_to_fortran()
{
}

static void
fortran_fields_from_fortran()
{
}

static void
fortran_fields_to_fortran()
{
}

static void
fortran_push_part_yz()
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_yz", 1., 0, psc.n_part * 11 * sizeof(double));
  }
  prof_start(pr);
 
  PIC_push_part_yz();

  prof_stop(pr);
}

static void
fortran_push_part_z()
{
  PIC_push_part_z();
}

static void
fortran_push_part_yz_a()
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_yz_a", 1., 0, psc.n_part * 11 * sizeof(double));
  }
  prof_start(pr);
  PIC_push_part_yz_a();
  prof_stop(pr);
}

static void
fortran_push_part_yz_b()
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_yz_b", 1., 0, psc.n_part * 11 * sizeof(double));
  }
  prof_start(pr);
  PIC_push_part_yz_b();
  prof_stop(pr);
}

struct psc_ops psc_ops_fortran = {
  .name = "fortran",
  .particles_from_fortran = fortran_particles_from_fortran,
  .particles_to_fortran   = fortran_particles_to_fortran,
  .fields_from_fortran    = fortran_fields_from_fortran,
  .fields_to_fortran      = fortran_fields_to_fortran,
  .push_part_yz           = fortran_push_part_yz,
  .push_part_z            = fortran_push_part_z,
  .push_part_yz_a         = fortran_push_part_yz_a,
  .push_part_yz_b         = fortran_push_part_yz_b,
};
