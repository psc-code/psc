
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
fortran_push_part_xz()
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_part_xz", 1., 0, psc.n_part * 11 * sizeof(double));
  }
  prof_start(pr);
  PIC_push_part_xz();
  prof_stop(pr);
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
  .push_part_xz           = fortran_push_part_xz,
  .push_part_yz           = fortran_push_part_yz,
  .push_part_z            = fortran_push_part_z,
  .push_part_yz_a         = fortran_push_part_yz_a,
  .push_part_yz_b         = fortran_push_part_yz_b,
};

// ======================================================================
// fortran sort

static void
fortran_sort()
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_sort", 1., 0, 0);
  }
  prof_start(pr);
  PIC_sort_1();
  prof_stop(pr);
}

struct psc_sort_ops psc_sort_ops_fortran = {
  .name = "fortran",
  .sort = fortran_sort,
};

// ======================================================================
// fortran push field

static void
fortran_push_field_a()
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_push_field_a", 1., 0, 0);
  }
  prof_start(pr);
  PIC_msa();
  prof_stop(pr);
}

static void
fortran_push_field_b()
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_push_field_b", 1., 0, 0);
  }
  prof_start(pr);
  PIC_msb();
  prof_stop(pr);
}

struct psc_push_field_ops psc_push_field_ops_fortran = {
  .name         = "fortran",
  .push_field_a = fortran_push_field_a,
  .push_field_b = fortran_push_field_b,
};

// ======================================================================
// fortran output

static void
fortran_out_field()
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_out_field", 1., 0, 0);
  }
  prof_start(pr);
  OUT_field_1();
  prof_stop(pr);
}

struct psc_output_ops psc_output_ops_fortran = {
  .name      = "fortran",
  .out_field = fortran_out_field,
};

// ======================================================================
// fortran bnd

static void
fortran_add_ghosts(int m)
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_add_ghosts", 1., 0, 0);
  }
  prof_start(pr);
  if (psc.domain.ihi[0] - psc.domain.ilo[0] > 1) {
    PIC_fax(m);
  }
  if (psc.domain.ihi[1] - psc.domain.ilo[1] > 1) {
    PIC_fay(m);
  }
  if (psc.domain.ihi[2] - psc.domain.ilo[2] > 1) {
    PIC_faz(m);
  }
  prof_stop(pr);
}

static void
fortran_fill_ghosts(int m)
{
  static int pr;
  if (!pr) {
    pr = prof_register("fort_fill_ghosts", 1., 0, 0);
  }
  prof_start(pr);
  if (psc.domain.ihi[0] - psc.domain.ilo[0] > 1) {
    PIC_fex(m);
  }
  if (psc.domain.ihi[1] - psc.domain.ilo[1] > 1) {
    PIC_fey(m);
  }
  if (psc.domain.ihi[2] - psc.domain.ilo[2] > 1) {
    PIC_fez(m);
  }
  prof_stop(pr);
}

struct psc_bnd_ops psc_bnd_ops_fortran = {
  .name        = "fortran",
  .add_ghosts  = fortran_add_ghosts,
  .fill_ghosts = fortran_fill_ghosts,
};
