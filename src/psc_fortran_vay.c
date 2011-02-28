
#include "psc.h"
#include "util/profile.h"

static void
fortran_vay_push_part_xy()
{
  particles_fortran_t pp;
  particles_fortran_get(&pp);
  fields_fortran_t pf;
  fields_fortran_get(&pf, EX, EX + 6);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_vay_part_xy", 1., 0, psc.pp.n_part * 11 * sizeof(double));
  }
  prof_start(pr);
  PIC_push_part_xy(&pp, &pf);            // here we can choose another particle pusher
  prof_stop(pr);

  particles_fortran_put(&pp);
  fields_fortran_put(&pf, JXI, JXI + 3);
}

static void
fortran_vay_push_part_xz()
{
  particles_fortran_t pp;
  particles_fortran_get(&pp);
  fields_fortran_t pf;
  fields_fortran_get(&pf, EX, EX + 6);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_vay_part_xz", 1., 0, psc.pp.n_part * 11 * sizeof(double));
  }
  prof_start(pr);
  PIC_push_part_xz(&pp, &pf);
  prof_stop(pr);

  particles_fortran_put(&pp);
  fields_fortran_put(&pf, JXI, JXI + 3);
}

static void
fortran_vay_push_part_yz()
{
  particles_fortran_t pp;
  particles_fortran_get(&pp);
  fields_fortran_t pf;
  fields_fortran_get(&pf, EX, EX + 6);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_vay_part_yz", 1., 0, psc.pp.n_part * 11 * sizeof(double));
  }
  prof_start(pr);
  PIC_push_part_yz(&pp, &pf);
  prof_stop(pr);

  particles_fortran_put(&pp);
  fields_fortran_put(&pf, JXI, JXI + 3);
}

static void
fortran_vay_push_part_xyz()
{
  particles_fortran_t pp;
  particles_fortran_get(&pp);
  fields_fortran_t pf;
  fields_fortran_get(&pf, EX, EX + 6);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_vay_part_xyz", 1., 0, psc.pp.n_part * 11 * sizeof(double));
  }
  prof_start(pr);
  PIC_push_part_xyz(&pp, &pf);
  prof_stop(pr);

  particles_fortran_put(&pp);
  fields_fortran_put(&pf, JXI, JXI + 3);
}

static void
fortran_vay_push_part_z()
{
  particles_fortran_t pp;
  particles_fortran_get(&pp);
  fields_fortran_t pf;
  fields_fortran_get(&pf, EX, EX + 6);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_vay_part_z", 1., 0, psc.pp.n_part * 11 * sizeof(double));
  }
  prof_start(pr);
  PIC_push_part_z_vay(&pp, &pf);
  prof_stop(pr);

  particles_fortran_put(&pp);
  fields_fortran_put(&pf, JXI, JXI + 3);
}

static void
fortran_vay_push_part_yz_a()
{
  particles_fortran_t pp;
  particles_fortran_get(&pp);
  fields_fortran_t pf;
  fields_fortran_get(&pf, EX, EX + 6);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_vay_part_yz_a", 1., 0, psc.pp.n_part * 11 * sizeof(double));
  }
  prof_start(pr);
  PIC_push_part_yz_a(&pp, &pf);
  prof_stop(pr);

  particles_fortran_put(&pp);
  fields_fortran_put(&pf, JXI, JXI + 3);
}

static void
fortran_vay_push_part_yz_b()
{
  particles_fortran_t pp;
  particles_fortran_get(&pp);
  fields_fortran_t pf;
  fields_fortran_get(&pf, EX, EX + 6);
  
  static int pr;
  if (!pr) {
    pr = prof_register("fort_vay_part_yz_b", 1., 0, psc.pp.n_part * 11 * sizeof(double));
  }
  prof_start(pr);
  PIC_push_part_yz_b(&pp, &pf);
  prof_stop(pr);

  particles_fortran_put(&pp);
  fields_fortran_put(&pf, JXI, JXI + 3);
}

struct psc_ops psc_ops_fortran_vay = {
  .name = "fortran_vay",
  .push_part_xy           = fortran_vay_push_part_xy,
  .push_part_xz           = fortran_vay_push_part_xz,
  .push_part_yz           = fortran_vay_push_part_yz,
  .push_part_xyz          = fortran_vay_push_part_xyz,
  .push_part_z            = fortran_vay_push_part_z,
  .push_part_yz_a         = fortran_vay_push_part_yz_a,
  .push_part_yz_b         = fortran_vay_push_part_yz_b,
};

