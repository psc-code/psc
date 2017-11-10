
#include "psc_push_particles_private.h"

#include "psc_push_fields_private.h"
#include "psc_bnd_fields.h"
#include "psc_bnd.h"

#include <mrc_profile.h>

// ======================================================================
// forward to subclass

extern int pr_time_step_no_comm;
extern double *psc_balance_comp_time_by_patch;

void
psc_push_particles_run(struct psc_push_particles *push,
		       struct psc_mparticles *mprts_base, struct psc_mfields *mflds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("push_particles_run", 1., 0, 0);
  }  

  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);

  struct psc_mparticles *mprts;
  if (ops->particles_type) {
    mprts = psc_mparticles_get_as(mprts_base, ops->particles_type, 0);
  } else {
    mprts = mprts_base;
  }
  struct psc_mfields *mflds;
  if (ops->fields_type) {
    mflds = psc_mfields_get_as(mflds_base, ops->fields_type, EX, EX + 6);
  } else {
    mflds = mflds_base;
  }
  
  prof_start(pr);
  prof_restart(pr_time_step_no_comm);
  psc_stats_start(st_time_particle);

  if (ops->push_mprts) {
    ops->push_mprts(push, mprts, mflds);
  } else {

  int *im = ppsc->domain.gdims;

  if (im[0] > 1 && im[1] == 1 && im[2] == 1) { // x
    assert(ops->push_mprts_x);
    ops->push_mprts_x(push, mprts, mflds);
  } else if (im[0] == 1 && im[1] > 1 && im[2] == 1) { // y
    assert(ops->push_mprts_y);
    ops->push_mprts_y(push, mprts, mflds);
  } else if (im[0] == 1 && im[1] == 1 && im[2] > 1) { // y
    assert(ops->push_mprts_z);
    ops->push_mprts_z(push, mprts, mflds);
  } else if (im[0] > 1 && im[1] > 1 && im[2] == 1) { // xy
    assert(ops->push_mprts_xy);
    ops->push_mprts_xy(push, mprts, mflds);
  } else if (im[0] > 1 && im[1] == 1 && im[2] > 1) { // xz
    assert(ops->push_mprts_xz);
    ops->push_mprts_xz(push, mprts, mflds);
  } else if (im[0] == 1 && im[1] > 1 && im[2] > 1) { // yz
    assert(ops->push_mprts_yz);
    ops->push_mprts_yz(push, mprts, mflds);
  } else if (im[0] > 1 && im[1] > 1 && im[2] > 1) { // xyz
    assert(ops->push_mprts_xyz);
    ops->push_mprts_xyz(push, mprts, mflds);
  } else {
    assert(ops->push_mprts_1);
    ops->push_mprts_1(push, mprts, mflds);
  }

  }
  
  psc_stats_stop(st_time_particle);
  prof_stop(pr_time_step_no_comm);
  prof_stop(pr);

  if (ops->particles_type) {
    psc_mparticles_put_as(mprts, mprts_base, 0);
  }
  if (ops->fields_type) {
    psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);
  }
}

void
psc_push_particles_stagger(struct psc_push_particles *push,
			   struct psc_mparticles *mprts_base, struct psc_mfields *mflds_base)
{
  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);

  struct psc_mparticles *mprts;
  if (ops->particles_type) {
    mprts = psc_mparticles_get_as(mprts_base, ops->particles_type, 0);
  } else {
    mprts = mprts_base;
  }
  struct psc_mfields *mflds;
  if (ops->fields_type) {
    mflds = psc_mfields_get_as(mflds_base, ops->fields_type, EX, EX + 6);
  } else {
    mflds = mflds_base;
  }
  
  int *im = ppsc->domain.gdims;

  if (im[0] == 1 && im[1] > 1 && im[2] > 1) { // yz
    if (!ops->stagger_mprts_yz) {
      mprintf("WARNING: no stagger_mprts() method!\n");
    } else {
      ops->stagger_mprts_yz(push, mprts, mflds);
    }
  } else if (im[0] == 1 && im[1] == 1 && im[2] == 1) { // 1
    if (!ops->stagger_mprts_1) {
      mprintf("WARNING: no stagger_mprts() method!\n");
    } else {
      ops->stagger_mprts_1(push, mprts, mflds);
    }
  } else {
    mprintf("WARNING: no stagger_mprts() case!\n");
  }

  if (ops->particles_type) {
    psc_mparticles_put_as(mprts, mprts_base, 0);
  }
  if (ops->fields_type) {
    psc_mfields_put_as(mflds, mflds_base, JXI, JXI + 3);
  }
}

unsigned int
psc_push_particles_get_mp_flags(struct psc_push_particles *push)
{
  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);
  return ops->mp_flags;
}

// ======================================================================
// psc_push_particles_init

static void
psc_push_particles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_generic_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_2nd_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1p5_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1st_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1sff_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb2_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec_single_by_block_ops);
#ifdef USE_FORTRAN
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_fortran_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_vay_ops);
#endif
#ifdef USE_SSE2
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_ps_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_ps2_ops);
#endif
#ifdef USE_CBE
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_cbe_ops);
#endif
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_4x4_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec3d_2x2_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec3d_4x4_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec3d_8x8_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec3d_2x2_gmem_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec3d_4x4_gmem_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec3d_8x8_gmem_cuda_ops);
#endif
#ifdef USE_CUDA2
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec_cuda2_host_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec_cuda2_ops);
#endif
#ifdef USE_ACC
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec_acc_ops);
#endif
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_vpic_ops);
#endif
}

// ======================================================================
// psc_push_particles class

struct mrc_class_psc_push_particles mrc_class_psc_push_particles = {
  .name             = "psc_push_particles",
  .size             = sizeof(struct psc_push_particles),
  .init             = psc_push_particles_init,
};

