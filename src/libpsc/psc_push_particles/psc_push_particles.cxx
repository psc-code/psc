
#include "psc_push_particles_private.h"

#include "psc_push_fields_private.h"
#include "psc_bnd_fields.h"
#include "psc_bnd.h"
#include "push_particles.hxx"

#include <mrc_profile.h>

// ======================================================================
// forward to subclass

extern int pr_time_step_no_comm;
extern double *psc_balance_comp_time_by_patch;

void
psc_push_particles_prep(struct psc_push_particles *push,
			struct psc_mparticles *mprts_base, struct psc_mfields *mflds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("push_particles_prep", 1., 0, 0);
  }  

  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);

  prof_start(pr);
  prof_restart(pr_time_step_no_comm);
  psc_stats_start(st_time_particle);

  if (ops->prep) {
    ops->prep(push, mprts_base, mflds_base);
  }

  psc_stats_stop(st_time_particle);
  prof_stop(pr_time_step_no_comm);
  prof_stop(pr);
}

void
psc_push_particles_run(struct psc_push_particles *push,
		       struct psc_mparticles *mprts_base, struct psc_mfields *mflds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("push_particles_run", 1., 0, 0);
  }  

  PscPushParticlesBase pushp(push);
  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);

  struct psc_mparticles *mprts;
  if (ops->particles_type) {
    mprts = psc_mparticles_get_as(mprts_base, ops->particles_type, 0);
  } else {
    mprts = mprts_base;
  }

  prof_start(pr);
  prof_restart(pr_time_step_no_comm);
  psc_stats_start(st_time_particle);

  if (ops->push_mprts) {
    pushp->push_mprts(mprts, mflds_base);
  } else {

  int *im = ppsc->domain.gdims;

  if (im[0] > 1 && im[1] == 1 && im[2] == 1) { // x
    pushp->push_mprts_x(mprts, mflds_base);
  } else if (im[0] == 1 && im[1] > 1 && im[2] == 1) { // y
    pushp->push_mprts_y(mprts, mflds_base);
  } else if (im[0] == 1 && im[1] == 1 && im[2] > 1) { // y
    pushp->push_mprts_z(mprts, mflds_base);
  } else if (im[0] > 1 && im[1] > 1 && im[2] == 1) { // xy
    pushp->push_mprts_xy(mprts, mflds_base);
  } else if (im[0] > 1 && im[1] == 1 && im[2] > 1) { // xz
    pushp->push_mprts_xz(mprts, mflds_base);
  } else if (im[0] == 1 && im[1] > 1 && im[2] > 1) { // yz
    pushp->push_mprts_yz(mprts, mflds_base);
  } else if (im[0] > 1 && im[1] > 1 && im[2] > 1) { // xyz
    pushp->push_mprts_xyz(mprts, mflds_base);
  } else {
    pushp->push_mprts_1(mprts, mflds_base);
  }

  }
  
  psc_stats_stop(st_time_particle);
  prof_stop(pr_time_step_no_comm);
  prof_stop(pr);

  if (ops->particles_type) {
    psc_mparticles_put_as(mprts, mprts_base, 0);
  }
}

void
psc_push_particles_stagger(struct psc_push_particles *push,
			   struct psc_mparticles *mprts_base, struct psc_mfields *mflds_base)
{
  PscPushParticlesBase pushp(push);
  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);

  struct psc_mparticles *mprts;
  if (ops->particles_type) {
    mprts = psc_mparticles_get_as(mprts_base, ops->particles_type, 0);
  } else {
    mprts = mprts_base;
  }

  if (ops->stagger_mprts) {
    pushp->stagger_mprts(mprts, mflds_base);
  } else {

    int *im = ppsc->domain.gdims;

    if (im[0] == 1 && im[1] > 1 && im[2] > 1) { // yz
      pushp->stagger_mprts_yz(mprts, mflds_base);
    } else if (im[0] == 1 && im[1] == 1 && im[2] == 1) { // 1
      pushp->stagger_mprts_1(mprts, mflds_base);
    } else {
      mprintf("WARNING: no stagger_mprts() case!\n");
    }
  }
  
  if (ops->particles_type) {
    psc_mparticles_put_as(mprts, mprts_base, 0);
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

extern struct psc_push_particles_ops psc_push_particles_generic_c_ops;
extern struct psc_push_particles_ops psc_push_particles_2nd_double_ops;
extern struct psc_push_particles_ops psc_push_particles_1st_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_single_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_double_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_ps_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_ps2_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb2_single_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec_single_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec_double_ops;
extern struct psc_push_particles_ops psc_push_particles_fortran_ops;
extern struct psc_push_particles_ops psc_push_particles_vay_ops;
extern struct psc_push_particles_ops psc_push_particles_sse2_ops;
extern struct psc_push_particles_ops psc_push_particles_cbe_ops;
extern struct psc_push_particles_ops psc_push_particles_1vb_cuda_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec3d_cuda_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec3d_gmem_cuda_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec_cuda2_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec_cuda2_host_ops;
extern struct psc_push_particles_ops psc_push_particles_1vbec_acc_ops;
extern struct psc_push_particles_ops psc_push_particles_vpic_ops;

static void
psc_push_particles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_generic_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_2nd_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1st_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb2_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec_double_ops);
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
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec3d_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vbec3d_gmem_cuda_ops);
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

struct mrc_class_psc_push_particles_ : mrc_class_psc_push_particles {
  mrc_class_psc_push_particles_() {
    name             = "psc_push_particles";
    size             = sizeof(struct psc_push_particles);
    init             = psc_push_particles_init;
  }
} mrc_class_psc_push_particles;

