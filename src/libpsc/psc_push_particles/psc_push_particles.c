
#include "psc_push_particles_private.h"

#include "psc_push_fields_private.h"
#include "psc_bnd_fields.h"
#include "psc_bnd.h"

#include <mrc_profile.h>

// ======================================================================
// forward to subclass

static void
push_a_yz(struct psc_push_particles *push, struct psc_mparticles *mprts,
	  struct psc_mfields *mflds)
{
  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);

  if (ops->push_a_yz) {
    for (int p = 0; p < mprts->nr_patches; p++) {
      ops->push_a_yz(push, psc_mparticles_get_patch(mprts, p),
		     psc_mfields_get_patch(mflds, p));
    }
  } else {
    assert(ops->push_mprts_yz);
    ops->push_mprts_yz(push, mprts, mflds);
  }
}

void
psc_push_particles_run(struct psc_push_particles *push,
		       mparticles_base_t *mprts, mfields_base_t *mflds)
{
  static int pr;
  if (!pr) {
    pr = prof_register("push_particles_run", 1., 0, 0);
  }  

  prof_start(pr);
  psc_stats_start(st_time_particle);
  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);
  int *im = ppsc->domain.gdims;

  if (im[0] == 1 && im[1] > 1 && im[2] > 1) { // yz
    push_a_yz(push, mprts, mflds);
    return;
  }

  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    if (im[0] > 1 && im[1] > 1 && im[2] > 1) { // xyz
      assert(ops->push_a_xyz);
      ops->push_a_xyz(push, prts, flds);
    } else if (im[0] > 1 && im[2] > 1) { // xz
      assert(ops->push_a_xz);
      ops->push_a_xz(push, prts, flds);
    } else if (im[0] > 1 && im[1] > 1) { // xy
      assert(ops->push_a_xy);
      ops->push_a_xy(push, prts, flds);
    } else if (im[2] > 1) { // z
      assert(ops->push_a_z);
      ops->push_a_z(push, prts, flds);
    } else {
      assert(0);
    }
  }
  psc_stats_stop(st_time_particle);
  prof_stop(pr);
}

void
psc_push_particles_run_b(struct psc_push_particles *push,
			 mparticles_base_t *particles, mfields_base_t *mflds)
{
  psc_stats_start(st_time_particle);
  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);
  int *im = ppsc->domain.gdims;
  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    if (im[0] > 1 && im[1] > 1 && im[2] > 1) { // xyz
    } else if (im[0] > 1 && im[2] > 1) { // xz
    } else if (im[0] > 1 && im[1] > 1) { // xy
    } else if (im[1] > 1 && im[2] > 1) { // yz
      if (ops->push_b_yz) {
	ops->push_b_yz(push, prts, flds);
      }
    } else if (im[2] > 1) { // z
    } else {
    }
  }
  psc_stats_stop(st_time_particle);
}

void
psc_push_particles_calc_j(struct psc_push_particles *push,
			  mparticles_base_t *particles, mfields_base_t *mflds)
{
  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);

  bool have_calc_j = false;
  int *im = ppsc->domain.gdims;
  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    if (im[0] > 1 && im[1] > 1 && im[2] > 1) { // xyz
      if (ops->calc_j_xyz) {
	ops->calc_j_xyz(push, prts, flds);
	have_calc_j = true;
      }
    } else if (im[0] > 1 && im[2] > 1) { // xz
      if (ops->calc_j_xz) {
	ops->calc_j_xz(push, prts, flds);
	have_calc_j = true;
      }
    } else if (im[0] > 1 && im[1] > 1) { // xy
      if (ops->calc_j_xy) {
	ops->calc_j_xy(push, prts, flds);
	have_calc_j = true;
      }
    } else if (im[1] > 1 && im[2] > 1) { // yz
      if (ops->calc_j_yz) {
	ops->calc_j_yz(push, prts, flds);
	have_calc_j = true;
      }
    } else if (im[2] > 1) { // z
      if (ops->calc_j_z) {
	ops->calc_j_z(push, prts, flds);
	have_calc_j = true;
      }
    } else {
      assert(0);
    }
  }

  if (have_calc_j) {
    // add and fill ghost for J
    psc_bnd_fields_add_ghosts_J(ppsc->push_fields->bnd_fields, mflds);
    psc_bnd_add_ghosts(ppsc->bnd, mflds, JXI, JXI + 3);
    psc_bnd_fill_ghosts(ppsc->bnd, mflds, JXI, JXI + 3);
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
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1st_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1sff_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_double_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb2_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_fortran_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_vay_ops);
#ifdef USE_SSE2
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_ps_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_ps2_ops);
#endif
#ifdef USE_CBE
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_cbe_ops);
#endif
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_cuda_1st_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_2x2_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_4x4_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_8x8_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_mix_ops);
#endif
}

// ======================================================================
// psc_push_particles class

struct mrc_class_psc_push_particles mrc_class_psc_push_particles = {
  .name             = "psc_push_particles",
  .size             = sizeof(struct psc_push_particles),
  .init             = psc_push_particles_init,
};

