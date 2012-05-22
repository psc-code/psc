
#include "psc_push_particles_private.h"

// ======================================================================
// forward to subclass

void
psc_push_particles_run(struct psc_push_particles *push,
		       mparticles_base_t *particles, mfields_base_t *mflds)
{
  psc_stats_start(st_time_particle);
  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);
  int *im = ppsc->domain.gdims;
  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
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
    } else if (im[1] > 1 && im[2] > 1) { // yz
      assert(ops->push_a_yz);
      ops->push_a_yz(push, prts, flds);
    } else if (im[2] > 1) { // z
      assert(ops->push_a_z);
      ops->push_a_z(push, prts, flds);
    } else {
      assert(0);
    }
  }
  psc_stats_stop(st_time_particle);
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
  int *im = ppsc->domain.gdims;
  for (int p = 0; p < particles->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
    if (im[0] > 1 && im[1] > 1 && im[2] > 1) { // xyz
      if (ops->calc_j_xyz) {
	ops->calc_j_xyz(push, prts, flds);
      }
    } else if (im[0] > 1 && im[2] > 1) { // xz
      if (ops->calc_j_xz) {
	ops->calc_j_xz(push, prts, flds);
      }
    } else if (im[0] > 1 && im[1] > 1) { // xy
      if (ops->calc_j_xy) {
	ops->calc_j_xy(push, prts, flds);
      }
    } else if (im[1] > 1 && im[2] > 1) { // yz
      if (ops->calc_j_yz) {
	ops->calc_j_yz(push, prts, flds);
      }
    } else if (im[2] > 1) { // z
      if (ops->calc_j_z) {
	ops->calc_j_z(push, prts, flds);
      }
    } else {
      assert(0);
    }
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
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1vb_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_1sff_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_single_1vb_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_double_1vb_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_fortran_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_vay_ops);
#ifdef USE_SSE2
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_ps_1vb_ops);
#endif
#ifdef USE_CBE
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_cbe_ops);
#endif
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_cuda_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_cuda_1st_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_cuda_1vb_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_cuda_2x2_1vb_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_cuda_8x8_1vb_ops);
#endif
}

// ======================================================================
// psc_push_particles class

struct mrc_class_psc_push_particles mrc_class_psc_push_particles = {
  .name             = "psc_push_particles",
  .size             = sizeof(struct psc_push_particles),
  .init             = psc_push_particles_init,
};

