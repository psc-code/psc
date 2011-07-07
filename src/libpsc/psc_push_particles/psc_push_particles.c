
#include "psc_push_particles_private.h"

// ======================================================================
// forward to subclass

void
psc_push_particles_run(struct psc_push_particles *push,
		       mparticles_base_t *particles, mfields_base_t *flds)
{
  psc_stats_start(st_time_particle);
  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);
  int *im = ppsc->domain.gdims;
  if (im[0] > 1 && im[1] > 1 && im[2] > 1) { // xyz
    assert(ops->push_xyz);
    ops->push_xyz(push, particles, flds);
  } else if (im[0] > 1 && im[2] > 1) { // xz
    assert(ops->push_xz);
    ops->push_xz(push, particles, flds);
  } else if (im[0] > 1 && im[1] > 1) { // xy
    assert(ops->push_xy);
    ops->push_xy(push, particles, flds);
  } else if (im[1] > 1 && im[2] > 1) { // yz
    assert(ops->push_yz);
    ops->push_yz(push, particles, flds);
  } else if (im[2] > 1) { // z
    assert(ops->push_z);
    ops->push_z(push, particles, flds);
  } else {
    assert(0);
  }
  psc_stats_stop(st_time_particle);
}

void
psc_push_particles_push_yz_a(struct psc_push_particles *push,
			     mparticles_base_t *particles, mfields_base_t *flds)
{
  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);
  assert(ops->push_yz_a);
  ops->push_yz_a(push, particles, flds);
}

void
psc_push_particles_push_yz_b(struct psc_push_particles *push,
			     mparticles_base_t *particles, mfields_base_t *flds)
{
  struct psc_push_particles_ops *ops = psc_push_particles_ops(push);
  assert(ops->push_yz_b);
  ops->push_yz_a(push, particles, flds);
}

// ======================================================================
// psc_push_particles_init

static void
psc_push_particles_init()
{
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_generic_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_fortran_ops);
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_vay_ops);
#ifdef USE_SSE2
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_sse2_ops);
#endif
#ifdef USE_CBE
  mrc_class_register_subclass(&mrc_class_psc_push_particles, &psc_push_particles_cbe_ops);
#endif

}

// ======================================================================
// psc_push_particles class

struct mrc_class_psc_push_particles mrc_class_psc_push_particles = {
  .name             = "psc_push_particles",
  .size             = sizeof(struct psc_push_particles),
  .init             = psc_push_particles_init,
};

