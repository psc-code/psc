
#include "psc_bnd_private.h"
#include "psc_bnd_fld.h"
#include "ddc_particles.h"
#include "psc_particles_single.h"
#include "psc_particles_cuda.h"

#include <mrc_profile.h>
#include <string.h>

// ----------------------------------------------------------------------
// ddcp_particles helpers

static void
ddcp_particles_realloc(void *_ctx, int p, int new_n_particles)
{
  struct psc_mparticles *particles = _ctx;
  struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
  if (psc_particles_ops(prts) == &psc_particles_single_ops) {
    particles_single_realloc(prts, new_n_particles);
  } else if (psc_particles_ops(prts) == &psc_particles_cuda_ops) {
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    assert(!cuda->bnd_prts);
    cuda->bnd_prts = malloc(new_n_particles * sizeof(*cuda->bnd_prts));
  } else {
    assert(0);
  }
}

static void *
ddcp_particles_get_addr(void *_ctx, int p, int n)
{
  struct psc_mparticles *particles = _ctx;
  struct psc_particles *prts = psc_mparticles_get_patch(particles, p);
  if (psc_particles_ops(prts) == &psc_particles_single_ops) {
    return particles_single_get_one(prts, n);
  } else if (psc_particles_ops(prts) == &psc_particles_cuda_ops) {
    struct psc_particles_cuda *cuda = psc_particles_cuda(prts);
    return &cuda->bnd_prts[n];
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_bnd_sub_setup

static void
psc_bnd_sub_setup(struct psc_bnd *bnd)
{
  psc_bnd_setup_super(bnd);
  bnd->ddcp = ddc_particles_create(bnd->ddc, sizeof(particle_single_t),
				   sizeof(particle_single_real_t),
				   MPI_PARTICLES_SINGLE_REAL,
				   ddcp_particles_realloc,
				   ddcp_particles_get_addr);
}

// ----------------------------------------------------------------------
// psc_bnd_sub_unsetup

static void
psc_bnd_sub_unsetup(struct psc_bnd *bnd)
{
  ddc_particles_destroy(bnd->ddcp);
  MHERE;
  // FIXME, the whole setup/unsetup business is broken in psc_bnd in general
}

static inline struct psc_bnd_ops *
get_ops(struct psc_particles *prts)
{
  if (psc_particles_ops(prts) == &psc_particles_single_ops) {
    return &psc_bnd_single2_ops;
  } else if (psc_particles_ops(prts) == &psc_particles_cuda_ops) {
    return &psc_bnd_cuda_ops;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_bnd_sub_exchange_particles

static void
psc_bnd_sub_exchange_particles(struct psc_bnd *bnd, struct psc_mparticles *mprts)
{
  static int pr_1, pr_2, pr_3;
  if (!pr_1) {
    pr_1 = prof_register("xchg_prep", 1., 0, 0);
    pr_2 = prof_register("xchg_comm", 1., 0, 0);
    pr_3 = prof_register("xchg_post", 1., 0, 0);
  }

  prof_start(pr_1);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_bnd_ops *ops = get_ops(prts);
    ops->exchange_particles_prep(bnd, prts);
  }
  prof_stop(pr_1);

  prof_start(pr_2);
  ddc_particles_comm(bnd->ddcp, mprts);
  prof_stop(pr_2);

  prof_start(pr_3);
  for (int p = 0; p < mprts->nr_patches; p++) {
    struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
    struct psc_bnd_ops *ops = get_ops(prts);
    ops->exchange_particles_post(bnd, prts);
  }
  prof_stop(pr_3);
}

// ======================================================================
// psc_bnd: subclass "mix"

struct psc_bnd_ops psc_bnd_mix_ops = {
  .name                    = "mix",
  .setup                   = psc_bnd_sub_setup,
  .unsetup                 = psc_bnd_sub_unsetup,
  .exchange_particles      = psc_bnd_sub_exchange_particles,

  .create_ddc              = psc_bnd_fld_mix_create,
  .add_ghosts              = psc_bnd_fld_mix_add_ghosts,
  .fill_ghosts             = psc_bnd_fld_mix_fill_ghosts,
};
