
#include "psc_bnd_particles_private.h"
#include "../psc_bnd/ddc_particles.h"
#include "psc_particles_mix.h"
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
    cuda->bnd_prts = realloc(cuda->bnd_prts, new_n_particles * sizeof(*cuda->bnd_prts));
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
// psc_bnd_particles_sub_setup

static void
psc_bnd_particles_sub_setup(struct psc_bnd_particles *bnd)
{
  psc_bnd_particles_setup_super(bnd);
  bnd->ddcp = ddc_particles_create(bnd->psc->mrc_domain, sizeof(particle_single_t),
				   sizeof(particle_single_real_t),
				   MPI_PARTICLES_SINGLE_REAL,
				   ddcp_particles_realloc,
				   ddcp_particles_get_addr);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_unsetup

static void
psc_bnd_particles_sub_unsetup(struct psc_bnd_particles *bnd)
{
  ddc_particles_destroy(bnd->ddcp);
  MHERE;
  // FIXME, the whole setup/unsetup business is broken in psc_bnd_particles in general
}

// ----------------------------------------------------------------------

static struct psc_mparticles_ops *mix_mprts_ops[] = {
  &psc_mparticles_cuda_ops,
  &psc_mparticles_single_ops,
  NULL,
};

static struct psc_bnd_particles_ops *mix_bnd_ops[] = {
  &psc_bnd_particles_cuda_ops,
  &psc_bnd_particles_single2_ops,
  NULL,
};

// ----------------------------------------------------------------------
// psc_bnd_particles_sub_exchange_particles

static void
psc_bnd_particles_sub_exchange_particles(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts)
{
  assert(psc_mparticles_ops(mprts) == &psc_mparticles_mix_ops);
  struct psc_mparticles_mix *mix = psc_mparticles_mix(mprts);

  static int pr_1, pr_2, pr_3;
  if (!pr_1) {
    pr_1 = prof_register("xchg_prep", 1., 0, 0);
    pr_2 = prof_register("xchg_comm", 1., 0, 0);
    pr_3 = prof_register("xchg_post", 1., 0, 0);
  }

  prof_start(pr_1);
  for (int i = 0; mix_mprts_ops[i]; i++) {
    if (psc_mparticles_ops(mix->sub) != mix_mprts_ops[i])
      continue;

    if (mix_bnd_ops[i]->exchange_mprts_prep) {
      // FIXME, not passing the right bnd object
      mix_bnd_ops[i]->exchange_mprts_prep(bnd, mix->sub);
    } else {
      assert(mix_bnd_ops[i]->exchange_particles_prep);
      for (int p = 0; p < mprts->nr_patches; p++) {
	struct psc_particles *prts = psc_mparticles_get_patch(mix->sub, p);
	mix_bnd_ops[i]->exchange_particles_prep(bnd, prts);
      }
    }
  }

  prof_stop(pr_1);

  prof_start(pr_2);
  ddc_particles_comm(bnd->ddcp, mprts);
  prof_stop(pr_2);

  prof_start(pr_3);
  for (int i = 0; mix_mprts_ops[i]; i++) {
    if (psc_mparticles_ops(mix->sub) != mix_mprts_ops[i])
      continue;

    if (mix_bnd_ops[i]->exchange_mprts_post) {
      // FIXME, not passing the right bnd object
      mix_bnd_ops[i]->exchange_mprts_post(bnd, mix->sub);
    } else {
      assert(mix_bnd_ops[i]->exchange_particles_post);
      for (int p = 0; p < mprts->nr_patches; p++) {
	struct psc_particles *prts = psc_mparticles_get_patch(mix->sub, p);
	mix_bnd_ops[i]->exchange_particles_post(bnd, prts);
      }
    }
  }
  prof_stop(pr_3);
}

// ======================================================================
// psc_bnd_particles: subclass "mix"

struct psc_bnd_particles_ops psc_bnd_particles_mix_ops = {
  .name                    = "mix",
  .setup                   = psc_bnd_particles_sub_setup,
  .unsetup                 = psc_bnd_particles_sub_unsetup,
  .exchange_particles      = psc_bnd_particles_sub_exchange_particles,
};
