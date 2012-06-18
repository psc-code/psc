
#include "psc_push_particles_private.h"

#include "psc_particles_mix.h"


#include <string.h>

static struct psc_mparticles_ops *mix_mprts_ops[] = {
  &psc_mparticles_cuda_ops,
  &psc_mparticles_single_ops,
  NULL,
};

static struct psc_push_particles_ops *mix_ops[] = {
  &psc_push_particles_1vb_4x4_cuda_ops,
  &psc_push_particles_1vb_ps2_ops,
  NULL,
};

static void
psc_push_particles_1vb_mix_push_mprts_yz(struct psc_push_particles *push,
					 struct psc_mparticles *mprts,
					 struct psc_mfields *mflds)
{
  assert(psc_mparticles_ops(mprts) == &psc_mparticles_mix_ops);
  struct psc_mparticles_mix *mix = psc_mparticles_mix(mprts);

  for (int i = 0; mix_mprts_ops[i]; i++) {
    if (psc_mparticles_ops(mix->sub) != mix_mprts_ops[i]) {
      continue;
    }
    if (mix_ops[i]->push_mprts_yz) {
      mix_ops[i]->push_mprts_yz(NULL, mix->sub, mflds);
    } else {
      assert(mix_ops[i]->push_a_yz);
      for (int p = 0; p < mprts->nr_patches; p++) {
	struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
	struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
	mix_ops[i]->push_a_yz(NULL, prts, flds);
      }
    }
  }
}

// ======================================================================
// psc_push_particles: subclass "1vb_mix"

struct psc_push_particles_ops psc_push_particles_1vb_mix_ops = {
  .name                  = "1vb_mix",
  .push_mprts_yz         = psc_push_particles_1vb_mix_push_mprts_yz,
};

