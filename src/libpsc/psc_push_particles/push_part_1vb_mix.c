
#include "psc_push_particles_private.h"

#include <string.h>

static struct psc_push_particles_ops *mix_ops[] = {
  &psc_push_particles_1vb_4x4_cuda_ops,
  &psc_push_particles_1vb_ps_ops,
  NULL,
};

static struct psc_particles_ops *mix_prts_ops[] = {
  &psc_particles_cuda_ops,
  &psc_particles_single_ops,
  NULL,
};

static void
psc_push_particles_1vb_mix_push_mprts_yz(struct psc_push_particles *push,
					 struct psc_mparticles *mprts,
					 struct psc_mfields *mflds)
{
  for (int i = 0; mix_ops[i]; i++) {
    if (mix_ops[i]->push_mprts_yz) {
      // FIXME, push isn't really the push we're calling...
      mix_ops[i]->push_mprts_yz(push, mprts, mflds);
    } else {
      assert(mix_ops[i]->push_a_yz);
      for (int p = 0; p < mprts->nr_patches; p++) {
	struct psc_particles *prts = psc_mparticles_get_patch(mprts, p);
	struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
	if (psc_particles_ops(prts) == mix_prts_ops[i]) {
	  mix_ops[i]->push_a_yz(push, prts, flds);
	}
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

