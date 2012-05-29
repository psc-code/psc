
#include "psc_push_particles_private.h"

#include <string.h>

static void
psc_push_particles_1vb_mix_push_a_yz(struct psc_push_particles *push,
				     struct psc_particles *prts_base,
				     struct psc_fields *flds_base)
{
  struct psc_push_particles_ops *ops;
  if (psc_particles_ops(prts_base) == &psc_particles_single_ops) {
    ops = &psc_push_particles_1vb_ps_ops;
  } else if (psc_particles_ops(prts_base) == &psc_particles_cuda_ops) {
    ops = &psc_push_particles_1vb_4x4_cuda_ops;
  } else {
    assert(0);
  }

  assert(ops->push_a_yz);
  ops->push_a_yz(push, prts_base, flds_base);
}

// ======================================================================
// psc_push_particles: subclass "1vb_mix"

struct psc_push_particles_ops psc_push_particles_1vb_mix_ops = {
  .name                  = "1vb_mix",
  .push_a_yz             = psc_push_particles_1vb_mix_push_a_yz,
};

