
#include "psc_cuda.h"

#include <mrc_profile.h>
#include <math.h>

void
psc_push_particles_cuda_push_yz_a(struct psc_push_particles *push,
				  mparticles_base_t *particles_base,
				  mfields_base_t *flds_base)
{
  mfields_cuda_t flds;
  mparticles_cuda_t particles;
  psc_mfields_cuda_get_from(&flds, EX, EX + 6, flds_base);
  psc_mparticles_cuda_get_from(&particles, particles_base);

  static int pr;
  if (!pr) {
    pr = prof_register("cuda_part_yz_a", 1., 0, 0);
  }
  prof_start(pr);
  psc_foreach_patch(ppsc, p) {
    particles_cuda_t *pp = &particles.p[p];
    fields_cuda_t *pf = &flds.f[p];
    yz_a_set_constants(pp, pf);
    __cuda_push_part_yz_a(pp, pf);
  }
  prof_stop(pr);

  psc_mfields_cuda_put_to(&flds, JXI, JXI + 3, flds_base);
  psc_mparticles_cuda_put_to(&particles, particles_base);
}
