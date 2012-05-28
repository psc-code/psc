
#include "psc_bnd_cuda.h"
#include "psc_cuda.h"
#include "psc_bnd_private.h"
#include "psc_bnd_cuda_fields.h"

#include <mrc_domain.h>
#include <mrc_ddc.h>
#include <mrc_profile.h>
#include <string.h>

// OPT opportunities:
// general: make fields really 2D in CUDA and in bnd xchg
// reduce packed buffers even more, e.g., add only on cuda side finally

void psc_bnd_cuda_setup(struct psc_bnd *bnd);
void psc_bnd_cuda_unsetup(struct psc_bnd *bnd);
void psc_bnd_cuda_exchange_particles(struct psc_bnd *bnd,
					  mparticles_base_t *particles_base);

// ======================================================================
// psc_bnd: subclass "cuda"

struct psc_bnd_ops psc_bnd_cuda_ops = {
  .name                  = "cuda",
  .size                  = sizeof(struct psc_bnd_cuda),
  .setup                 = psc_bnd_cuda_setup,
  .unsetup               = psc_bnd_cuda_unsetup,
  .exchange_particles    = psc_bnd_cuda_exchange_particles,

  .create_ddc            = psc_bnd_fields_cuda_create,
  .add_ghosts            = psc_bnd_fields_cuda_add_ghosts,
  .fill_ghosts           = psc_bnd_fields_cuda_fill_ghosts,
};

