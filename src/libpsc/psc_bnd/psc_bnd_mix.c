
#include "psc_bnd_private.h"
#include "psc_bnd_fld.h"
#include "psc_particles_mix.h"
#include "psc_particles_single.h"
#include "psc_particles_cuda.h"

#include <mrc_profile.h>
#include <string.h>

// ======================================================================
// psc_bnd: subclass "mix"

struct psc_bnd_ops psc_bnd_mix_ops = {
  .name                    = "mix",
  .create_ddc              = psc_bnd_fld_mix_create,
  .add_ghosts              = psc_bnd_fld_mix_add_ghosts,
  .fill_ghosts             = psc_bnd_fld_mix_fill_ghosts,
};
