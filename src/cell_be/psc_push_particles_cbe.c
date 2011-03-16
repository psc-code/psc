#include "psc_push_particles_private.h"
#include "psc_ppu.h"

// ======================================================================
// psc_push_particles: subclass "cbe"

struct psc_push_particles_ops psc_push_particles_cbe_ops = {
  .name                   = "cbe",
  .push_xy                = psc_push_particles_cbe_push_xy,
};
