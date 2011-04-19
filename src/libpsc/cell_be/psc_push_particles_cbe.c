#include "psc_push_particles_private.h"
#include "psc_ppu.h"

// ======================================================================
// psc_push_particles: subclass "cbe"

static void 
cbe_push_setup(struct psc_push_particles *push)
{
  // Initialize the spes and create the context.
  printf("Particle setup called\n");
  psc_init_spes();
}

static void 
cbe_push_destroy(struct psc_push_particles *push)
{
  // The spes and free the context
  psc_kill_spes();
}

struct psc_push_particles_ops psc_push_particles_cbe_ops = {
  .name                   = "cbe",
  .push_xy                = psc_push_particles_cbe_push_xy,
  .setup                  = cbe_push_setup,
  .destroy                = cbe_push_destroy,
};
