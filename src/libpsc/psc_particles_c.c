
#include "psc.h"
#include "psc_particles_as_c.h"

#include "psc_particles_common.c"

void
particles_c_realloc(struct psc_particles *prts, int new_n_part)
{
  struct psc_particles_c *c = psc_particles_c(prts);
  if (new_n_part <= psc_mparticles_n_alloced(prts->mprts, prts->p))
    return;

  int n_alloced = new_n_part * 1.2;
  psc_mparticles_set_n_alloced(prts->mprts, prts->p, n_alloced);
  c->particles = realloc(c->particles, n_alloced * sizeof(*c->particles));
}

// ======================================================================
// psc_mparticles: subclass "c"
  
struct psc_mparticles_ops psc_mparticles_c_ops = {
  .name                    = "c",
};

// ======================================================================
// psc_particles: subclass "c"

struct psc_particles_ops psc_particles_c_ops = {
  .name                    = "c",
  .size                    = sizeof(struct psc_particles_c),
  .setup                   = psc_particles_c_setup,
  .destroy                 = psc_particles_c_destroy,
};
