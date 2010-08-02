
#include "psc.h"
#include "psc_particles_fortran.h"

#include <stdlib.h>

static bool __alloced;

void
psc_particles_fortran_alloc(psc_particles_fortran_t *pp, int n_part)
{
  if (!__alloced) {
    __alloced = true;
    pp->particles = ALLOC_particles(n_part);
  } else {
    // FIXME, realloc also copies all particles over,
    // which is not needed when fortran is not the particles base type,
    // but we just use this function to alloc temp storage
    pp->particles = REALLOC_particles(n_part);
  }
}

void
psc_particles_fortran_realloc(psc_particles_fortran_t *pp, int new_n_part)
{
  assert(__alloced);
  pp->particles = REALLOC_particles(new_n_part);
}

void
psc_particles_fortran_free(psc_particles_fortran_t *pp)
{
  assert(__alloced);
  FREE_particles();
  pp->particles = NULL;
  __alloced = false;
}

void
psc_particles_fortran_get(psc_particles_fortran_t *pp)
{
  pp->particles = psc.pp.particles;
  pp->n_part = psc.pp.n_part;
}

void
psc_particles_fortran_put(psc_particles_fortran_t *pp)
{
  psc.pp.n_part = pp->n_part;
  psc.pp.particles = pp->particles;
}

