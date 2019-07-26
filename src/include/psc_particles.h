
#ifndef PSC_PARTICLES_H
#define PSC_PARTICLES_H

#include <mrc_obj.h>
#include <psc_bits.h>
#include <assert.h>

#include "grid.hxx"

#include "particle.h"

struct particle_inject
{
  using real_t = double;
  using Real3 = Vec3<real_t>;
  
  Real3 x;
  Real3 u;
  real_t w;
  int kind;
  psc::particle::Id id;
};

#define MP_DONT_COPY (0x1)

#endif
