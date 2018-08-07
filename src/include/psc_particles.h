
#ifndef PSC_PARTICLES_H
#define PSC_PARTICLES_H

#include <mrc_obj.h>
#include <psc_bits.h>
#include <assert.h>

#include "grid.hxx"
#include "particles_traits.hxx"

struct particle_inject
{
  using real_t = double; // FIXME, do we want/need to keep this?
  
  double x[3];
  double u[3];
  double w;
  int kind;
};

#define MP_DONT_COPY (0x1)
#define MP_DONT_RESIZE (0x2)

#endif
