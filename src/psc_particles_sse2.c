
#include "psc.h"
#include "psc_particles_sse2.h"

#include <math.h>
#include <stdlib.h>
#include <assert.h>

static size_t __sse2_part_allocated;
static particle_sse2_t *__sse2_part_data;

/// Copy particles from base data structures to an SSE2 friendly format.
void
psc_particles_sse2_get(psc_particles_sse2_t *particles)
{
  int n_part = psc.pp.n_part;
  int pad = 0;
  if((n_part % VEC_SIZE) != 0){
    pad = VEC_SIZE - (n_part % VEC_SIZE);
  }
  
  if(n_part > __sse2_part_allocated) {
    free(__sse2_part_data);
    __sse2_part_allocated = n_part * 1.2;
    if(n_part*0.2 < pad){
      __sse2_part_allocated += pad;
    }
    __sse2_part_data = calloc(__sse2_part_allocated, sizeof(*__sse2_part_data));
  }
  particles->particles = __sse2_part_data;
  particles->n_part = n_part;

  for (int n = 0; n < n_part; n++) {
    particle_base_t *base_part = psc_particles_base_get_one(&psc.pp, n);
    particle_sse2_t *part = &particles->particles[n];

    part->xi  = base_part->xi;
    part->yi  = base_part->yi;
    part->zi  = base_part->zi;
    part->pxi = base_part->pxi;
    part->pyi = base_part->pyi;
    part->pzi = base_part->pzi;
    part->qni = base_part->qni;
    part->mni = base_part->mni;
    part->wni = base_part->wni;
    assert(round(part->xi) == 0); ///< \FIXME This assert only fits for the yz pusher. Idealy, no assert would be needed here, but until we can promise 'true 2D' some check is needed.
  }
  // We need to give the padding a non-zero mass to avoid NaNs
  for(int n = n_part; n < (n_part + pad); n++){
    particle_base_t *base_part = psc_particles_base_get_one(&psc.pp, n_part - 1);
    particle_sse2_t *part = &particles->particles[n];
    part->xi  = base_part->xi; //We need to be sure the padding loads fields inside the local domain
    part->yi  = base_part->yi;
    part->zi  = base_part->zi;
    part->mni = 1.0;
  }
}

/// Copy particles from SSE2 data structures to base structures.
void
psc_particles_sse2_put(psc_particles_sse2_t *particles)
{
   for(int n = 0; n < psc.pp.n_part; n++) {
     particle_base_t *base_part = psc_particles_base_get_one(&psc.pp, n);
     particle_sse2_t *part = &particles->particles[n];
     
     base_part->xi  = part->xi;
     base_part->yi  = part->yi;
     base_part->zi  = part->zi;
     base_part->pxi = part->pxi;
     base_part->pyi = part->pyi;
     base_part->pzi = part->pzi;
     base_part->qni = part->qni;
     base_part->mni = part->mni;
     base_part->wni = part->wni;
   }
}

