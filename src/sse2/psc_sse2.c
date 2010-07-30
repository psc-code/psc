
#include "psc_sse2.h"

#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "psc_sse2.h"
#include "simd_wrap.h"

//---------------------------------------------
/// Fills the global constant numerical vectors used by sse2. 

/// **Must** be called before they are used
void 
init_vec_numbers(void) {		
  ones.r = pv_set1_real(1.0);			
  half.r = pv_set1_real(.5);			
  onepfive.r = pv_set1_real(1.5);		
  threefourths.r = pv_set1_real(.75);		
  third.r = pv_set1_real(1./3.);		
  ione.r = pv_set1_int(1);			
}

/// Allocate mememory for the SSE2 context.
static void
sse2_create(void)
{
  struct psc_sse2 *sse2 = malloc(sizeof(*sse2));
  memset(sse2, 0, sizeof(*sse2));
  psc.c_ctx = sse2;
}

/// Cleanup the SSE2 context.
static void
sse2_destroy(void)
{
  struct psc_sse2 *sse2 = psc.c_ctx;
  free(sse2);
}

// For now this is all more or less identical to kai's generic_c. 
/// Copy particles from Fortran data structures to an SSE2 friendly format.
void
sse2_particles_from_fortran(struct psc_sse2 *sse2)
{
  int n_part = psc.pp.n_part;
  int pad = 0;
  if((n_part % VEC_SIZE) != 0){
    pad = VEC_SIZE - (n_part % VEC_SIZE);
  }
  
  if(n_part > sse2->part_allocated) {
    free(sse2->part.particles); // Is this safe? ie, does free(NULL) do nothing as part of the standard?
    sse2->part_allocated = n_part * 1.2;
    if(n_part*0.2 < pad){
      sse2->part_allocated += pad;
    }
    sse2->part.particles = calloc(sse2->part_allocated, sizeof(*sse2->part.particles));
  }

  for (int n = 0; n < n_part; n++) {
    particle_base_t *base_part = psc_particles_base_get_one(&psc.pp, n);
    struct sse2_particle *part = &sse2->part.particles[n];

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
    struct sse2_particle *part = &sse2->part.particles[n];
    part->xi  = base_part->xi; //We need to be sure the padding loads fields inside the local domain
    part->yi  = base_part->yi;
    part->zi  = base_part->zi;
    part->mni = 1.0;
  }
}

/// Copy particles from SSE2 data structures to fortran structures.
void
sse2_particles_to_fortran(struct psc_sse2 *sse2)
{
   for(int n = 0; n < psc.pp.n_part; n++) {
     particle_base_t *base_part = psc_particles_base_get_one(&psc.pp, n);
     struct sse2_particle *part = &sse2->part.particles[n];
     
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

/// Copy fields from Fortran data structures to an SSE2 friendly format.
void
sse2_fields_from_fortran(struct psc_sse2 *sse2)
{
  sse2->flds.flds = _mm_malloc(NR_FIELDS*psc.fld_size*sizeof(sse2_real), 16);

  int *ilg = psc.ilg;
  for(int m = 0; m < NR_FIELDS; m++){
    for(int n = 0; n < psc.fld_size; n++){
      //preserve Fortran ordering for now
      sse2->flds.flds[m * psc.fld_size + n] =
	(sse2_real) ((&F3_BASE(m, ilg[0],ilg[1],ilg[2]))[n]);
    }
  }
}

/// Copy fields from SSE2 data structures into Fortran structures.
void
sse2_fields_to_fortran(struct psc_sse2 *sse2)
{
  assert(sse2->flds.flds != NULL);

  int *ilg = psc.ilg;
  for(int m = 0; m < NR_FIELDS; m++){
    for(int n = 0; n < psc.fld_size; n++){
      ((&F3_BASE(m, ilg[0],ilg[1],ilg[2]))[n]) = 
	sse2->flds.flds[m * psc.fld_size + n];
    }
  }
  _mm_free(sse2->flds.flds);
}

/// Pointers to functions optimized for SSE2

struct psc_ops psc_ops_sse2 = {
  .name                   = "sse2",
  .create                 = sse2_create,
  .destroy                = sse2_destroy,
  .push_part_yz_a         = sse2_push_part_yz_a,
  .push_part_yz_b         = sse2_push_part_yz_b,
  .push_part_yz           = sse2_push_part_yz,
}; 

/// \file psc_sse2.c Backend functions for SSE2 implementation.
