
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
static void
sse2_particles_from_fortran(void)
{
  struct psc_sse2 *sse2 = psc.c_ctx;
  int n_part = psc.pp.n_part;
  int pad = 0;
  if((n_part % VEC_SIZE) != 0){
    pad = VEC_SIZE - (n_part % VEC_SIZE);
  }
  
  if(n_part > sse2->part_allocated) {
    free(sse2->part); // Is this safe? ie, does free(NULL) do nothing as part of the standard?
    sse2->part_allocated = n_part * 1.2;
    if(n_part*0.2 < pad){
      sse2->part_allocated += pad;
    }
    sse2->part = calloc(sse2->part_allocated, sizeof(*sse2->part));
  }

  for (int n = 0; n < n_part; n++) {
    particle_base_t *f_part = psc_particles_base_get_one(&psc.pp, n);
    struct sse2_particle *part = &sse2->part[n];

    part->xi  = f_part->xi;
    part->yi  = f_part->yi;
    part->zi  = f_part->zi;
    part->pxi = f_part->pxi;
    part->pyi = f_part->pyi;
    part->pzi = f_part->pzi;
    part->qni = f_part->qni;
    part->mni = f_part->mni;
    part->wni = f_part->wni;
    assert(round(part->xi) == 0); ///< \FIXME This assert only fits for the yz pusher. Idealy, no assert would be needed here, but until we can promise 'true 2D' some check is needed.
  }
  // We need to give the padding a non-zero mass to avoid NaNs
  for(int n = n_part; n < (n_part + pad); n++){
    particle_base_t *f_part = psc_particles_base_get_one(&psc.pp, n_part-1);
    struct sse2_particle *part = &sse2->part[n];
    part->xi  = f_part->xi; //We need to be sure the padding loads fields inside the local domain
    part->yi  = f_part->yi;
    part->zi  = f_part->zi;
    part->mni = 1.0;
  }
}

/// Copy particles from SSE2 data structures to fortran structures.
static void
sse2_particles_to_fortran()
{
   struct psc_sse2 *sse2 = psc.c_ctx;
   
   for(int n = 0; n < psc.pp.n_part; n++) {
     particle_base_t *f_part = psc_particles_base_get_one(&psc.pp, n);
     struct sse2_particle *part = &sse2->part[n];
     
     f_part->xi  = part->xi;
     f_part->yi  = part->yi;
     f_part->zi  = part->zi;
     f_part->pxi = part->pxi;
     f_part->pyi = part->pyi;
     f_part->pzi = part->pzi;
     f_part->qni = part->qni;
     f_part->mni = part->mni;
     f_part->wni = part->wni;
   }
}

/// Copy fields from Fortran data structures to an SSE2 friendly format.
static void
sse2_fields_from_fortran(){
  struct psc_sse2 *sse2 = psc.c_ctx;
  sse2->fields = _mm_malloc(NR_FIELDS*psc.fld_size*sizeof(sse2_real), 16);

  for(int m = 0; m < NR_FIELDS; m++){
    for(int n = 0; n < psc.fld_size; n++){
      //preserve Fortran ordering for now
      sse2->fields[m * psc.fld_size + n] = (sse2_real)psc.pf.flds[m][n];
    }
  }
}

/// Copy fields from SSE2 data structures into Fortran structures.
static void
sse2_fields_to_fortran(){
  struct psc_sse2 *sse2 = psc.c_ctx;
  assert(sse2->fields != NULL);
  for(int m = 0; m < NR_FIELDS; m++){
    for(int n = 0; n < psc.fld_size; n++){
      psc.pf.flds[m][n] = sse2->fields[m * psc.fld_size + n];
    }
  }
  _mm_free(sse2->fields);
}

/// Pointers to functions optimized for SSE2
struct psc_ops psc_ops_sse2 = {
  .name = "sse2",
  .create                 = sse2_create,
  .destroy                = sse2_destroy,
  .particles_from_fortran = sse2_particles_from_fortran,
  .particles_to_fortran   = sse2_particles_to_fortran,
  .fields_from_fortran    = sse2_fields_from_fortran,
  .fields_to_fortran      = sse2_fields_to_fortran,
  .push_part_yz_a         = sse2_push_part_yz_a,
  .push_part_yz_b         = sse2_push_part_yz_b,
  .push_part_yz           = sse2_push_part_yz,
}; 

/// \file psc_sse2.c Backend functions for SSE2 implementation.
