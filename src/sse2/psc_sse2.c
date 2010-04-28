
#include "psc_sse2.h"

#include <stdlib.h>
#include <assert.h>

static void
sse2_create()
{
  struct psc_sse2 *sse2 = malloc(sizeof(*sse2));
  psc.c_ctx = sse2;
}

static void
sse2_destroy()
{
  struct psc_sse2 *sse2 = psc.c_ctx;
  free(sse2);
}

// For now this is all identical to kai's generic_c. 
static void
sse2_particles_from_fortran()
{
  struct psc_sse2 *sse2 = psc.c_ctx;
  sse2->part = calloc(psc.n_part, sizeof(*sse2->part));
  for (int n = 0; n < psc.n_part; n++) {
    struct f_particle *f_part = &psc.f_part[n];
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
  }
}

static void
sse2_particles_to_fortran()
{
   struct psc_sse2 *sse2 = psc.c_ctx;
   
   for(int n = 0; n < psc.n_part; n++) {
     struct f_particle *f_part = &psc.f_part[n];
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

static void
sse2_fields_from_fortran(){
  struct psc_sse2 *sse2 = psc.c_ctx;
  sse2->fields = calloc(NR_FIELDS*psc.fld_size, sizeof(real));
  assert(sse2->fields != NULL);
  //It's a dirty job, but somebody's gotta do it
  for(int m = 0; m < NR_FIELDS; m++){
    for(int n = 0; n < psc.fld_size; n++){
      //preserve Fortran ordering for now
      sse2->fields[m * psc.fld_size + n] = (float)psc.f_fields[m][n];
    }
  }
}


struct psc_ops psc_ops_sse2 = {
  .name = "sse2",
  .create                 = sse2_create,
  .destroy                = sse2_destroy,
  .particles_from_fortran = sse2_particles_from_fortran,
  .particles_to_fortran   = sse2_particles_to_fortran,
  .fields_from_fortran    = sse2_fields_from_fortran,
  .push_part_yz_a         = sse2_push_part_yz_a,
  .push_part_yz_b         = sse2_push_part_yz_b,
}; 
