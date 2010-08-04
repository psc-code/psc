
#include "psc_generic_c.h"

#include <stdlib.h>
#include <string.h>

static particle_c_t *__arr;
static int __arr_size;

void
genc_particles_from_fortran(psc_particles_c_t *pp)
{
  if (psc.pp.n_part > __arr_size) {
    free(__arr);
    __arr = NULL;
  }
  if (!__arr) {
    __arr_size = psc.pp.n_part * 1.2;
    __arr = calloc(__arr_size, sizeof(*__arr));
  }

  pp->particles = __arr;
  for (int n = 0; n < psc.pp.n_part; n++) {
    particle_base_t *f_part = psc_particles_base_get_one(&psc.pp, n);
    particle_c_t *part = psc_particles_c_get_one(pp, n);

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

void
genc_particles_to_fortran(psc_particles_c_t *pp)
{
  for (int n = 0; n < psc.pp.n_part; n++) {
    particle_base_t *f_part = psc_particles_base_get_one(&psc.pp, n);
    particle_c_t *part = psc_particles_c_get_one(pp, n);

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

void
genc_fields_from_fortran(psc_fields_c_t *pf)
{
  pf->flds = calloc(NR_FIELDS * psc.fld_size, sizeof(*pf->flds));

  for (int m = EX; m <= BZ; m++) {
    for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
      for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
	for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	  F3_C(pf, m, jx,jy,jz) = F3_BASE(m, jx,jy,jz);
	}
      }
    }
  }
}

void
genc_fields_to_fortran(psc_fields_c_t *pf)
{
  for (int m = JXI; m <= JZI; m++) {
    for (int jz = psc.ilg[2]; jz < psc.ihg[2]; jz++) {
      for (int jy = psc.ilg[1]; jy < psc.ihg[1]; jy++) {
	for (int jx = psc.ilg[0]; jx < psc.ihg[0]; jx++) {
	  F3_BASE(m, jx,jy,jz) = F3_C(pf, m, jx,jy,jz);
	}
      }
    }
  }

  free(pf->flds);
}

struct psc_ops psc_ops_generic_c = {
  .name = "generic_c",
  .push_part_xz           = genc_push_part_xz,
  .push_part_yz_a         = genc_push_part_yz_a,
  .push_part_yz_b         = genc_push_part_yz_b,
};
