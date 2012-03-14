
#include "psc.h"
#include "psc_fields_fortran.h"
#include "psc_glue.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>

void
fields_fortran_alloc(fields_fortran_t *pf, int ib[3], int ie[3], int nr_comp,
		     int first_comp)
{
  pf->flds = calloc(nr_comp, sizeof(*pf->flds));

  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    pf->ib[d] = ib[d];
    pf->im[d] = ie[d] - ib[d];
    size *= pf->im[d];
  }
  pf->nr_comp = nr_comp;
  pf->first_comp = first_comp;

  pf->flds[0] = calloc(size * nr_comp, sizeof(*pf->flds[0]));
  for (int i = 1; i < nr_comp; i++) {
    pf->flds[i] = pf->flds[0] + i * size;
  }
}

void
fields_fortran_free(fields_fortran_t *pf)
{
  free(pf->flds[0]);

  for (int i = 0; i < pf->nr_comp; i++) {
    pf->flds[i] = NULL;
  }
  free(pf->flds);
}

void
fields_fortran_zero(fields_fortran_t *pf, int m)
{
  memset(pf->flds[m], 0, 
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_fortran_real_t));
}

void
fields_fortran_set(fields_fortran_t *pf, int m, fields_fortran_real_t val)
{
  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_FORTRAN(pf, m, jx, jy, jz) = val;
      }
    }
  }
}

void
fields_fortran_copy(fields_fortran_t *pf, int m_to, int m_from)
{
  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_FORTRAN(pf, m_to, jx,jy,jz) = F3_FORTRAN(pf, m_from, jx,jy,jz);
      }
    }
  }
}

void
fields_fortran_axpy(fields_fortran_t *y, fields_fortran_real_t a,
		    fields_fortran_t *x)
{
  assert(y->nr_comp == x->nr_comp);
  for (int m = 0; m < y->nr_comp; m++) {
    for (int jz = y->ib[2]; jz < y->ib[2] + y->im[2]; jz++) {
      for (int jy = y->ib[1]; jy < y->ib[1] + y->im[1]; jy++) {
	for (int jx = y->ib[0]; jx < y->ib[0] + y->im[0]; jx++) {
	  F3_FORTRAN(y, m, jx, jy, jz) += a * F3_FORTRAN(x, m, jx, jy, jz);
	}
      }
    }
  }
}

void
fields_fortran_scale(fields_fortran_t *pf, fields_fortran_real_t val)
{
  for (int m = 0; m < pf->nr_comp; m++) {
    for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
      for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
	for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	  F3_FORTRAN(pf, m, jx, jy, jz) *= val;
	}
      }
    }
  }
}

static void
_psc_mfields_fortran_setup(mfields_fortran_t *flds)
{
  struct mrc_patch *patches = mrc_domain_get_patches(flds->domain,
						     &flds->nr_patches);
  flds->data = calloc(flds->nr_patches, sizeof(fields_fortran_t));
  for (int p = 0; p < flds->nr_patches; p++) {			
    int ilg[3] = { -flds->ibn[0], -flds->ibn[1], -flds->ibn[2] };
    int ihg[3] = { patches[p].ldims[0] + flds->ibn[0],
		   patches[p].ldims[1] + flds->ibn[1],
		   patches[p].ldims[2] + flds->ibn[2] };
    fields_fortran_alloc(psc_mfields_get_patch_fortran(flds, p), ilg, ihg,
			 flds->nr_fields, flds->first_comp);
  }
}

static void
_psc_mfields_fortran_destroy(mfields_fortran_t *flds)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    fields_fortran_free(psc_mfields_get_patch_fortran(flds, p));
  }
  free(flds->data);
}									
									
void
psc_mfields_fortran_copy_to_c(mfields_fortran_t *flds_fortran, mfields_c_t *flds_c, int mb, int me)
{
  psc_foreach_patch(ppsc, p) {
    fields_c_t *pf_c = psc_mfields_get_patch_c(flds_c, p);
    fields_fortran_t *pf_base = psc_mfields_get_patch_fortran(flds_fortran, p);
    for (int m = mb; m < me; m++) {
      psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
	F3_C(pf_c, m, jx,jy,jz) = F3_FORTRAN(pf_base, m, jx,jy,jz);
      } foreach_3d_g_end;
    }
  }
}

void
psc_mfields_fortran_copy_from_c(mfields_fortran_t *flds, mfields_c_t *flds_base, int mb, int me)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    fields_fortran_t *pf = psc_mfields_get_patch_fortran(flds, p);
    fields_c_t *pf_c = psc_mfields_get_patch_c(flds_base, p);
    for (int m = mb; m < me; m++) {
      psc_foreach_3d_g(ppsc, p, jx, jy, jz) {
	F3_FORTRAN(pf, m, jx,jy,jz) = F3_C(pf_c, m, jx,jy,jz);
      } foreach_3d_g_end;
    }
  }
}

// ======================================================================
// psc_mfields: subclass "fortran"
  
struct psc_mfields_ops psc_mfields_fortran_ops = {
  .name                  = "fortran",
  .setup                 = _psc_mfields_fortran_setup,
  .destroy               = _psc_mfields_fortran_destroy,
  .copy_to_c             = psc_mfields_fortran_copy_to_c,
  .copy_from_c           = psc_mfields_fortran_copy_from_c,
};

