
#include "psc.h"
#include "psc_fields_fortran.h"
#include "psc_fields_c.h"
#include "psc_glue.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>

void
fields_fortran_alloc(struct psc_fields *pf, int ib[3], int ie[3], int nr_comp,
		     int first_comp)
{
  pf->data = calloc(nr_comp, sizeof(fields_fortran_real_t *));

  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    pf->ib[d] = ib[d];
    pf->im[d] = ie[d] - ib[d];
    size *= pf->im[d];
  }
  pf->nr_comp = nr_comp;
  pf->first_comp = first_comp;

  fields_fortran_real_t **flds = pf->data;
  flds[0] = calloc(size * nr_comp, sizeof(flds[0]));
  for (int i = 1; i < nr_comp; i++) {
    flds[i] = flds[0] + i * size;
  }
}

static void
psc_fields_fortran_setup(struct psc_fields *pf)
{
}

static void
psc_fields_fortran_destroy(struct psc_fields *pf)
{
  fields_fortran_real_t **flds = pf->data;
  free(flds[0]);

  for (int i = 0; i < pf->nr_comp; i++) {
    flds[i] = NULL;
  }
  free(flds);
}

void
fields_fortran_zero(struct psc_fields *pf, int m)
{
  fields_fortran_real_t **flds = pf->data;
  memset(flds[m], 0, 
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_fortran_real_t));
}

void
fields_fortran_set(struct psc_fields *pf, int m, fields_fortran_real_t val)
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
fields_fortran_copy(struct psc_fields *pf, int m_to, int m_from)
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
fields_fortran_axpy(struct psc_fields *y, fields_fortran_real_t a,
		    struct psc_fields *x)
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
fields_fortran_scale(struct psc_fields *pf, fields_fortran_real_t val)
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
  psc_mfields_setup_super(flds);

  struct mrc_patch *patches = mrc_domain_get_patches(flds->domain,
						     &flds->nr_patches);
  flds->flds = calloc(flds->nr_patches, sizeof(*flds->flds));
  for (int p = 0; p < flds->nr_patches; p++) {
    struct psc_fields *pf = psc_fields_create(psc_mfields_comm(flds));
    psc_fields_set_type(pf, "fortran");
    psc_fields_setup(pf);
    flds->flds[p] = pf;
    int ilg[3] = { -flds->ibn[0], -flds->ibn[1], -flds->ibn[2] };
    int ihg[3] = { patches[p].ldims[0] + flds->ibn[0],
		   patches[p].ldims[1] + flds->ibn[1],
		   patches[p].ldims[2] + flds->ibn[2] };
    fields_fortran_alloc(pf, ilg, ihg, flds->nr_fields, flds->first_comp);
  }
}

static void
_psc_mfields_fortran_destroy(mfields_fortran_t *flds)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    struct psc_fields *pf = psc_mfields_get_patch_fortran(flds, p);
    psc_fields_destroy(pf);
  }
  free(flds->flds);
}									
									
void
psc_mfields_fortran_copy_to_c(mfields_fortran_t *flds_fortran, mfields_c_t *flds_c, int mb, int me)
{
  psc_foreach_patch(ppsc, p) {
    struct psc_fields *pf_c = psc_mfields_get_patch_c(flds_c, p);
    struct psc_fields *pf_base = psc_mfields_get_patch_fortran(flds_fortran, p);
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
    struct psc_fields *pf = psc_mfields_get_patch_fortran(flds, p);
    struct psc_fields *pf_c = psc_mfields_get_patch_c(flds_base, p);
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

// ======================================================================
// psc_fields: subclass "fortran"
  
struct psc_fields_ops psc_fields_fortran_ops = {
  .name                  = "fortran",
  .setup                 = psc_fields_fortran_setup,
  .destroy               = psc_fields_fortran_destroy,
};

