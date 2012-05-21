
#include "psc.h"
#include "psc_fields_single.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// FIXME, very duplicated from psc_fields_c.c

void
fields_single_alloc(fields_single_t *pf, int ib[3], int ie[3], int nr_comp, int first_comp)
{
  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    pf->ib[d] = ib[d];
    pf->im[d] = ie[d] - ib[d];
    size *= pf->im[d];
  }
  pf->nr_comp = nr_comp;
  pf->first_comp = first_comp;
  pf->flds = calloc(nr_comp * size, sizeof(*pf->flds));
}

void
fields_single_free(fields_single_t *pf)
{
  free(pf->flds);
}

void
fields_single_zero(fields_single_t *pf, int m)
{
  memset(&F3_S(pf, m, pf->ib[0], pf->ib[1], pf->ib[2]), 0,
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_single_real_t));
}

void
fields_single_set(fields_single_t *pf, int m, fields_single_real_t val)
{
  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_S(pf, m, jx, jy, jz) = val;
      }
    }
  }
}

void
fields_single_copy_comp(fields_single_t *pto, int m_to, fields_single_t *pfrom, int m_from)
{
  for (int jz = pto->ib[2]; jz < pto->ib[2] + pto->im[2]; jz++) {
    for (int jy = pto->ib[1]; jy < pto->ib[1] + pto->im[1]; jy++) {
      for (int jx = pto->ib[0]; jx < pto->ib[0] + pto->im[0]; jx++) {
	F3_S(pto, m_to, jx, jy, jz) = F3_S(pfrom, m_from, jx, jy, jz);
      }
    }
  }
}

void
fields_single_axpy(fields_single_t *y, fields_single_real_t a, fields_single_t *x)
{
  assert(y->nr_comp == x->nr_comp);
  for (int m = 0; m < y->nr_comp; m++) {
    for (int jz = y->ib[2]; jz < y->ib[2] + y->im[2]; jz++) {
      for (int jy = y->ib[1]; jy < y->ib[1] + y->im[1]; jy++) {
	for (int jx = y->ib[0]; jx < y->ib[0] + y->im[0]; jx++) {
	  F3_S(y, m, jx, jy, jz) += a * F3_S(x, m, jx, jy, jz);
	}
      }
    }
  }
}

void
fields_single_axpy_comp(fields_single_t *y, int ym, fields_single_real_t a, fields_single_t *x, int xm)
{
  for (int jz = y->ib[2]; jz < y->ib[2] + y->im[2]; jz++) {
    for (int jy = y->ib[1]; jy < y->ib[1] + y->im[1]; jy++) {
      for (int jx = y->ib[0]; jx < y->ib[0] + y->im[0]; jx++) {
	F3_S(y, ym, jx, jy, jz) += a * F3_S(x, xm, jx, jy, jz);
      }
    }
  }
}

void
fields_single_scale(fields_single_t *pf, fields_single_real_t val)
{
  for (int m = 0; m < pf->nr_comp; m++) {
    for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
      for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
	for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	  F3_S(pf, m, jx, jy, jz) *= val;
	}
      }
    }
  }
}

static void
_psc_mfields_single_zero_comp(mfields_single_t *flds, int m)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    fields_single_zero(psc_mfields_get_patch_single(flds, p), m);
  }
}

static void
_psc_mfields_single_axpy(mfields_single_t *yf, double alpha, mfields_single_t *xf)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    fields_single_axpy(psc_mfields_get_patch_single(yf, p), alpha, psc_mfields_get_patch_single(xf, p));
  }
}

static void
_psc_mfields_single_axpy_comp(mfields_single_t *yf, int ym, double alpha,
			      mfields_single_t *xf, int xm)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    fields_single_axpy_comp(psc_mfields_get_patch_single(yf, p), ym, alpha,
		       psc_mfields_get_patch_single(xf, p), xm);
  }
}

static void
_psc_mfields_single_scale(mfields_single_t *yf, double alpha)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    fields_single_scale(psc_mfields_get_patch_single(yf, p), alpha);
  }
}

static void
_psc_mfields_single_set_comp(mfields_single_t *yf, int m, double alpha)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    fields_single_set(psc_mfields_get_patch_single(yf, p), m, alpha);
  }
}

static void
_psc_mfields_single_copy_comp(mfields_single_t *to, int mto, mfields_single_t *from, int mfrom)
{
  for (int p = 0; p < to->nr_patches; p++) {
    fields_single_copy_comp(psc_mfields_get_patch_single(to, p), mto,
		       psc_mfields_get_patch_single(from, p), mfrom);
  }
}

// ======================================================================
// psc_mfields_single

static void
_psc_mfields_single_setup(mfields_single_t *flds)
{
  psc_mfields_setup_super(flds);

  struct mrc_patch *patches = mrc_domain_get_patches(flds->domain,
						     &flds->nr_patches);
  flds->flds = calloc(flds->nr_patches, sizeof(*flds->flds));
  for (int p = 0; p < flds->nr_patches; p++) {
    fields_single_t *pf = calloc(1, sizeof(*pf));
    flds->flds[p] = (struct psc_fields *) pf;
    int ilg[3] = { -flds->ibn[0], -flds->ibn[1], -flds->ibn[2] };
    int ihg[3] = { patches[p].ldims[0] + flds->ibn[0],
		   patches[p].ldims[1] + flds->ibn[1],
		   patches[p].ldims[2] + flds->ibn[2] };
    fields_single_alloc(pf, ilg, ihg, flds->nr_fields, flds->first_comp);
  }
}

static void
_psc_mfields_single_destroy(mfields_single_t *flds)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    fields_single_t *pf = psc_mfields_get_patch_single(flds, p);
    fields_single_free(pf);
    free(pf);
  }
  free(flds->flds);
}

// ======================================================================
// psc_mfields: subclass "single"
  
struct psc_mfields_ops psc_mfields_single_ops = {
  .name                  = "single",
  .setup                 = _psc_mfields_single_setup,
  .destroy               = _psc_mfields_single_destroy,
  .zero_comp             = _psc_mfields_single_zero_comp,
  .set_comp              = _psc_mfields_single_set_comp,
  .scale                 = _psc_mfields_single_scale,
  .copy_comp             = _psc_mfields_single_copy_comp,
  .axpy                  = _psc_mfields_single_axpy,
  .axpy_comp             = _psc_mfields_single_axpy_comp,
};

