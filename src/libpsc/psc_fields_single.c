
#include "psc.h"
#include "psc_fields_single.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// FIXME, very duplicated from psc_fields_c.c

static void
psc_fields_single_setup(struct psc_fields *pf)
{
  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    size *= pf->im[d];
  }
  pf->data = calloc(pf->nr_comp * size, sizeof(fields_single_real_t));
}

static void
psc_fields_single_destroy(struct psc_fields *pf)
{
  free(pf->data);
}

void
fields_single_zero(struct psc_fields *pf, int m)
{
  memset(&F3_S(pf, m, pf->ib[0], pf->ib[1], pf->ib[2]), 0,
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_single_real_t));
}

void
fields_single_set(struct psc_fields *pf, int m, fields_single_real_t val)
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
fields_single_copy_comp(struct psc_fields *pto, int m_to, struct psc_fields *pfrom, int m_from)
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
fields_single_axpy(struct psc_fields *y, fields_single_real_t a, struct psc_fields *x)
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
fields_single_axpy_comp(struct psc_fields *y, int ym, fields_single_real_t a, struct psc_fields *x, int xm)
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
fields_single_scale(struct psc_fields *pf, fields_single_real_t val)
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
// psc_mfields: subclass "single"
  
struct psc_mfields_ops psc_mfields_single_ops = {
  .name                  = "single",
  .zero_comp             = _psc_mfields_single_zero_comp,
  .set_comp              = _psc_mfields_single_set_comp,
  .scale                 = _psc_mfields_single_scale,
  .copy_comp             = _psc_mfields_single_copy_comp,
  .axpy                  = _psc_mfields_single_axpy,
  .axpy_comp             = _psc_mfields_single_axpy_comp,
};

// ======================================================================
// psc_fields: subclass "single"
  
struct psc_fields_ops psc_fields_single_ops = {
  .name                  = "single",
  .setup                 = psc_fields_single_setup,
  .destroy               = psc_fields_single_destroy,
};

