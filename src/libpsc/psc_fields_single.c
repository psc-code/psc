
#include "psc.h"
#include "psc_fields_single.h"
#include "psc_fields_c.h"

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

static void
psc_fields_single_zero_comp(struct psc_fields *pf, int m)
{
  memset(&F3_S(pf, m, pf->ib[0], pf->ib[1], pf->ib[2]), 0,
	 pf->im[0] * pf->im[1] * pf->im[2] * sizeof(fields_single_real_t));
}

static void
psc_fields_single_set_comp(struct psc_fields *pf, int m, double _val)
{
  fields_single_real_t val = _val;

  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_S(pf, m, jx, jy, jz) = val;
      }
    }
  }
}

static void
psc_fields_single_scale_comp(struct psc_fields *pf, int m, double _val)
{
  fields_single_real_t val = _val;
  for (int jz = pf->ib[2]; jz < pf->ib[2] + pf->im[2]; jz++) {
    for (int jy = pf->ib[1]; jy < pf->ib[1] + pf->im[1]; jy++) {
      for (int jx = pf->ib[0]; jx < pf->ib[0] + pf->im[0]; jx++) {
	F3_S(pf, m, jx, jy, jz) *= val;
      }
    }
  }
}

static void
psc_fields_single_copy_comp(struct psc_fields *pto, int m_to, struct psc_fields *pfrom, int m_from)
{
  for (int jz = pto->ib[2]; jz < pto->ib[2] + pto->im[2]; jz++) {
    for (int jy = pto->ib[1]; jy < pto->ib[1] + pto->im[1]; jy++) {
      for (int jx = pto->ib[0]; jx < pto->ib[0] + pto->im[0]; jx++) {
	F3_S(pto, m_to, jx, jy, jz) = F3_S(pfrom, m_from, jx, jy, jz);
      }
    }
  }
}

static void
psc_fields_single_axpy_comp(struct psc_fields *y, int ym, double _a, struct psc_fields *x, int xm)
{
  fields_single_real_t a = _a;

  for (int jz = y->ib[2]; jz < y->ib[2] + y->im[2]; jz++) {
    for (int jy = y->ib[1]; jy < y->ib[1] + y->im[1]; jy++) {
      for (int jx = y->ib[0]; jx < y->ib[0] + y->im[0]; jx++) {
	F3_S(y, ym, jx, jy, jz) += a * F3_S(x, xm, jx, jy, jz);
      }
    }
  }
}

// ======================================================================
// convert to c

static void
psc_fields_single_copy_from_c(struct psc_fields *flds_single, struct psc_fields *flds_c,
			      int mb, int me)
{
  for (int m = mb; m < me; m++) {
    for (int jz = flds_single->ib[2]; jz < flds_single->ib[2] + flds_single->im[2]; jz++) {
      for (int jy = flds_single->ib[1]; jy < flds_single->ib[1] + flds_single->im[1]; jy++) {
	for (int jx = flds_single->ib[0]; jx < flds_single->ib[0] + flds_single->im[0]; jx++) {
	  F3_S(flds_single, m, jx,jy,jz) = F3_C(flds_c, m, jx,jy,jz);
	}
      }
    }
  }
}

void
psc_fields_single_copy_to_c(struct psc_fields *flds_single, struct psc_fields *flds_c,
			    int mb, int me)
{
  for (int m = mb; m < me; m++) {
    for (int jz = flds_single->ib[2]; jz < flds_single->ib[2] + flds_single->im[2]; jz++) {
      for (int jy = flds_single->ib[1]; jy < flds_single->ib[1] + flds_single->im[1]; jy++) {
	for (int jx = flds_single->ib[0]; jx < flds_single->ib[0] + flds_single->im[0]; jx++) {
	  F3_C(flds_c, m, jx,jy,jz) = F3_S(flds_single, m, jx,jy,jz);
	}
      }
    }
  }
}

// ======================================================================
// psc_mfields: subclass "single"
  
struct psc_mfields_ops psc_mfields_single_ops = {
  .name                  = "single",
};

// ======================================================================
// psc_fields: subclass "single"
  
static struct mrc_obj_method psc_fields_single_methods[] = {
  MRC_OBJ_METHOD("copy_to_c",   psc_fields_single_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c", psc_fields_single_copy_from_c),
  {}
};

struct psc_fields_ops psc_fields_single_ops = {
  .name                  = "single",
  .methods               = psc_fields_single_methods,
  .setup                 = psc_fields_single_setup,
  .destroy               = psc_fields_single_destroy,
  .zero_comp             = psc_fields_single_zero_comp,
  .set_comp              = psc_fields_single_set_comp,
  .scale_comp            = psc_fields_single_scale_comp,
  .copy_comp             = psc_fields_single_copy_comp,
  .axpy_comp             = psc_fields_single_axpy_comp,
};

