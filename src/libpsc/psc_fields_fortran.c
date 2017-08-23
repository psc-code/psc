
#include "psc.h"
#include "psc_fields_as_fortran.h"
#include "psc_fields_c.h"

#include <mrc_params.h>
#include <stdlib.h>
#include <string.h>

#include "psc_fields_inc.h"

// ======================================================================
// convert to/from "c"

static void
psc_fields_fortran_copy_to_c(struct psc_fields *flds_fortran, struct psc_fields *flds_c,
			     int mb, int me)
{
  for (int m = mb; m < me; m++) {
    for (int jz = flds_c->ib[2]; jz < flds_c->ib[2] + flds_c->im[2]; jz++) {
      for (int jy = flds_c->ib[1]; jy < flds_c->ib[1] + flds_c->im[1]; jy++) {
	for (int jx = flds_c->ib[0]; jx < flds_c->ib[0] + flds_c->im[0]; jx++) {
	  F3_C(flds_c, m, jx,jy,jz) = F3_FORTRAN(flds_fortran, m, jx,jy,jz);
	}
      }
    }
  }
}

static void
psc_fields_fortran_copy_from_c(struct psc_fields *flds_fortran, struct psc_fields *flds_c,
			       int mb, int me)
{
  for (int m = mb; m < me; m++) {
    for (int jz = flds_c->ib[2]; jz < flds_c->ib[2] + flds_c->im[2]; jz++) {
      for (int jy = flds_c->ib[1]; jy < flds_c->ib[1] + flds_c->im[1]; jy++) {
	for (int jx = flds_c->ib[0]; jx < flds_c->ib[0] + flds_c->im[0]; jx++) {
	  F3_FORTRAN(flds_fortran, m, jx,jy,jz) = F3_C(flds_c, m, jx,jy,jz);
	}
      }
    }
  }
}

static struct mrc_obj_method psc_fields_fortran_methods[] = {
  MRC_OBJ_METHOD("copy_to_c",   psc_fields_fortran_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c", psc_fields_fortran_copy_from_c),
  {}
};

#include "psc_fields_common.c"

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

// ======================================================================
// psc_mfields: subclass "fortran"
  
struct psc_mfields_ops psc_mfields_fortran_ops = {
  .name                  = "fortran",
};

// ======================================================================
// psc_fields: subclass "fortran"
  
struct psc_fields_ops psc_fields_fortran_ops = {
  .name                  = "fortran",
  .methods               = psc_fields_fortran_methods,
  .setup                 = psc_fields_fortran_setup,
  .destroy               = psc_fields_fortran_destroy,
  .zero_comp             = psc_fields_fortran_zero_comp,
  .set_comp              = psc_fields_fortran_set_comp,
  .scale_comp            = psc_fields_fortran_scale_comp,
  .copy_comp             = psc_fields_fortran_copy_comp,
  .axpy_comp             = psc_fields_fortran_axpy_comp,
};

