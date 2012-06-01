
#include "psc_push_fields_private.h"
#include "psc_fields_cuda.h"
#include "psc.h"
#include "psc_bnd.h"
#include "psc_bnd_fields.h"

EXTERN_C void cuda_push_fields_E_yz(int p, struct psc_fields *flds);
EXTERN_C void cuda_push_fields_H_yz(int p, struct psc_fields *flds);

// ----------------------------------------------------------------------
// psc_push_fields_cuda_push_E

static void
psc_push_fields_cuda_push_E(struct psc_push_fields *push, struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", JXI, HX + 3);
  if (ppsc->domain.gdims[0] == 1) {
    cuda_push_fields_E_yz(flds->p, flds);
  } else {
    assert(0);
  }
  psc_fields_put_as(flds, flds_base, EX, EX + 3);
}

// ----------------------------------------------------------------------
// psc_push_fields_cuda_push_H

static void
psc_push_fields_cuda_push_H(struct psc_push_fields *push, struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", EX, HX + 3);
  if (ppsc->domain.gdims[0] == 1) {
    cuda_push_fields_H_yz(flds->p, flds);
  } else {
    assert(0);
  }
  psc_fields_put_as(flds, flds_base, HX, HX + 3);
}

// ======================================================================
// psc_push_fields: subclass "cuda"

struct psc_push_fields_ops psc_push_fields_cuda_ops = {
  .name                  = "cuda",
  .push_E                = psc_push_fields_cuda_push_E,
  .push_H                = psc_push_fields_cuda_push_H,
};
