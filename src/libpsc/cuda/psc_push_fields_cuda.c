
#include "psc_push_fields_private.h"
#include "psc_fields_cuda.h"
#include "psc.h"

EXTERN_C void cuda_push_fields_a_E_yz(int p, struct psc_fields *flds);
EXTERN_C void cuda_push_fields_a_H_yz(int p, struct psc_fields *flds);
EXTERN_C void cuda_push_fields_b_H_yz(int p, struct psc_fields *flds);
EXTERN_C void cuda_push_fields_b_E_yz(int p, struct psc_fields *flds);

// ----------------------------------------------------------------------
// E-field propagation E^(n)    , H^(n), j^(n) 
//                  -> E^(n+0.5), H^(n), j^(n)

static void
psc_push_fields_cuda_push_a_E(struct psc_push_fields *push, struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", JXI, HX + 3);

  if (ppsc->domain.gdims[0] == 1) {
    cuda_push_fields_a_E_yz(flds->p, flds);
  } else {
    assert(0);
  }

  psc_fields_put_as(flds, flds_base, EX, EX + 3);
}

// ----------------------------------------------------------------------
// B-field propagation E^(n+0.5), H^(n    ), j^(n), m^(n+0.5)
//                  -> E^(n+0.5), H^(n+0.5), j^(n), m^(n+0.5)

static void
psc_push_fields_cuda_push_a_H(struct psc_push_fields *push, struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", EX, HX + 3);

  if (ppsc->domain.gdims[0] == 1) {
    cuda_push_fields_a_H_yz(flds->p, flds);
  } else {
    assert(0);
  }

  psc_fields_put_as(flds, flds_base, HX, HX + 3);
}

// ----------------------------------------------------------------------
// B-field propagation E^(n+0.5), B^(n+0.5), j^(n+1.0), m^(n+0.5)
//                  -> E^(n+0.5), B^(n+1.0), j^(n+1.0), m^(n+0.5)

static void
psc_push_fields_cuda_push_b_H(struct psc_push_fields *push, struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", EX, HX + 3);

  if (ppsc->domain.gdims[0] == 1) {
    cuda_push_fields_b_H_yz(flds->p, flds);
  } else {
    assert(0);
  }

  psc_fields_put_as(flds, flds_base, HX, HX + 3);
}

// ----------------------------------------------------------------------
// E-field propagation E^(n+0.5), B^(n+1.0), j^(n+1.0) 
//                  -> E^(n+1.0), B^(n+1.0), j^(n+1.0)

static void
psc_push_fields_cuda_push_b_E(struct psc_push_fields *push, struct psc_fields *flds_base)
{
  struct psc_fields *flds = psc_fields_get_as(flds_base, "cuda", JXI, HX + 3);

  if (ppsc->domain.gdims[0] == 1) {
    cuda_push_fields_b_E_yz(flds->p, flds);
  } else {
    assert(0);
  }

  psc_fields_put_as(flds, flds_base, EX, EX + 3);
}

// ======================================================================
// psc_push_fields: subclass "cuda"

struct psc_push_fields_ops psc_push_fields_cuda_ops = {
  .name                  = "cuda",
  .push_a_E              = psc_push_fields_cuda_push_a_E,
  .push_a_H              = psc_push_fields_cuda_push_a_H,
  .push_b_H              = psc_push_fields_cuda_push_b_H,
  .push_b_E              = psc_push_fields_cuda_push_b_E,
};
