
#include "psc.h"
#include "psc_fields_private.h"

#include <mrc_profile.h>
#include <mrc_params.h>

// ======================================================================
// psc_fields_size

unsigned int
psc_fields_size(struct psc_fields *pf)
{
  return pf->im[0] * pf->im[1] * pf->im[2];
}

// ======================================================================
// psc_fields_zero_comp

void
psc_fields_zero_comp(struct psc_fields *pf, int m)
{
  struct psc_fields_ops *ops = psc_fields_ops(pf);

  assert(ops && ops->zero_comp);
  ops->zero_comp(pf, m);
}

// ======================================================================
// psc_fields_set_comp

void
psc_fields_set_comp(struct psc_fields *pf, int m, double alpha)
{
  struct psc_fields_ops *ops = psc_fields_ops(pf);

  assert(ops && ops->set_comp);
  ops->set_comp(pf, m, alpha);
}

// ======================================================================
// psc_fields_scale_comp

void
psc_fields_scale_comp(struct psc_fields *pf, int m, double alpha)
{
  struct psc_fields_ops *ops = psc_fields_ops(pf);

  assert(ops && ops->scale_comp);
  ops->scale_comp(pf, m, alpha);
}

// ======================================================================
// psc_fields_copy_comp

void
psc_fields_copy_comp(struct psc_fields *to, int mto,
		     struct psc_fields *from, int mfrom)
{
  struct psc_fields_ops *ops = psc_fields_ops(to);
  assert(ops == psc_fields_ops(from));

  assert(ops && ops->copy_comp);
  ops->copy_comp(to, mto, from, mfrom);
}

// ======================================================================
// psc_fields_axpy_comp

void
psc_fields_axpy_comp(struct psc_fields *yf, int ym, double alpha,
		     struct psc_fields *xf, int xm)
{
  struct psc_fields_ops *ops = psc_fields_ops(yf);
  assert(ops == psc_fields_ops(xf));

  assert(ops && ops->axpy_comp);
  ops->axpy_comp(yf, ym, alpha, xf, xm);
}

// ======================================================================
// psc_fields_descr

#define VAR(x) (void *)offsetof(struct psc_fields, x)

static struct param psc_fields_descr[] = {
  { "ib"            , VAR(ib)                     , PARAM_INT3(0, 0, 0)  },
  { "im"            , VAR(im)                     , PARAM_INT3(0, 0, 0)  },
  { "nr_comp"       , VAR(nr_comp)                , PARAM_INT(1)         },
  { "first_comp"    , VAR(first_comp)             , PARAM_INT(0)         },
  { "p"             , VAR(p)                      , PARAM_INT(0)         },
  {}
};

// ======================================================================
// psc_fields_init

static void
psc_fields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_fields, &psc_fields_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_fields, &psc_fields_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_fields, &psc_fields_fortran_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_fields, &psc_fields_cuda_ops);
#endif
}

// ======================================================================
// psc_fields class

struct mrc_class_psc_fields mrc_class_psc_fields = {
  .name             = "psc_fields",
  .size             = sizeof(struct psc_fields),
  .param_descr      = psc_fields_descr,
  .init             = psc_fields_init,
};

