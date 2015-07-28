
#include "psc.h"
#include "psc_fields_private.h"

#include <mrc_profile.h>
#include <mrc_params.h>
#include <string.h>

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
// psc_fields_zero_range

void
psc_fields_zero_range(struct psc_fields *pf, int mb, int me)
{
  for (int m = mb; m < me; m++) {
    psc_fields_zero_comp(pf, m);
  }
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
// psc_fields_get_as

struct psc_fields *
psc_fields_get_as(struct psc_fields *flds_base, const char *type,
		  int mb, int me)
{
  const char *type_base = psc_fields_type(flds_base);
  // If we're already the subtype, nothing to be done
  if (strcmp(type_base, type) == 0)
    return flds_base;

  static int pr;
  if (!pr) {
    pr = prof_register("fields_get_as", 1., 0, 0);
  }
  prof_start(pr);

  struct psc_fields *flds = psc_fields_create(psc_fields_comm(flds_base));
  psc_fields_set_type(flds, type);
  psc_fields_set_param_int3(flds, "ib", flds_base->ib);
  psc_fields_set_param_int3(flds, "im", flds_base->im);
  psc_fields_set_param_int(flds, "nr_comp", flds_base->nr_comp);
  psc_fields_set_param_int(flds, "first_comp", flds_base->first_comp);
  psc_fields_set_param_int(flds, "p", flds_base->p);
  psc_fields_setup(flds);

  char s[strlen(type) + 12]; sprintf(s, "copy_to_%s", type);
  psc_fields_copy_to_func_t copy_to = (psc_fields_copy_to_func_t)
    psc_fields_get_method(flds_base, s);
  if (copy_to) {
    copy_to(flds_base, flds, mb, me);
  } else {
    sprintf(s, "copy_from_%s", type_base);
    psc_fields_copy_to_func_t copy_from = (psc_fields_copy_from_func_t)
      psc_fields_get_method(flds, s);
    if (copy_from) {
      copy_from(flds, flds_base, mb, me);
    } else {
      fprintf(stderr, "ERROR: no 'copy_to_%s' in psc_fields '%s' and "
	      "no 'copy_from_%s' in '%s'!\n",
	      type, psc_fields_type(flds_base), type_base, psc_fields_type(flds));
      assert(0);
    }
  }

  prof_stop(pr);
  return flds;
}

// ======================================================================
// psc_fields_put_as

void
psc_fields_put_as(struct psc_fields *flds, struct psc_fields *flds_base,
		  int mb, int me)
{
  // If we're already the subtype, nothing to be done
  const char *type = psc_fields_type(flds);
  const char *type_base = psc_fields_type(flds_base);
  if (strcmp(type_base, type) == 0)
    return;

  static int pr;
  if (!pr) {
    pr = prof_register("fields_put_as", 1., 0, 0);
  }
  prof_start(pr);

  char s[strlen(type) + 12]; sprintf(s, "copy_from_%s", type);
  psc_fields_copy_from_func_t copy_from = (psc_fields_copy_from_func_t)
    psc_fields_get_method(flds_base, s);
  if (copy_from) {
    copy_from(flds_base, flds, mb, me);
  } else {
    sprintf(s, "copy_to_%s", type_base);
    psc_fields_copy_from_func_t copy_to = (psc_fields_copy_from_func_t)
      psc_fields_get_method(flds, s);
    if (copy_to) {
      copy_to(flds, flds_base, mb, me);
    } else {
      fprintf(stderr, "ERROR: no 'copy_from_%s' in psc_fields '%s' and "
	      "no 'copy_to_%s' in '%s'!\n",
	      type, psc_fields_type(flds_base), type_base, psc_fields_type(flds));
      assert(0);
    }
  }
  psc_fields_destroy(flds);

  prof_stop(pr);
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
#ifdef USE_CUDA2
  mrc_class_register_subclass(&mrc_class_psc_fields, &psc_fields_cuda2_ops);
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

