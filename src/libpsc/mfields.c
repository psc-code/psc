
#include "psc.h"

#include <mrc_params.h>
#include <mrc_profile.h>
#include <stdlib.h>
#include <string.h>

// ======================================================================
// _psc_mfields_setup

static void
_psc_mfields_setup(struct psc_mfields *flds)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(flds);

  flds->comp_name = calloc(flds->nr_fields, sizeof(*flds->comp_name));

  struct mrc_patch *patches = mrc_domain_get_patches(flds->domain,
						     &flds->nr_patches);
  flds->flds = calloc(flds->nr_patches, sizeof(*flds->flds));
  for (int p = 0; p < flds->nr_patches; p++) {
    struct psc_fields *pf = psc_fields_create(psc_mfields_comm(flds));
    psc_fields_set_type(pf, ops->name);
    char name[20]; sprintf(name, "flds%d", p);
    psc_fields_set_name(pf, name);
    for (int d = 0; d < 3; d++) {
      pf->ib[d] = -flds->ibn[d];
      pf->im[d] = patches[p].ldims[d] + 2 * flds->ibn[d];
    }
    pf->nr_comp = flds->nr_fields;
    pf->first_comp = flds->first_comp;
    pf->p = p;
    psc_fields_setup(pf);
    flds->flds[p] = pf;
  }
}

// ======================================================================
// psc_mfields_destroy

static void
_psc_mfields_destroy(struct psc_mfields *flds)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    struct psc_fields *pf = flds->flds[p];
    psc_fields_destroy(pf);
  }
  free(flds->flds);

  // sub-destroy has already been called
  for (int m = 0; m < flds->nr_fields; m++) {
    free(flds->comp_name[m]);
  }
  free(flds->comp_name);
}

void
psc_mfields_set_domain(struct psc_mfields *flds, struct mrc_domain *domain)
{
  flds->domain = domain;
}

void
psc_mfields_set_comp_name(struct psc_mfields *flds, int m, const char *s)
{
  assert(m >= flds->first_comp && m < flds->first_comp + flds->nr_fields);
  flds->comp_name[m - flds->first_comp] = strdup(s);
}

const char *
psc_mfields_comp_name(struct psc_mfields *flds, int m)
{
  assert(m >= flds->first_comp && m < flds->first_comp + flds->nr_fields);
  return flds->comp_name[m - flds->first_comp];
}

void
psc_mfields_zero_comp(struct psc_mfields *flds, int m)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    psc_fields_zero_comp(psc_mfields_get_patch(flds, p), m);
  }
}

void
psc_mfields_zero_range(struct psc_mfields *flds, int mb, int me)
{
  for (int m = mb; m < me; m++) {
    psc_mfields_zero_comp(flds, m);
  }
}

void
psc_mfields_set_comp(struct psc_mfields *flds, int m, double alpha)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    psc_fields_set_comp(psc_mfields_get_patch(flds, p), m, alpha);
  }
}

void
psc_mfields_scale_comp(struct psc_mfields *flds, int m, double alpha)
{
  for (int p = 0; p < flds->nr_patches; p++) {
    psc_fields_scale_comp(psc_mfields_get_patch(flds, p), m, alpha);
  }
}

void
psc_mfields_scale(struct psc_mfields *flds, double alpha)
{
  for (int m = flds->first_comp; m < flds->first_comp + flds->nr_fields; m++) {
    psc_mfields_scale_comp(flds, m, alpha);
  }
}

void
psc_mfields_copy_comp(struct psc_mfields *to, int mto,
		      struct psc_mfields *from, int mfrom)
{
  for (int p = 0; p < to->nr_patches; p++) {
    psc_fields_copy_comp(psc_mfields_get_patch(to, p), mto,
			 psc_mfields_get_patch(from, p), mfrom);
  }
}

void
psc_mfields_axpy_comp(struct psc_mfields *yf, int ym, double alpha,
		      struct psc_mfields *xf, int xm)
{
  for (int p = 0; p < yf->nr_patches; p++) {
    psc_fields_axpy_comp(psc_mfields_get_patch(yf, p), ym, alpha,
			 psc_mfields_get_patch(xf, p), xm);
  }
}

void
psc_mfields_axpy(struct psc_mfields *yf, double alpha,
		 struct psc_mfields *xf)
{
  for (int m = yf->first_comp; m < yf->first_comp + yf->nr_fields; m++) {
    psc_mfields_axpy_comp(yf, m, alpha, xf, m);
  }
}

struct psc_mfields *
psc_mfields_get_as(struct psc_mfields *mflds_base, const char *type,
		   int mb, int me)
{
  const char *type_base = psc_mfields_type(mflds_base);
  // If we're already the subtype, nothing to be done
  if (strcmp(type_base, type) == 0)
    return mflds_base;

  static int pr;
  if (!pr) {
    pr = prof_register("mfields_get_as", 1., 0, 0);
  }
  prof_start(pr);

  struct psc_mfields *mflds = psc_mfields_create(psc_mfields_comm(mflds_base));
  psc_mfields_set_type(mflds, type);
  psc_mfields_set_domain(mflds, mflds_base->domain);
  psc_mfields_set_param_int(mflds, "nr_fields", mflds_base->nr_fields);
  psc_mfields_set_param_int3(mflds, "ibn", mflds_base->ibn);
  psc_mfields_set_param_int(mflds, "first_comp", mflds_base->first_comp);
  psc_mfields_setup(mflds);

  for (int p = 0; p < mflds_base->nr_patches; p++) {
    struct psc_fields *flds_base = psc_mfields_get_patch(mflds_base, p);
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
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
  }

  prof_stop(pr);
  return mflds;
}

void
psc_mfields_put_as(struct psc_mfields *mflds, struct psc_mfields *mflds_base,
		   int mb, int me)
{
  // If we're already the subtype, nothing to be done
  const char *type = psc_mfields_type(mflds);
  const char *type_base = psc_mfields_type(mflds_base);
  if (strcmp(type_base, type) == 0)
    return;

  static int pr;
  if (!pr) {
    pr = prof_register("mfields_put_as", 1., 0, 0);
  }
  prof_start(pr);

  for (int p = 0; p < mflds_base->nr_patches; p++) {
    struct psc_fields *flds_base = psc_mfields_get_patch(mflds_base, p);
    struct psc_fields *flds = psc_mfields_get_patch(mflds, p);
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
  }
  psc_mfields_destroy(mflds);

  prof_stop(pr);
}

#define MAKE_MFIELDS_GET_PUT(type)					\
									\
struct psc_mfields *							\
psc_mfields_get_##type(struct psc_mfields *mflds_base, int mb, int me)	\
{									\
  return psc_mfields_get_as(mflds_base, #type, mb, me);			\
}									\
									\
void									\
psc_mfields_put_##type(struct psc_mfields *mflds,			\
		       struct psc_mfields *mflds_base, int mb, int me)	\
{									\
  psc_mfields_put_as(mflds, mflds_base, mb, me);			\
}									\

MAKE_MFIELDS_GET_PUT(c)
MAKE_MFIELDS_GET_PUT(single)
MAKE_MFIELDS_GET_PUT(fortran)
#ifdef USE_CUDA
MAKE_MFIELDS_GET_PUT(cuda)
#endif

// ======================================================================

void
psc_mfields_list_add(list_t *head, struct psc_mfields **flds_p)
{
  struct psc_mfields_list_entry *p = malloc(sizeof(*p));
  p->flds_p = flds_p;
  list_add_tail(&p->entry, head);
}

void
psc_mfields_list_del(list_t *head, struct psc_mfields **flds_p)
{
  struct psc_mfields_list_entry *p;
  __list_for_each_entry(p, head, entry, struct psc_mfields_list_entry) {
    if (p->flds_p == flds_p) {
      list_del(&p->entry);
      free(p);
      return;
    }
  }
  assert(0);
}

// ======================================================================

static void
psc_mfields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_single_ops);
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_fortran_ops);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_mix_ops);
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_cuda_ops);
#endif
}

#define VAR(x) (void *)offsetof(struct psc_mfields, x)
static struct param psc_mfields_descr[] = {
  { "nr_fields"      , VAR(nr_fields)       , PARAM_INT(1)        },
  { "ibn"            , VAR(ibn)             , PARAM_INT3(0, 0, 0) },
  { "first_comp"     , VAR(first_comp)      , PARAM_INT(0)        },
  {},
};
#undef VAR

struct mrc_class_psc_mfields mrc_class_psc_mfields = {
  .name             = "psc_mfields",
  .size             = sizeof(struct psc_mfields),
  .init             = psc_mfields_init,
  .param_descr      = psc_mfields_descr,
  .setup            = _psc_mfields_setup,
  .destroy          = _psc_mfields_destroy,
};

