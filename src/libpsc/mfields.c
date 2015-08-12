
#include "psc.h"
#include "psc_fields_as_c.h"

#include <mrc_params.h>
#include <mrc_profile.h>
#include <mrc_io.h>
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
    struct psc_fields *pf = psc_fields_create(MPI_COMM_NULL);
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

  for (int m = 0; m < flds->nr_fields; m++) {
    free(flds->comp_name[m]);
  }
  free(flds->comp_name);
}

static void
_psc_mfields_write(struct psc_mfields *mflds, struct mrc_io *io)
{
  mrc_io_write_ref(io, mflds, "domain", mflds->domain);
  for (int m = 0; m < mflds->nr_fields; m++) {
    char name[20]; sprintf(name, "comp_name_%d", m);
    mrc_io_write_string(io, mflds, name, psc_mfields_comp_name(mflds, m));
  }

  for (int p = 0; p < mflds->nr_patches; p++) {
    char name[20]; sprintf(name, "flds%d", p);
    mrc_io_write_ref(io, mflds, name, mflds->flds[p]);
  }
}

static void
_psc_mfields_read(struct psc_mfields *mflds, struct mrc_io *io)
{
  mflds->domain = mrc_io_read_ref(io, mflds, "domain", mrc_domain);
  mrc_domain_get_patches(mflds->domain, &mflds->nr_patches);

  mflds->comp_name = calloc(mflds->nr_fields, sizeof(*mflds->comp_name));
  for (int m = 0; m < mflds->nr_fields; m++) {
    char name[20]; sprintf(name, "comp_name_%d", m);
    char *s;
    mrc_io_read_string(io, mflds, name, &s);
    if (s) {
      psc_mfields_set_comp_name(mflds, m, s);
    }
  }

  mflds->flds = calloc(mflds->nr_patches, sizeof(*mflds->flds));
  mprintf("nr_p %d\n", mflds->nr_patches);
  for (int p = 0; p < mflds->nr_patches; p++) {
    char name[20]; sprintf(name, "flds%d", p);
    mflds->flds[p] = mrc_io_read_ref_comm(io, mflds, name, psc_fields,
					  MPI_COMM_NULL);
  }
  // FIXME mark as set up?
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

static void
copy_to_mrc_fld(struct mrc_fld *m3, struct psc_mfields *mflds_base)
{
  struct psc_mfields *mflds = 
    psc_mfields_get_as(mflds_base, FIELDS_TYPE, 0, mflds_base->nr_fields);
  psc_foreach_patch(ppsc, p) {
    struct psc_fields *pf = psc_mfields_get_patch(mflds, p);
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);
    mrc_fld_foreach(m3, ix,iy,iz, 0,0) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	MRC_M3(m3p,m, ix,iy,iz) = F3(pf,m, ix,iy,iz);
      }
    } mrc_fld_foreach_end;
    mrc_fld_patch_put(m3);
  }

  psc_mfields_put_as(mflds, mflds_base, 0, 0);
}

void
psc_mfields_write_as_mrc_fld(struct psc_mfields *mflds, struct mrc_io *io)
{
  struct mrc_fld *fld = mrc_domain_m3_create(ppsc->mrc_domain);
  mrc_fld_set_name(fld, psc_mfields_name(mflds));
  mrc_fld_set_param_int(fld, "nr_ghosts", 2);
  mrc_fld_set_param_int(fld, "nr_comps", mflds->nr_fields);
  mrc_fld_setup(fld);
  for (int m = 0; m < mrc_fld_nr_comps(fld); m++) {
      mrc_fld_set_comp_name(fld, m, psc_mfields_comp_name(mflds, m));
  }
  copy_to_mrc_fld(fld, mflds);
  mrc_fld_write(fld, io);
  mrc_fld_destroy(fld);
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
    char s[strlen(type) + strlen(type_base) + 20]; sprintf(s, "copy_to_%s", type);
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
#ifdef USE_CUDA2
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_cuda2_ops);
#endif
#ifdef USE_ACC
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_acc_ops);
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
  .read             = _psc_mfields_read,
  .write            = _psc_mfields_write,
};

