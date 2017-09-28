
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
psc_mfields_zero_comp(struct psc_mfields *mflds, int m)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(mflds);

  assert(ops && ops->zero_comp);
  ops->zero_comp(mflds, m);
}

void
psc_mfields_set_comp(struct psc_mfields *mflds, int m, double alpha)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(mflds);

  assert(ops && ops->set_comp);
  ops->set_comp(mflds, m, alpha);
}

void
psc_mfields_scale_comp(struct psc_mfields *mflds, int m, double alpha)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(mflds);

  assert(ops && ops->scale_comp);
  ops->scale_comp(mflds, m, alpha);
}

void
psc_mfields_copy_comp(struct psc_mfields *to, int mto,
		      struct psc_mfields *from, int mfrom)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(to);
  assert(ops == psc_mfields_ops(from));

  assert(ops && ops->copy_comp);
  ops->copy_comp(to, mto, from, mfrom);
}

void
psc_mfields_axpy_comp(struct psc_mfields *yf, int ym, double alpha,
		      struct psc_mfields *xf, int xm)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(yf);
  assert(ops == psc_mfields_ops(xf));

  assert(ops && ops->axpy_comp);
  ops->axpy_comp(yf, ym, alpha, xf, xm);
}

void
psc_mfields_zero_range(struct psc_mfields *mflds, int mb, int me)
{
  for (int m = mb; m < me; m++) {
    psc_mfields_zero_comp(mflds, m);
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
    fields_t flds = fields_t_mflds(mflds, p);
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);
    mrc_fld_foreach(m3, ix,iy,iz, 0,0) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	MRC_M3(m3p ,m, ix,iy,iz) = _F3(flds, m, ix,iy,iz);
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

static void
copy(struct psc_mfields *mflds_from, struct psc_mfields *mflds_to,
     const char *type_from, const char *type_to, int mb, int me)
{
  psc_mfields_copy_func_t copy_to, copy_from;

  char s[strlen(type_to) + 12];
  sprintf(s, "copy_to_%s", type_to);
  copy_to = (psc_mfields_copy_func_t) psc_mfields_get_method(mflds_from, s);
  if (!copy_to) {
    sprintf(s, "copy_from_%s", type_from);
    copy_from = (psc_mfields_copy_func_t) psc_mfields_get_method(mflds_to, s);
  }
  if (!copy_to && !copy_from) {
    fprintf(stderr, "ERROR: no 'copy_to_%s' in psc_mfields '%s' and "
	    "no 'copy_from_%s' in '%s'!\n",
	    type_to, psc_mfields_type(mflds_from), type_from, psc_mfields_type(mflds_to));
    assert(0);
  }

  if (copy_to) {
    copy_to(mflds_from, mflds_to, mb, me);
  } else {
    copy_from(mflds_to, mflds_from, mb, me);
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

  /* static int cnt; */
  /* mprintf("get_as %s (%s) %d %d\n", type, psc_mfields_type(mflds_base), mb, me); */
  /* if (cnt++ == 10) assert(0); */
  
  struct psc_mfields *mflds = psc_mfields_create(psc_mfields_comm(mflds_base));
  psc_mfields_set_type(mflds, type);
  psc_mfields_set_domain(mflds, mflds_base->domain);
  psc_mfields_set_param_int(mflds, "nr_fields", mflds_base->nr_fields);
  psc_mfields_set_param_int3(mflds, "ibn", mflds_base->ibn);
  psc_mfields_set_param_int(mflds, "first_comp", mflds_base->first_comp);
  psc_mfields_setup(mflds);

  copy(mflds_base, mflds, type_base, type, mb, me);

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

  copy(mflds, mflds_base, type, type_base, mb, me);
  psc_mfields_destroy(mflds);

  prof_stop(pr);
}

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

