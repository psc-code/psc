
#include "psc.h"
#include "fields.hxx"
#include "psc_fields_as_c.h"

#include <mrc_params.h>
#include <mrc_profile.h>
#include <mrc_io.h>
#include <stdlib.h>
#include <string.h>

using Fields = Fields3d<fields_t>;

// ======================================================================
// _psc_mfields_setup

static void
_psc_mfields_setup(struct psc_mfields *mflds)
{
  assert(mflds->domain);
  
  mflds->comp_name = new char* [mflds->nr_fields]();
}

// ======================================================================
// psc_mfields_destroy

static void
_psc_mfields_destroy(struct psc_mfields *mflds)
{
  if (mflds->comp_name) {
    for (int m = 0; m < mflds->nr_fields; m++) {
      free(mflds->comp_name[m]);
    }
    delete[] mflds->comp_name;
  }
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

  mflds->comp_name = new char* [mflds->nr_fields];
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
psc_mfields_clear_rhof(struct psc_mfields *mflds)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(mflds);

  assert(ops && ops->clear_rhof);
  ops->clear_rhof(mflds);
}

void
psc_mfields_accumulate_rho_p(struct psc_mfields *mflds,
			     struct psc_mparticles *mprts)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(mflds);

  assert(ops && ops->accumulate_rho_p);
  ops->accumulate_rho_p(mflds, mprts);
}

void
psc_mfields_synchronize_rho(struct psc_mfields *mflds)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(mflds);

  assert(ops && ops->synchronize_rho);
  ops->synchronize_rho(mflds);
}

void
psc_mfields_compute_rhob(struct psc_mfields *mflds)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(mflds);

  assert(ops && ops->compute_rhob);
  ops->compute_rhob(mflds);
}

void
psc_mfields_compute_curl_b(struct psc_mfields *mflds)
{
  struct psc_mfields_ops *ops = psc_mfields_ops(mflds);

  assert(ops && ops->compute_curl_b);
  ops->compute_curl_b(mflds);
}

// ----------------------------------------------------------------------

static void
copy_to_mrc_fld(struct mrc_fld *m3, struct psc_mfields *mflds_base)
{
  mfields_t mf = mflds_base->get_as<mfields_t>(0, mflds_base->nr_fields);

  psc_foreach_patch(ppsc, p) {
    Fields F(mf[p]);
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);
    mrc_fld_foreach(m3, ix,iy,iz, 0,0) {
      for (int m = 0; m < mflds_base->nr_fields; m++) {
	MRC_M3(m3p ,m, ix,iy,iz) = F(m, ix,iy,iz);
      }
    } mrc_fld_foreach_end;
    mrc_fld_patch_put(m3);
  }

  mf.put_as(mflds_base, 0, 0);
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
  psc_mfields_set_param_obj(mflds, "domain", mflds_base->domain);
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
  struct psc_mfields_list_entry *p = new struct psc_mfields_list_entry();
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
      delete p;
      return;
    }
  }
  assert(0);
}

// ======================================================================

extern struct psc_mfields_ops psc_mfields_c_ops;
extern struct psc_mfields_ops psc_mfields_single_ops;
extern struct psc_mfields_ops psc_mfields_vpic_ops;
extern struct psc_mfields_ops psc_mfields_cuda_ops;

static void
psc_mfields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_c_ops);
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_single_ops);
#ifdef USE_VPIC
  mrc_class_register_subclass(&mrc_class_psc_mfields, &psc_mfields_vpic_ops);
#endif
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
  { "domain"         , VAR(domain)          , PARAM_OBJ(mrc_domain) },
  { "nr_fields"      , VAR(nr_fields)       , PARAM_INT(1)        },
  { "ibn"            , VAR(ibn)             , PARAM_INT3(0, 0, 0) },
  { "first_comp"     , VAR(first_comp)      , PARAM_INT(0)        },

  { "nr_patches"     , VAR(nr_patches)      , MRC_VAR_INT         },
  {},
};
#undef VAR

struct mrc_class_psc_mfields_ : mrc_class_psc_mfields {
  mrc_class_psc_mfields_() {
    name             = "psc_mfields";
    size             = sizeof(struct psc_mfields);
    init             = psc_mfields_init;
    param_descr      = psc_mfields_descr;
    setup            = _psc_mfields_setup;
    destroy          = _psc_mfields_destroy;
    read             = _psc_mfields_read;
    write            = _psc_mfields_write;
  }
} mrc_class_psc_mfields;

