
#include "psc.h"
#include "fields.hxx"
#include "psc_fields_as_c.h"

#include <mrc_params.h>
#include <mrc_profile.h>
#include <mrc_io.h>
#include <stdlib.h>
#include <string.h>

std::list<MfieldsBase*> MfieldsBase::instances;

// ======================================================================
// _psc_mfields_setup

static void
_psc_mfields_setup(struct psc_mfields *mflds)
{
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
  for (int m = 0; m < mflds->nr_fields; m++) {
    char name[20]; sprintf(name, "comp_name_%d", m);
    mrc_io_write_string(io, mflds, name, psc_mfields_comp_name(mflds, m));
  }
}

static void
_psc_mfields_read(struct psc_mfields *mflds, struct mrc_io *io)
{
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
  assert(m >= 0 && m < flds->nr_fields);
  flds->comp_name[m] = strdup(s);
}

const char *
psc_mfields_comp_name(struct psc_mfields *flds, int m)
{
  assert(m >= 0 && m < flds->nr_fields);
  return flds->comp_name[m];
}

// ----------------------------------------------------------------------

static void
copy_to_mrc_fld(struct mrc_fld *m3, struct psc_mfields *_mflds_base)
{
  using Fields = Fields3d<mfields_t::fields_t>;

  auto mflds_base = PscMfieldsBase{_mflds_base};
  mfields_t mf = mflds_base.get_as<mfields_t>(0, mflds_base->n_comps());

  psc_foreach_patch(ppsc, p) {
    Fields F(mf[p]);
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);
    mrc_fld_foreach(m3, ix,iy,iz, 0,0) {
      for (int m = 0; m < mflds_base->n_comps(); m++) {
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
  struct mrc_fld *fld = mrc_domain_m3_create(ppsc->mrc_domain_);
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

inline void MfieldsBase::convert(MfieldsBase& mf_from, MfieldsBase& mf_to, int mb, int me)
{
  // FIXME, implementing == wouldn't hurt
  assert(&mf_from.grid() == &mf_to.grid());
  
  auto convert_to = mf_from.convert_to().find(std::type_index(typeid(mf_to)));
  if (convert_to != mf_from.convert_to().cend()) {
    convert_to->second(mf_from, mf_to, mb, me);
    return;
  }
  
  auto convert_from = mf_to.convert_from().find(std::type_index(typeid(mf_from)));
  if (convert_from != mf_to.convert_from().cend()) {
    convert_from->second(mf_to, mf_from, mb, me);
    return;
  }

  fprintf(stderr, "ERROR: no conversion known from %s to %s!\n",
	  typeid(mf_from).name(), typeid(mf_to).name());
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
  { "nr_fields"      , VAR(nr_fields)       , PARAM_INT(1)        },
  { "ibn"            , VAR(ibn)             , PARAM_INT3(0, 0, 0) },
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

