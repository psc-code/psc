
#include "psc_inject_private.h"
#include "fields.hxx"

#include <mrc_profile.h>
#include <stdlib.h>

// ======================================================================
// psc_inject

// ----------------------------------------------------------------------
// _psc_inject_setup

static void
_psc_inject_setup(struct psc_inject *inject)
{
  // set up necessary bits for calculating / averaging density moment

  psc_inject_setup_member_objs(inject);

  inject->mflds_n = psc_output_fields_item_create_mfields(inject->item_n);
  psc_mfields_set_name(inject->mflds_n, "mflds_n");
  psc_mfields_list_add(&psc_mfields_base_list, &inject->mflds_n);
}

// ----------------------------------------------------------------------
// _psc_inject_destroy

static void
_psc_inject_destroy(struct psc_inject *inject)
{
  psc_mfields_destroy(inject->mflds_n);
}

#if 0
// ----------------------------------------------------------------------
// debug_dump

static void
copy_to_mrc_fld(struct mrc_fld *m3, struct psc_mfields *mflds)
{
  mfields_t mf(mflds);
  psc_foreach_patch(ppsc, p) {
    Fields3d<fields_t> F(mf[p]);
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);
    mrc_fld_foreach(m3, ix,iy,iz, 0,0) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	MRC_M3(m3p, m, ix,iy,iz) = F(m, ix,iy,iz);
      }
    } mrc_fld_foreach_end;
    mrc_fld_patch_put(m3);
  }
}

static void _mrc_unused
debug_dump(struct mrc_io *io, struct psc_mfields *mflds)
{
  /* if (ppsc->timestep % debug_every_step != 0) { */
  /*   return; */
  /* } */

  struct mrc_fld *mrc_fld = mrc_domain_m3_create(ppsc->mrc_domain);
  mrc_fld_set_name(mrc_fld, psc_mfields_name(mflds));
  mrc_fld_set_param_int(mrc_fld, "nr_ghosts", 2);
  mrc_fld_set_param_int(mrc_fld, "nr_comps", mflds->nr_fields);
  mrc_fld_setup(mrc_fld);
  for (int m = 0; m < mflds->nr_fields; m++) {
    mrc_fld_set_comp_name(mrc_fld, m, psc_mfields_comp_name(mflds, m));
  }
  copy_to_mrc_fld(mrc_fld, mflds);
  mrc_fld_write(mrc_fld, io);
  mrc_fld_destroy(mrc_fld);
}
#endif

// ----------------------------------------------------------------------
// psc_inject_run
//
// FIXME mostly duplicated from psc_setup_particles

void
psc_inject_run(struct psc_inject *inject, struct psc_mparticles *mprts_base,
	       struct psc_mfields *mflds_base)
{
  static int pr;
  if (!pr) {
    pr = prof_register("inject_run", 1., 0, 0);
  }  

  if (!inject->do_inject ||
      ppsc->timestep % inject->every_step != 0) {
    return;
  }

  struct psc_inject_ops *ops = psc_inject_ops(inject);
  assert(ops && ops->run);

  prof_start(pr);
  ops->run(inject, mprts_base, mflds_base);
  prof_stop(pr);
}

// ----------------------------------------------------------------------
// psc_inject_init

extern struct psc_inject_ops psc_inject_ops_single;
extern struct psc_inject_ops psc_inject_ops_double;
extern struct psc_inject_ops psc_inject_ops_cuda;

static void
psc_inject_init(void)
{
  mrc_class_register_subclass(&mrc_class_psc_inject, &psc_inject_ops_single);
  mrc_class_register_subclass(&mrc_class_psc_inject, &psc_inject_ops_double);
#ifdef USE_CUDA
  mrc_class_register_subclass(&mrc_class_psc_inject, &psc_inject_ops_cuda);
#endif
}

// ----------------------------------------------------------------------
// psc_inject class

struct mrc_class_psc_inject_ : mrc_class_psc_inject {
  mrc_class_psc_inject_() {
    name             = "psc_inject";
    size             = sizeof(struct psc_inject);
    param_descr      = psc_inject_descr;
    init             = psc_inject_init;
    setup            = _psc_inject_setup;
    destroy          = _psc_inject_destroy;
  }
} mrc_class_psc_inject;

