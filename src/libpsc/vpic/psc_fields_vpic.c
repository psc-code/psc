
#include "psc_fields_vpic.h"

#include "psc_fields_c.h"

#include "mrc_domain.h"

#include "vpic_iface.h"

// ----------------------------------------------------------------------
// psc_mfields_vpic_setup

static void
psc_mfields_vpic_setup(struct psc_mfields *mflds)
{
  struct psc_mfields_vpic *sub = psc_mfields_vpic(mflds);

  psc_mfields_setup_super(mflds);

  mrc_domain_get_patches(mflds->domain,
			 &mflds->nr_patches);
  assert(mflds->nr_patches == 1);

  /* struct vpic_info info = { */
  /*   /\* .dx = dx[0], *\/ */
  /*   /\* .dy = dx[1], *\/ */
  /*   /\* .dz = dx[2], *\/ */
  /*   /\* .dt = , *\/ */
  /*   .c = 1., */
  /*   .eps0 = 1., */
  /* }; */

  /* vpic_base_init(&info); */

  sub->vmflds = vpic_mfields_create();
  //vpic_mfields_ctor(sub->vmflds);
}

// ======================================================================
// convert from/to "c"

static void
psc_mfields_vpic_copy_from_c(struct psc_mfields *mflds_vpic, struct psc_mfields *mflds_c,
			    int mb, int me)
{
  for (int p = 0; p < mflds_vpic->nr_patches; p++) {
#if 0
    fields_c_t flds_c = fields_c_t_mflds(mflds_c, p);
    fields_vpic_t flds_vpic = fields_vpic_t_mflds(mflds_vpic, p);

    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
	for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
	  for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	    _F3_VPIC(flds_vpic, m, jx,jy,jz) = _F3_C(flds_c, m, jx,jy,jz);
	  }
	}
      }
    }
#endif
  }
}

static void
psc_mfields_vpic_copy_to_c(struct psc_mfields *mflds_vpic, struct psc_mfields *mflds_c,
			  int mb, int me)
{
  for (int p = 0; p < mflds_vpic->nr_patches; p++) {
#if 0
    fields_c_t flds_c = fields_c_t_mflds(mflds_c, p);
    fields_vpic_t flds_vpic = fields_vpic_t_mflds(mflds_vpic, p);
  
    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
	for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
	  for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	    _F3_C(flds_c, m, jx,jy,jz) = _F3_VPIC(flds_vpic, m, jx,jy,jz);
	  }
	}
      }
    }
#endif
  }
}

// ======================================================================
// psc_mfields: subclass "vpic"

static struct mrc_obj_method psc_mfields_vpic_methods[] = {
  MRC_OBJ_METHOD("copy_to_c"       , psc_mfields_vpic_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c"     , psc_mfields_vpic_copy_from_c),
#if 0
  MRC_OBJ_METHOD("copy_to_single"  , psc_mfields_vpic_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_mfields_vpic_copy_from_single),
#endif
  {}
};

struct psc_mfields_ops psc_mfields_vpic_ops = {
  .name                  = "vpic",
  .size                  = sizeof(struct psc_mfields_vpic),
  .methods               = psc_mfields_vpic_methods,
  .setup                 = psc_mfields_vpic_setup,
#if 0
  .destroy               = psc_mfields_vpic_destroy,
#ifdef HAVE_LIBHDF5_HL
  .write                 = psc_mfields_vpic_write,
  .read                  = psc_mfields_vpic_read,
#endif
  .zero_comp             = psc_mfields_vpic_zero_comp,
  .axpy_comp             = psc_mfields_vpic_axpy_comp,
#endif
};


