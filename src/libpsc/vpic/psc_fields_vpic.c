
#include "psc_fields_vpic.h"

#include "psc_fields_c.h"
#include "psc_fields_single.h"
#include "psc.h"
#include "psc_method.h"
#include "psc_particles_vpic.h"

#include "mrc_domain.h"
#include "mrc_bits.h"

#include "vpic_iface.h"

static const int map_psc2vpic[VPIC_MFIELDS_N_COMP] = {
  [JXI] = VPIC_MFIELDS_JFX, [JYI] = VPIC_MFIELDS_JFY, [JZI] = VPIC_MFIELDS_JFZ,
  [EX]  = VPIC_MFIELDS_EX , [EY]  = VPIC_MFIELDS_EY , [EZ]  = VPIC_MFIELDS_EZ,
  [HX]  = VPIC_MFIELDS_BX , [HY]  = VPIC_MFIELDS_BY , [HZ]  = VPIC_MFIELDS_BZ,

  [9]   = VPIC_MFIELDS_TCAX, [10]  = VPIC_MFIELDS_TCAY, [11]  = VPIC_MFIELDS_TCAZ,
  [12]  = VPIC_MFIELDS_DIV_E_ERR, [13] = VPIC_MFIELDS_DIV_B_ERR,
  [14]  = VPIC_MFIELDS_RHOB     , [15] = VPIC_MFIELDS_RHOF,
  [16]  = 16, [17] = 17, [18] = 18, [19] = 19,
};

static int ref_count_fields, ref_count_hydro;

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
  assert(mflds->ibn[0] == 1);
  assert(mflds->ibn[1] == 1);
  assert(mflds->ibn[2] == 1);
  assert(mflds->first_comp == 0);

  struct Simulation *sim;
  psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sim);

  if (mflds->nr_fields == VPIC_MFIELDS_N_COMP) {
    // make sure we notice if we create a second psc_mfields
    // which would share its memory with the first
    assert(ref_count_fields == 0);
    ref_count_fields++;

    sub->vmflds_fields = vpic_mfields_new_fields_array(sim);
  } else if (mflds->nr_fields == VPIC_HYDRO_N_COMP) {
    // make sure we notice if we create a second psc_mfields
    // which would share its memory with the first
    assert(ref_count_hydro == 0);
    ref_count_hydro++;

    sub->vmflds_hydro = vpic_mfields_new_hydro_array(sim);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_vpic_destroy

static void
psc_mfields_vpic_destroy(struct psc_mfields *mflds)
{
  if (mflds->nr_fields == VPIC_MFIELDS_N_COMP) {
    ref_count_fields--;
  } else if (mflds->nr_fields == VPIC_HYDRO_N_COMP) {
    ref_count_hydro--;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_vpic_get_field_t

fields_vpic_t
psc_mfields_vpic_get_field_t(struct psc_mfields *mflds, int p)
{
  assert((struct psc_mfields_ops *) mflds->obj.ops == &psc_mfields_vpic_ops);
  struct psc_mfields_vpic *sub = mrc_to_subobj(mflds, struct psc_mfields_vpic);
  fields_vpic_t flds;

  // FIXME hacky...
  if (mflds->nr_fields == VPIC_MFIELDS_N_COMP) {
    flds.data = vpic_mfields_get_data(sub->vmflds_fields, flds.ib, flds.im);
    flds.nr_comp = VPIC_MFIELDS_N_COMP;
  } else if (mflds->nr_fields == VPIC_HYDRO_N_COMP) {
    flds.data = vpic_mfields_hydro_get_data(sub->vmflds_hydro, flds.ib, flds.im);
    flds.nr_comp = VPIC_HYDRO_N_COMP;
  } else {
    assert(0);
  }
  flds.first_comp = 0;
  assert(mflds->first_comp == 0);

  return flds;
}

// ----------------------------------------------------------------------
// forward to vmflds

static double
psc_mfields_vpic_synchronize_tang_e_norm_b(struct psc_mfields *mflds)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  return vpic_mfields_synchronize_tang_e_norm_b(vmflds);
}

static void
psc_mfields_vpic_compute_div_b_err(struct psc_mfields *mflds)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  vpic_mfields_compute_div_b_err(vmflds);
}

static double
psc_mfields_vpic_compute_rms_div_b_err(struct psc_mfields *mflds)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  return vpic_mfields_compute_rms_div_b_err(vmflds);
}

static void
psc_mfields_vpic_clean_div_b(struct psc_mfields *mflds)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  vpic_mfields_clean_div_b(vmflds);
}

static void
psc_mfields_vpic_compute_div_e_err(struct psc_mfields *mflds)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  vpic_mfields_compute_div_e_err(vmflds);
}

static double
psc_mfields_vpic_compute_rms_div_e_err(struct psc_mfields *mflds)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  return vpic_mfields_compute_rms_div_e_err(vmflds);
}

static void
psc_mfields_vpic_clean_div_e(struct psc_mfields *mflds)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  vpic_mfields_clean_div_e(vmflds);
}

static void
psc_mfields_vpic_clear_rhof(struct psc_mfields *mflds)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  vpic_mfields_clear_rhof(vmflds);
}

static void
psc_mfields_vpic_accumulate_rho_p(struct psc_mfields *mflds,
				  struct psc_mparticles *mprts)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;
  struct vpic_mparticles *vmprts = psc_mparticles_vpic(mprts)->vmprts;

  vpic_mfields_accumulate_rho_p(vmflds, vmprts);
}

static void
psc_mfields_vpic_synchronize_rho(struct psc_mfields *mflds)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  vpic_mfields_synchronize_rho(vmflds);
}

static void
psc_mfields_vpic_compute_rhob(struct psc_mfields *mflds)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  vpic_mfields_compute_rhob(vmflds);
}

static void
psc_mfields_vpic_compute_curl_b(struct psc_mfields *mflds)
{
  struct vpic_mfields *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  vpic_mfields_compute_curl_b(vmflds);
}

// ======================================================================
// convert from/to "single"

static void
psc_mfields_vpic_copy_from_single(struct psc_mfields *mflds, struct psc_mfields *mflds_single,
				  int mb, int me)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_single_t flds_single = fields_single_t_mflds(mflds_single, p);
    fields_vpic_t flds = fields_vpic_t_mflds(mflds, p);

    // FIXME, hacky way to distinguish whether we want
    // to copy the field into the standard PSC component numbering
    if (mflds->nr_fields == VPIC_MFIELDS_N_COMP) {
      for (int m = mb; m < me; m++) {
	int m_vpic = map_psc2vpic[m];
	for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
	  for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
	    for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	      _F3_VPIC(flds, m_vpic, jx,jy,jz) = _F3_S(flds_single, m, jx,jy,jz);
	    }
	  }
	}
      }
    } else {
      assert(0);
    }
  }
}

static void
psc_mfields_vpic_copy_to_single(struct psc_mfields *mflds, struct psc_mfields *mflds_single,
			   int mb, int me)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_single_t flds_single = fields_single_t_mflds(mflds_single, p);
    fields_vpic_t flds = fields_vpic_t_mflds(mflds, p);

    int ib[3], ie[3];
    for (int d = 0; d < 3; d++) {
      ib[d] = MAX(flds.ib[d], flds_single.ib[d]);
      ie[d] = MIN(flds.ib[d] + flds.im[d], flds_single.ib[d] + flds_single.im[d]);
    }

    // FIXME, hacky way to distinguish whether we want
    // to copy the field into the standard PSC component numbering or,
    // as in this case, just copy one-to-one
    if (mflds->nr_fields == VPIC_HYDRO_N_COMP) {
      for (int m = mb; m < me; m++) {
	for (int jz = ib[2]; jz < ie[2]; jz++) {
	  for (int jy = ib[1]; jy < ie[1]; jy++) {
	    for (int jx = ib[0]; jx < ie[0]; jx++) {
	      _F3_S(flds_single, m, jx,jy,jz) = _F3_VPIC(flds, m, jx,jy,jz);
	    }
	  }
	}
      }
    } else if (mflds->nr_fields == VPIC_MFIELDS_N_COMP) {
      for (int m = mb; m < me; m++) {
	int m_vpic = map_psc2vpic[m];
	for (int jz = ib[2]; jz < ie[2]; jz++) {
	  for (int jy = ib[1]; jy < ie[1]; jy++) {
	    for (int jx = ib[0]; jx < ie[0]; jx++) {
	      _F3_S(flds_single, m, jx,jy,jz) = _F3_VPIC(flds, m_vpic, jx,jy,jz);
	    }
	  }
	}
      }
    } else {
      assert(0);
    }
  }
}

// ======================================================================
// convert from/to "c"

static void
psc_mfields_vpic_copy_from_c(struct psc_mfields *mflds, struct psc_mfields *mflds_c,
			    int mb, int me)
{
  if (me > mb) {
    assert(0);
  }
  for (int p = 0; p < mflds->nr_patches; p++) {
#if 0
    fields_c_t flds_c = fields_c_t_mflds(mflds_c, p);
    fields_vpic_t flds = fields_vpic_t_mflds(mflds, p);

    for (int m = mb; m < me; m++) {
      for (int jz = flds.ib[2]; jz < flds.ib[2] + flds.im[2]; jz++) {
	for (int jy = flds.ib[1]; jy < flds.ib[1] + flds.im[1]; jy++) {
	  for (int jx = flds.ib[0]; jx < flds.ib[0] + flds.im[0]; jx++) {
	    _F3_VPIC(flds, m, jx,jy,jz) = _F3_C(flds_c, m, jx,jy,jz);
	  }
	}
      }
    }
#endif
  }
}

static void
psc_mfields_vpic_copy_to_c(struct psc_mfields *mflds, struct psc_mfields *mflds_c,
			   int mb, int me)
{
  for (int p = 0; p < mflds->nr_patches; p++) {
    fields_c_t flds_c = fields_c_t_mflds(mflds_c, p);
    fields_vpic_t flds = fields_vpic_t_mflds(mflds, p);

    int ib[3], ie[3];
    for (int d = 0; d < 3; d++) {
      ib[d] = MAX(flds.ib[d], flds_c.ib[d]);
      ie[d] = MIN(flds.ib[d] + flds.im[d], flds_c.ib[d] + flds_c.im[d]);
    }

    // FIXME, hacky way to distinguish whether we want
    // to copy the field into the standard PSC component numbering or,
    // as in this case, just copy one-to-one
    if (mflds->nr_fields == VPIC_HYDRO_N_COMP) {
      for (int m = mb; m < me; m++) {
	for (int jz = ib[2]; jz < ie[2]; jz++) {
	  for (int jy = ib[1]; jy < ie[1]; jy++) {
	    for (int jx = ib[0]; jx < ie[0]; jx++) {
	      _F3_C(flds_c, m, jx,jy,jz) = _F3_VPIC(flds, m, jx,jy,jz);
	    }
	  }
	}
      }
    } else if (mflds->nr_fields == VPIC_MFIELDS_N_COMP) {
      for (int m = mb; m < me; m++) {
	int m_vpic = map_psc2vpic[m];
	for (int jz = ib[2]; jz < ie[2]; jz++) {
	  for (int jy = ib[1]; jy < ie[1]; jy++) {
	    for (int jx = ib[0]; jx < ie[0]; jx++) {
	      _F3_C(flds_c, m, jx,jy,jz) = _F3_VPIC(flds, m_vpic, jx,jy,jz);
	    }
	  }
	}
      }
    } else {
      assert(0);
    }
  }
}

// ======================================================================
// psc_mfields: subclass "vpic"

static struct mrc_obj_method psc_mfields_vpic_methods[] = {
  MRC_OBJ_METHOD("copy_to_c"       , psc_mfields_vpic_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c"     , psc_mfields_vpic_copy_from_c),
  MRC_OBJ_METHOD("copy_to_single"  , psc_mfields_vpic_copy_to_single),
  MRC_OBJ_METHOD("copy_from_single", psc_mfields_vpic_copy_from_single),
  {}
};

struct psc_mfields_ops psc_mfields_vpic_ops = {
  .name                  = "vpic",
  .size                  = sizeof(struct psc_mfields_vpic),
  .methods               = psc_mfields_vpic_methods,
  .setup                 = psc_mfields_vpic_setup,
  .destroy               = psc_mfields_vpic_destroy,
#if 0
#ifdef HAVE_LIBHDF5_HL
  .write                 = psc_mfields_vpic_write,
  .read                  = psc_mfields_vpic_read,
#endif
  .zero_comp             = psc_mfields_vpic_zero_comp,
  .axpy_comp             = psc_mfields_vpic_axpy_comp,
#endif
  .synchronize_tang_e_norm_b = psc_mfields_vpic_synchronize_tang_e_norm_b,
  .compute_div_b_err         = psc_mfields_vpic_compute_div_b_err,
  .compute_rms_div_b_err     = psc_mfields_vpic_compute_rms_div_b_err,
  .clean_div_b               = psc_mfields_vpic_clean_div_b,
  .compute_div_e_err         = psc_mfields_vpic_compute_div_e_err,
  .compute_rms_div_e_err     = psc_mfields_vpic_compute_rms_div_e_err,
  .clean_div_e               = psc_mfields_vpic_clean_div_e,
  .clear_rhof                = psc_mfields_vpic_clear_rhof,
  .accumulate_rho_p          = psc_mfields_vpic_accumulate_rho_p,
  .synchronize_rho           = psc_mfields_vpic_synchronize_rho,
  .compute_rhob              = psc_mfields_vpic_compute_rhob,
  .compute_curl_b            = psc_mfields_vpic_compute_curl_b,
};


