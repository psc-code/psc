
#include "vpic_config.h"
#include "psc_fields_vpic.h"

#include <psc_method.h>

static int ref_count_fields, ref_count_hydro;

// ----------------------------------------------------------------------
// psc_mfields_vpic_setup

void psc_mfields_vpic_setup(struct psc_mfields *mflds)
{
  struct psc_mfields_vpic *sub = psc_mfields_vpic(mflds);

  psc_mfields_setup_super(mflds);

  mrc_domain_get_patches(mflds->domain, &mflds->nr_patches);
  assert(mflds->nr_patches == 1);
  assert(mflds->ibn[0] == 1);
  assert(mflds->ibn[1] == 1);
  assert(mflds->ibn[2] == 1);
  assert(mflds->first_comp == 0);

  psc_method_get_param_ptr(ppsc->method, "sim", (void **) &sub->sim);

  if (mflds->nr_fields == VPIC_MFIELDS_N_COMP) {
    // make sure we notice if we create a second psc_mfields
    // which would share its memory with the first
    assert(ref_count_fields == 0);
    ref_count_fields++;

    assert(FieldArray::N_COMP == VPIC_MFIELDS_N_COMP);
    sub->vmflds_fields = sub->sim->field_array_;
  } else if (mflds->nr_fields == VPIC_HYDRO_N_COMP) {
    // make sure we notice if we create a second psc_mfields
    // which would share its memory with the first
    assert(ref_count_hydro == 0);
    ref_count_hydro++;

    //assert(HydroArray::N_COMP == VPIC_HYDRO_N_COMP); FIXME
    sub->vmflds_hydro = sub->sim->hydro_array_;
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// psc_mfields_vpic_destroy

void psc_mfields_vpic_destroy(struct psc_mfields *mflds)
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

fields_vpic_t psc_mfields_vpic_get_field_t(struct psc_mfields *mflds, int p)
{
  struct psc_mfields_vpic *sub = mrc_to_subobj(mflds, struct psc_mfields_vpic);
  fields_vpic_t flds;

  // FIXME hacky...
  if (mflds->nr_fields == VPIC_MFIELDS_N_COMP) {
    flds.data = sub->vmflds_fields->getData(flds.ib, flds.im);
    flds.nr_comp = VPIC_MFIELDS_N_COMP;
  } else if (mflds->nr_fields == VPIC_HYDRO_N_COMP) {
    flds.data = sub->vmflds_hydro->getData(flds.ib, flds.im);
    flds.nr_comp = VPIC_HYDRO_N_COMP;
  } else {
    assert(0);
  }
  flds.first_comp = 0;
  assert(mflds->first_comp == 0);

  return flds;
}

