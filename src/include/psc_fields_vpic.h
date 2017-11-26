
#ifndef PSC_FIELDS_VPIC_H
#define PSC_FIELDS_VPIC_H

#include "psc_fields_private.h"

#include <mrc_common.h>

#define FTYPE FTYPE_VPIC
#include "psc_fields_common.h"
#undef FTYPE

#include "vpic_iface.h"

BEGIN_C_DECLS

// ----------------------------------------------------------------------

struct psc_mfields_vpic {
  Simulation *sim;
  FieldArray *vmflds_fields;
  struct HydroArray *vmflds_hydro;
};

#define psc_mfields_vpic(mflds) ({					\
      assert((struct psc_mfields_ops *) mflds->obj.ops == &psc_mfields_vpic_ops); \
      mrc_to_subobj(mflds, struct psc_mfields_vpic);			\
    })


// These should otherwise be static, but they are defined in a separate C++ source file

double psc_mfields_vpic_synchronize_tang_e_norm_b(struct psc_mfields *mflds);

void psc_mfields_vpic_compute_div_b_err(struct psc_mfields *mflds);
double psc_mfields_vpic_compute_rms_div_b_err(struct psc_mfields *mflds);
void psc_mfields_vpic_clean_div_b(struct psc_mfields *mflds);

void psc_mfields_vpic_compute_div_e_err(struct psc_mfields *mflds);
double psc_mfields_vpic_compute_rms_div_e_err(struct psc_mfields *mflds);
void psc_mfields_vpic_clean_div_e(struct psc_mfields *mflds);

void psc_mfields_vpic_clear_rhof(struct psc_mfields *mflds);
void psc_mfields_vpic_synchronize_rho(struct psc_mfields *mflds);

void psc_mfields_vpic_compute_rhob(struct psc_mfields *mflds);
void psc_mfields_vpic_compute_curl_b(struct psc_mfields *mflds);

END_C_DECLS

#endif
