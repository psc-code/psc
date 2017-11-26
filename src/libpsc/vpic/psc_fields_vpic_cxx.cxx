
#include "psc_fields_vpic.h"

double psc_mfields_vpic_synchronize_tang_e_norm_b(struct psc_mfields *mflds)
{
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  double err;
  TIC err = vmflds->synchronize_tang_e_norm_b(); TOC(synchronize_tang_e_norm_b, 1);
  return err;
}

void psc_mfields_vpic_compute_div_b_err(struct psc_mfields *mflds)
{
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  TIC vmflds->compute_div_b_err(); TOC(compute_div_b_err, 1);
}

double psc_mfields_vpic_compute_rms_div_b_err(struct psc_mfields *mflds)
{
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  double err;
  TIC err = vmflds->compute_rms_div_b_err(); TOC(compute_rms_div_b_err, 1);
  return err;
}

void psc_mfields_vpic_clean_div_b(struct psc_mfields *mflds)
{
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  TIC vmflds->clean_div_b(); TOC(clean_div_b, 1);
}

void psc_mfields_vpic_compute_div_e_err(struct psc_mfields *mflds)
{
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  TIC vmflds->compute_div_e_err(); TOC(compute_div_e_err, 1);
}

double psc_mfields_vpic_compute_rms_div_e_err(struct psc_mfields *mflds)
{
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  double err;
  TIC err = vmflds->compute_rms_div_e_err(); TOC(compute_rms_div_e_err, 1);
  return err;
}

void psc_mfields_vpic_clean_div_e(struct psc_mfields *mflds)
{
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  TIC vmflds->clean_div_e(); TOC(clean_div_e, 1);
}

void psc_mfields_vpic_clear_rhof(struct psc_mfields *mflds)
{
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  TIC vmflds->clear_rhof(); TOC(clear_jf, 1);
}

void psc_mfields_vpic_synchronize_rho(struct psc_mfields *mflds)
{
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  TIC vmflds->synchronize_rho(); TOC(compute_curl_b, 1);
}

void psc_mfields_vpic_compute_rhob(struct psc_mfields *mflds)
{
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  TIC vmflds->compute_rhob(); TOC(compute_rhob, 1);
}

void psc_mfields_vpic_compute_curl_b(struct psc_mfields *mflds)
{
  FieldArray *vmflds = psc_mfields_vpic(mflds)->vmflds_fields;

  TIC vmflds->compute_curl_b(); TOC(compute_curl_b, 1);
}

