
#ifndef GGCM_MHD_GKEYLL_H
#define GGCM_MHD_GKEYLL_H

#include <ggcm_mhd_private.h>
#include <ggcm_mhd_defs.h>

// shifts for a species in a five-moment fluid
enum {
  G5M_RRS,
  G5M_RVXS,
  G5M_RVYS,
  G5M_RVZS,
  G5M_UUS,
  G5M_NRS,
};

// shifts within the em fields
enum {
  GK_EX,
  GK_EY,
  GK_EZ,
  GK_BX,
  GK_BY,
  GK_BZ,
  GK_PHI, // correction potentials
  GK_PSI,
  GK_NR_EM,
};


void ggcm_mhd_gkeyll_set_nr_moments(struct ggcm_mhd *mhd, int nr_moments);
void ggcm_mhd_gkeyll_set_nr_fluids(struct ggcm_mhd *mhd, int nr_fluids);

int ggcm_mhd_gkeyll_nr_moments(struct ggcm_mhd *mhd);
int ggcm_mhd_gkeyll_nr_fluids(struct ggcm_mhd *mhd);
float *ggcm_mhd_gkeyll_mass(struct ggcm_mhd *mhd);
float *ggcm_mhd_gkeyll_charge(struct ggcm_mhd *mhd);
float *ggcm_mhd_gkeyll_pressure_ratios(struct ggcm_mhd *mhd);

int ggcm_mhd_gkeyll_fluid_species_index(struct ggcm_mhd *mhd, int species);
void ggcm_mhd_gkeyll_fluid_species_index_all(struct ggcm_mhd *mhd, int indices[]);
int ggcm_mhd_gkeyll_em_fields_index(struct ggcm_mhd *mhd);
void ggcm_mhd_gkeyll_fluid_species_q_m_all(struct ggcm_mhd *mhd, float q_m[]);
void ggcm_mhd_gkeyll_fluid_species_mass_ratios_all(struct ggcm_mhd *mhd, float mass_weight[]);

// FIMXE: better place to declare the convert_primitive funcs
void
ggcm_mhd_convert_primitive_gkeyll_5m_point(struct mrc_fld *fld, int nr_fluids,
    int idx[], float mass_ratios[], float momentum_ratios[],
    float pressure_ratios[], float gamma_m1, int idx_em, int dx, int dy, int dz,
    int ix, int iy, int iz, int p);

#endif
