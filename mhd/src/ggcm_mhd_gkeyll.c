#include <ggcm_mhd_gkeyll.h>

// parameter setters

void
ggcm_mhd_gkeyll_set_nr_moments(struct ggcm_mhd *mhd,
    int nr_moments)
{
  assert(nr_moments == 5 || nr_moments == 10);
  mhd->par.gk_nr_moments = nr_moments;
}

void
ggcm_mhd_gkeyll_set_nr_fluids(struct ggcm_mhd *mhd,
    int nr_fluids)
{
  assert(nr_fluids > 0 && nr_fluids <= GK_NR_FLUIDS_MAX);
  mhd->par.gk_nr_fluids = nr_fluids;
}

// parameter getters

int
ggcm_mhd_gkeyll_nr_moments(struct ggcm_mhd *mhd)
{
  int nr_moments = mhd->par.gk_nr_moments;
  assert(nr_moments == 5 || nr_moments == 10);
  return nr_moments;
}

int
ggcm_mhd_gkeyll_nr_fluids(struct ggcm_mhd *mhd)
{
  int nr_fluids = mhd->par.gk_nr_fluids;
  assert(nr_fluids > 0 && nr_fluids <= GK_NR_FLUIDS_MAX);
  return nr_fluids;
}

float *
ggcm_mhd_gkeyll_mass(struct ggcm_mhd *mhd)
{
  return mhd->par.gk_mass.vals;
}

float *
ggcm_mhd_gkeyll_charge(struct ggcm_mhd *mhd)
{
  return mhd->par.gk_charge.vals;
}

float *
ggcm_mhd_gkeyll_pressure_ratios(struct ggcm_mhd *mhd)
{
  return mhd->par.gk_pressure_ratios.vals;
}

// index calculators

int
ggcm_mhd_gkeyll_fluid_species_index(struct ggcm_mhd *mhd, int species)
{
  // species starts from 0 to nr_fluids-1
  assert(species >=0 && species < ggcm_mhd_gkeyll_nr_fluids(mhd));
  return ggcm_mhd_gkeyll_nr_moments(mhd) * species;
}

void
ggcm_mhd_gkeyll_fluid_species_index_all(struct ggcm_mhd *mhd, int indices[])
{
  for ( int s = 0; s < ggcm_mhd_gkeyll_nr_fluids(mhd); s++) {
    indices[s] = ggcm_mhd_gkeyll_fluid_species_index(mhd, s);
  }
}

int
ggcm_mhd_gkeyll_em_fields_index(struct ggcm_mhd *mhd)
{
  return ggcm_mhd_gkeyll_nr_moments(mhd) * ggcm_mhd_gkeyll_nr_fluids(mhd);
}

void
ggcm_mhd_gkeyll_fluid_species_q_m_all(struct ggcm_mhd *mhd, float q_m[])
{
  for (int sp = 0; sp < mhd->par.gk_nr_fluids; sp++)
    q_m[sp] = mhd->par.gk_charge.vals[sp] / mhd->par.gk_mass.vals[sp];
}

void
ggcm_mhd_gkeyll_fluid_species_mass_ratios_all(struct ggcm_mhd *mhd, float mass_weight[])
{
  float mass_total = 0.;
  for (int sp = 0; sp < mhd->par.gk_nr_fluids; sp++)
    mass_total += mhd->par.gk_mass.vals[sp];

  for (int sp = 0; sp < mhd->par.gk_nr_fluids; sp++)
    mass_weight[sp] = mhd->par.gk_mass.vals[sp] / mass_total;
}

