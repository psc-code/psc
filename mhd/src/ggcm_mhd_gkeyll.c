#include <ggcm_mhd_gkeyll.h>

// parameter setters

void
ggcm_mhd_gkeyll_set_nr_moments(struct ggcm_mhd *mhd,
    int nr_moments)
{
  assert(nr_moments == 5 || nr_moments == 10);
  mhd->par.nr_moments = nr_moments;
}

void
ggcm_mhd_gkeyll_set_nr_fluids(struct ggcm_mhd *mhd,
    int nr_fluids)
{
  assert(nr_fluids > 0 && nr_fluids <= GK_NR_FLUIDS_MAX);
  mhd->par.nr_fluids = nr_fluids;
}

void
ggcm_mhd_gkeyll_set_mass_ratios(struct ggcm_mhd *mhd,
    double mass_ratios_in[])
{
  for (int s = 0; s < ggcm_mhd_gkeyll_nr_fluids(mhd); s++)
    mhd->par.mass_ratios[s] = mass_ratios_in[s];
}

void
ggcm_mhd_gkeyll_set_momentum_ratios(struct ggcm_mhd *mhd,
    double momentum_ratios_in[])
{
  for (int s = 0; s < ggcm_mhd_gkeyll_nr_fluids(mhd); s++)
    mhd->par.momentum_ratios[s] = momentum_ratios_in[s];
}

void
ggcm_mhd_gkeyll_set_pressure_ratios(struct ggcm_mhd *mhd,
    double pressure_ratios_in[])
{
  for (int s = 0; s < ggcm_mhd_gkeyll_nr_fluids(mhd); s++)
    mhd->par.pressure_ratios[s] = pressure_ratios_in[s];
}

// parameter getters

int
ggcm_mhd_gkeyll_nr_moments(struct ggcm_mhd *mhd)
{
  int nr_moments = mhd->par.nr_moments;
  assert(nr_moments == 5 || nr_moments == 10);
  return nr_moments;
}

int
ggcm_mhd_gkeyll_nr_fluids(struct ggcm_mhd *mhd)
{
  int nr_fluids = mhd->par.nr_fluids;
  assert(nr_fluids > 0 && nr_fluids <= GK_NR_FLUIDS_MAX);
  return nr_fluids;
}

double *
ggcm_mhd_gkeyll_mass_ratios(struct ggcm_mhd *mhd)
{
  return mhd->par.mass_ratios;
}

double *
ggcm_mhd_gkeyll_momentum_ratios(struct ggcm_mhd *mhd)
{
  return mhd->par.momentum_ratios;
}

double *
ggcm_mhd_gkeyll_pressure_ratios(struct ggcm_mhd *mhd)
{
  return mhd->par.pressure_ratios;
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

