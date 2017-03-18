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

// ----------------------------------------------------------------------	
// convert_primitive_5m_point_comove

void
convert_primitive_5m_point_comove(float state[], float vals[], int nr_fluids, int nr_moments,
    float mass[], float charge[], float pressure_ratios[], float gamm)
{
  float mass_ratios[nr_fluids];
  float mass_total = 0.;
  int idx[nr_fluids];

  for (int s = 0; s < nr_fluids; s++) {
    mass_total += mass[s];
    idx[s] = s * nr_moments;
  }
  for (int s = 0; s < nr_fluids; s++)
    mass_ratios[s] = mass[s] / mass_total;

  int idx_em = nr_fluids * nr_moments;

  float rr = vals[RR];
  float vx = vals[VX];
  float vy = vals[VY];
  float vz = vals[VZ];
  float pp = vals[PP];
  float bx = vals[BX];
  float by = vals[BY];
  float bz = vals[BZ];

  state[idx_em + GK_EX] = - vy * bz + vz * by;
  state[idx_em + GK_EY] = - vz * bx + vx * bz;
  state[idx_em + GK_EZ] = - vx * by + vy * bx;

  state[idx_em + GK_BX] = bx;
  state[idx_em + GK_BY] = by;
  state[idx_em + GK_BZ] = bz;

  for (int s = 0; s < nr_fluids; s++) {
    state[idx[s] + G5M_RRS ] = rr * mass_ratios[s];
    state[idx[s] + G5M_RVXS] = rr * mass_ratios[s] * vx;
    state[idx[s] + G5M_RVYS] = rr * mass_ratios[s] * vy;
    state[idx[s] + G5M_RVZS] = rr * mass_ratios[s] * vz;
    state[idx[s] + G5M_UUS ] = pp * pressure_ratios[s] / (gamm - 1.)
      + .5 * (sqr(state[idx[s] + G5M_RVXS])
            + sqr(state[idx[s] + G5M_RVYS])
            + sqr(state[idx[s] + G5M_RVZS])) / state[idx[s] + G5M_RRS];
  }

  state[idx_em + GK_PHI] = 0.;
  state[idx_em + GK_PSI] = 0.;
} 

