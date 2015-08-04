/* 
 * utility functions relevant to an mrc_fld with gkeyll data
 */

#include <mrc_fld.h>
#include <ggcm_mhd_gkeyll.h>

// ----------------------------------------------------------------------
// mrc_fld_gkeyll_nr_moments

int
mrc_fld_gkeyll_nr_moments(struct mrc_fld *f)
{
  int nr_moments = 0;
  mrc_fld_get_param_int(f, "nr_moments", &nr_moments);
  return nr_moments;
}

// ----------------------------------------------------------------------
// mrc_fld_gkeyll_nr_moments

int
mrc_fld_gkeyll_nr_fluids(struct mrc_fld *f)
{
  int nr_fluids = 0;
  mrc_fld_get_param_int(f, "nr_fluids", &nr_fluids);
  return nr_fluids;
}

// ----------------------------------------------------------------------
// mrc_fld_gkeyll_set_nr_moments

void
mrc_fld_gkeyll_set_nr_moments(struct mrc_fld *f, int nr_moments)
{
  // TODO: overwrite existing value
  mrc_fld_dict_add_int(f, "nr_moments", nr_moments);
}

// ----------------------------------------------------------------------
// mrc_fld_gkeyll_set_nr_moments

void
mrc_fld_gkeyll_set_nr_fluids(struct mrc_fld *f, int nr_fluids)
{
  mrc_fld_dict_add_int(f, "nr_fluids", nr_fluids);
}

// ----------------------------------------------------------------------
// mrc_fld_gkeyll_electron_index_two_fluids
//
// first index of the electron moments for a two-fluid case

int
mrc_fld_gkeyll_electron_index_two_fluids(struct mrc_fld *f, int m_beg)
{
  assert(mrc_fld_gkeyll_nr_fluids(f) == 2);
  return m_beg;
}

// ----------------------------------------------------------------------
// mrc_fld_gkeyll_ion_index_two_fluids
//
// first index of the electron moments for a two-fluid case

int
mrc_fld_gkeyll_ion_index_two_fluids(struct mrc_fld *f, int m_beg)
{
  assert(mrc_fld_gkeyll_nr_fluids(f) == 2);
  return m_beg + mrc_fld_gkeyll_nr_moments(f);
}

// ----------------------------------------------------------------------
// mrc_fld_gkeyll_ion_index_two_fluids
//
// first index of the electron moments for a two-fluid case

int
mrc_fld_gkeyll_species_index(struct mrc_fld *f, int m_beg, int species)
{
  assert(species < mrc_fld_gkeyll_nr_fluids(f));
  return m_beg + mrc_fld_gkeyll_nr_moments(f) * (species - 1);
}

// ----------------------------------------------------------------------
// mrc_fld_gkeyll_em_index
//
// first index of the EM fields

int
mrc_fld_gkeyll_em_index(struct mrc_fld *f, int m_beg)
{
  return m_beg + mrc_fld_gkeyll_nr_moments(f) * mrc_fld_gkeyll_nr_fluids(f);
}

