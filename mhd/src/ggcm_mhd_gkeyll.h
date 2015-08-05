
enum {
  G5M_RRE,
  G5M_RVXE,
  G5M_RVYE,
  G5M_RVZE,
  G5M_UUE,
  G5M_RRI,
  G5M_RVXI,
  G5M_RVYI,
  G5M_RVZI,
  G5M_UUI,
  G5M_EX,
  G5M_EY,
  G5M_EZ,
  G5M_BX,
  G5M_BY,
  G5M_BZ,
  G5M_PHI,
  G5M_PSI,
  G5M_NR,
};

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

int
mrc_fld_gkeyll_nr_moments(struct mrc_fld *f);

int
mrc_fld_gkeyll_nr_fluids(struct mrc_fld *f);

void
mrc_fld_gkeyll_mass_ratios(struct mrc_fld *f, float mass_ratios[]);

void
mrc_fld_gkeyll_momentum_ratios(struct mrc_fld *f, float momentum_ratios[]);

void
mrc_fld_gkeyll_temperature_ratios(struct mrc_fld *f, float temperature_ratios[]);

void
mrc_fld_gkeyll_set_nr_moments(struct mrc_fld *f, int nr_moments);

void
mrc_fld_gkeyll_get_nr_fluids(struct mrc_fld *f, int nr_fluids);

int
mrc_fld_gkeyll_electron_index_two_fluids(struct mrc_fld *f, int m_beg);

int
mrc_fld_gkeyll_ion_index_two_fluids(struct mrc_fld *f, int m_beg);

int
mrc_fld_gkeyll_species_index(struct mrc_fld *f, int m_beg, int species);

void
mrc_fld_gkeyll_species_index_all(struct mrc_fld *f, int m_beg, int indices[]);

int
mrc_fld_gkeyll_em_index(struct mrc_fld *f, int m_beg);
