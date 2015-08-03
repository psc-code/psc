
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

int
mrc_fld_gkeyll_nr_moments(struct mrc_fld *f);

int
mrc_fld_gkeyll_nr_fluids(struct mrc_fld *f);

int
mrc_fld_gkeyll_electron_index_two_fluids(struct mrc_fld *f, int m_beg);

int
mrc_fld_gkeyll_ion_index_two_fluids(struct mrc_fld *f, int m_beg);

int
mrc_fld_gkeyll_species_index(struct mrc_fld *f, int m_beg, int species);

int
mrc_fld_gkeyll_em_index(struct mrc_fld *f, int m_beg);
