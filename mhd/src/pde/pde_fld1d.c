
// ======================================================================
// the following provides a (optional) basis for mrc_fld based 
//
// line fields (fld1d_*_t)

#define MRC_FLD_F1(f, m, i) MRC_FLD(f, mrc_fld_data_t, m,i,0,0,0)

// ----------------------------------------------------------------------
// mrc_fld_create_1d

static inline struct mrc_fld *
mrc_fld_create_1d(int nr_comps, int size)
{
  struct mrc_fld *f = mrc_fld_create(MPI_COMM_NULL);
  mrc_fld_set_type(f, FLD_TYPE);
  mrc_fld_set_param_int_array(f , "dims", 2, (int []) { nr_comps, size });
  mrc_fld_set_param_int_array(f , "offs", 2, (int []) { 0, -s_n_ghosts });
  mrc_fld_setup(f);

  return f;
}

// ======================================================================
// fld1d_state_t
//
// 1d line of s_n_comp doubles, access with F1S

typedef struct {
  struct mrc_fld *mrc_fld;
} fld1d_state_t;

#define F1S(f, m, i) MRC_FLD_F1((f).mrc_fld, m, i)

static inline void
fld1d_state_setup(fld1d_state_t *f)
{
  f->mrc_fld = mrc_fld_create_1d(s_n_comps, s_size_1d);
}

static inline bool
fld1d_state_is_setup(fld1d_state_t f)
{
  return f.mrc_fld;
}

