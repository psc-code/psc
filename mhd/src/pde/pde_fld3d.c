
// ======================================================================
// fld3d_t
//
// 3d patch of mrc_fld_data_t, access with F3

typedef struct {
  struct mrc_fld *mrc_fld;
  int p;
} fld3d_t;

#define F3S(f, m, i,j,k) M3((f).mrc_fld, m, i,j,k, (f).p)

static inline void
fld3d_setup(fld3d_t *f)
{
  f->mrc_fld = NULL;
}

static inline void
fld3d_setup_tmp(fld3d_t *f, int n_comps)
{
  f->mrc_fld = mrc_fld_create(MPI_COMM_SELF);
  mrc_fld_set_type(f->mrc_fld, FLD_TYPE);
  mrc_fld_set_param_int_array(f->mrc_fld, "dims", 5,
			      (int[]) { s_ldims[0], s_ldims[1], s_ldims[2], n_comps, 1 });
  mrc_fld_set_param_int_array(f->mrc_fld, "sw", 5,
			      (int[]) { s_sw[0], s_sw[1], s_sw[2], 0, 0 });
  mrc_fld_setup(f->mrc_fld);
  f->p = 0;
}

static inline void
fld3d_get(fld3d_t *f, struct mrc_fld *fld, int p)
{
  f->mrc_fld = fld;
  f->p = p;
}

static inline void
fld3d_put(fld3d_t *f, struct mrc_fld *fld, int p)
{
  f->mrc_fld = NULL;
}

static inline bool
fld3d_is_setup(fld3d_t f)
{
  return f.mrc_fld;
}
