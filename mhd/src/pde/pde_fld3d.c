
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
