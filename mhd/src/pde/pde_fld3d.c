
// ======================================================================
// fld3d_t
//
// 3d patch of mrc_fld_data_t, access with F3S

typedef struct {
  mrc_fld_data_t *arr_off;
  struct mrc_fld *mrc_fld;
} fld3d_t;

//#define F3S(f, m, i,j,k) M3((f).mrc_fld, m, i,j,k, (f).p)

#define F3S(f, m, i,j,k)						\
  (((f).arr_off)[((((m)) * s_lgdims[2] +				\
		   (k)) * s_lgdims[1] +					\
		  (j)) * s_lgdims[0] +					\
		 (i)])

static inline void
fld3d_setup(fld3d_t *f, struct mrc_fld *fld)
{
  f->mrc_fld = fld;
}

static inline void
fld3d_setup_tmp(fld3d_t *f, int n_comps)
{
  f->mrc_fld = NULL;
  f->arr_off = calloc(s_lgdims[0] * s_lgdims[1] * s_lgdims[2] * n_comps,
		      sizeof(*f->arr_off));
  f->arr_off += (((0
		   * s_lgdims[2] + s_sw[2] )
		  * s_lgdims[1] + s_sw[1])
		 * s_lgdims[0] + s_sw[0]);
}

static inline void
fld3d_setup_view(fld3d_t *f, fld3d_t base, int m)
{
  assert(base.arr_off);
  f->arr_off = base.arr_off + m * s_lgdims[0] * s_lgdims[1] * s_lgdims[2];
  f->mrc_fld = NULL;
}

static inline fld3d_t
fld3d_make_view(fld3d_t base, int m)
{
  fld3d_t f;
  assert(base.arr_off);
  f.arr_off = base.arr_off + m * s_lgdims[0] * s_lgdims[1] * s_lgdims[2];
  f.mrc_fld = NULL;
  return f;
}

static inline void
fld3d_get(fld3d_t *f, int p)
{
  assert(f->mrc_fld);
  f->arr_off = f->mrc_fld->nd.arr_off;
  assert(f->mrc_fld->_stride[0] == 1);
  assert(f->mrc_fld->_stride[1] == s_lgdims[0]);
  assert(f->mrc_fld->_stride[2] == s_lgdims[0] * s_lgdims[1]);
  assert(f->mrc_fld->_stride[3] == s_lgdims[0] * s_lgdims[1] * s_lgdims[2]);
  f->arr_off += p * s_lgdims[0] * s_lgdims[1] * s_lgdims[2];
}

// FIXME, do we require fld3d_put to close fld3d_get?
// If so, we should do it consistently -- if not we should just get rid of it

static inline void
fld3d_put(fld3d_t *f, int p)
{
  assert(f->mrc_fld);
  f->arr_off = NULL;
}

static inline bool
fld3d_is_setup(fld3d_t f)
{
  return f.mrc_fld;
}

static void _mrc_unused
fld3d_get_list(int p, fld3d_t *list[])
{
  for (fld3d_t **f = list; *f; f++) {
    fld3d_get(*f, p);
  }
}

static void _mrc_unused
fld3d_put_list(int p, fld3d_t *list[])
{
  for (fld3d_t **f = list; *f; f++) {
    fld3d_put(*f, p);
  }
}

#if OPT_TMP == OPT_TMP_COMPAT

#define fld3d_setup_tmp_compat(p_tmp, n_comps, m) do {	\
    if (!(p_tmp)->arr_off) {				\
      fld3d_setup_view(p_tmp, s_p_f, m);		\
    }							\
  } while (0)

#else

#define fld3d_setup_tmp_compat(p_tmp, n_comps, m) do {	\
    if (!(p_tmp)->arr_off) {				\
      fld3d_setup_tmp(p_tmp, n_comps);			\
    }							\
  } while (0)

#endif

