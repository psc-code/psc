
// ======================================================================
// fld1d_t
//
// 1d line of TYPE, based on simple TYPE[] array

typedef struct {
  TYPE *arr;
} FLD1D_(t);

// ----------------------------------------------------------------------
// fld1d_setup

static inline void
FLD1D_(setup)(FLD1D_(t) *f)
{
  assert(!f->arr);

  f->arr = calloc(s_size_1d, sizeof(*f->arr));
  f->arr += s_n_ghosts;
}

// ----------------------------------------------------------------------
// fld1d_is_setup

static inline bool
FLD1D_(is_setup)(FLD1D_(t) f)
{
  return f.arr;
}

// ----------------------------------------------------------------------
// fld1d_free

static inline void
FLD1D_(free)(FLD1D_(t) *f)
{
  f->arr -= s_n_ghosts;
  free(f->arr);
  f->arr = NULL;
}

