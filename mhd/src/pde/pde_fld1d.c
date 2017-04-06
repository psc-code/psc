
// ======================================================================
// defaults for which implementation to use for the various fld1d_*_t

#ifndef OPT_FLD1D
#define OPT_FLD1D OPT_FLD1D_C_ARRAY
#endif

#ifndef OPT_FLD1D_STATE
#define OPT_FLD1D_STATE OPT_FLD1D_PTR_ARRAY
#endif

#ifndef OPT_FLD1D_VEC
#define OPT_FLD1D_VEC OPT_FLD1D_PTR_ARRAY
#endif


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
// fld1d_t
//
// 1d line of mrc_fld_data_t, access with F1()

#if OPT_FLD1D == OPT_FLD1D_MRC_FLD

typedef struct {
  struct mrc_fld *mrc_fld;
} fld1d_t;

#define F1(f, i) MRC_FLD_F1((f).mrc_fld, 0, i)

static inline void
fld1d_setup(fld1d_t *f)
{
  assert(!f->mrc_fld);

  f->mrc_fld = mrc_fld_create_1d(1, s_size_1d);
}

static inline bool
fld1d_is_setup(fld1d_t f)
{
  return f.mrc_fld;
}

static inline void
fld1d_free(fld1d_t *f)
{
  assert(f->mrc_fld);

  mrc_fld_destroy(f->mrc_fld);
  f->mrc_fld = NULL;
}

#elif OPT_FLD1D == OPT_FLD1D_C_ARRAY

#define TYPE mrc_fld_data_t
#define FLD1D_(s) fld1d_ ## s
#include "pde_fld1d_array.c"
#undef TYPE
#undef FLD1D_

#define F1(f, i) ((f).arr[i])

#endif

// ======================================================================
// fld1d_vec_t
//
// 1d line of mrc_fld_data_t, access with F1()

#if OPT_FLD1D_VEC == OPT_FLD1D_MRC_FLD

typedef struct {
  struct mrc_fld *mrc_fld;
} fld1d_vec_t;

#define F1V(f, m, i) MRC_FLD_F1((f).mrc_fld, m, i)

static inline void
fld1d_vec_setup(fld1d_vec_t *f)
{
  assert(!f->mrc_fld);

  f->mrc_fld = mrc_fld_create_1d(3, s_size_1d);
}

static inline bool
fld1d_vec_is_setup(fld1d_vec_t f)
{
  return f.mrc_fld;
}

#elif OPT_FLD1D_VEC == OPT_FLD1D_PTR_ARRAY

typedef struct {
  mrc_fld_data_t **ptr;
} fld1d_vec_t;

#define F1V(f, m, i) ((f).ptr[i][m])

static inline void
fld1d_vec_setup(fld1d_vec_t *f)
{
  mrc_fld_data_t *arr = calloc(s_size_1d * 3, sizeof(*arr));
  f->ptr = calloc(s_size_1d, sizeof(*f->ptr));
  f->ptr += s_n_ghosts;
  for (int i = -s_n_ghosts; i < -s_n_ghosts + s_size_1d; i++) {
    f->ptr[i] = arr;
    arr += 3;
  }
}

static inline bool
fld1d_vec_is_setup(fld1d_vec_t f)
{
  return f.ptr;
}

#endif

// ======================================================================
// fld1d_state_t
//
// 1d line of s_n_comp mrc_fld_data_t, access with F1S

#if OPT_FLD1D_STATE == OPT_FLD1D_MRC_FLD

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

#elif OPT_FLD1D_STATE == OPT_FLD1D_PTR_ARRAY

typedef struct {
  mrc_fld_data_t **ptr;
} fld1d_state_t;

#define F1S(f, m, i) ((f).ptr[i][m])

static inline void
fld1d_state_setup(fld1d_state_t *f)
{
  mrc_fld_data_t *arr = calloc(s_size_1d * s_n_comps, sizeof(*arr));
  f->ptr = calloc(s_size_1d, sizeof(*f->ptr));
  f->ptr += s_n_ghosts;
  for (int i = -s_n_ghosts; i < -s_n_ghosts + s_size_1d; i++) {
    f->ptr[i] = arr;
    arr += s_n_comps;
  }
}

static inline bool
fld1d_state_is_setup(fld1d_state_t f)
{
  return f.ptr;
}

// FIXME, fld1d_state_t and fld1d_vec_t could be consolidated with macro/include hacks

#endif


