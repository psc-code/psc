
#include <mrc_fld.h>
#include <mrc_io.h>
#include <mrc_params.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// ======================================================================
// mrc_fld

// ----------------------------------------------------------------------
// mrc_fld_destroy

static void
_mrc_fld_destroy(struct mrc_fld *fld)
{
  if (!fld->_with_array) {
    free(fld->_arr);
  }
  fld->_arr = NULL;
}

static void
_mrc_f3_destroy(struct mrc_f3 *f3)
{
  if (!f3->_with_array) {
    free(f3->_arr);
  }
  f3->_arr = NULL;

  for (int m = 0; m < f3->_nr_allocated_comp_name; m++) {
    free(f3->_comp_name[m]);
  }
  free(f3->_comp_name);
}

// ----------------------------------------------------------------------
// mrc_fld_setup

static void
_mrc_fld_setup(struct mrc_fld *fld)
{
  if (fld->_offs.nr_vals == 0) {
    fld->_offs.nr_vals = fld->_dims.nr_vals;
  }
  if (fld->_sw.nr_vals == 0) {
    fld->_sw.nr_vals = fld->_dims.nr_vals;
  }
  assert(fld->_dims.nr_vals == fld->_offs.nr_vals &&
	 fld->_dims.nr_vals == fld->_sw.nr_vals);
  assert(fld->_dims.nr_vals <= MRC_FLD_MAXDIMS);

  fld->_len = 1;
  for (int d = 0; d < fld->_dims.nr_vals; d++) {
    fld->_ghost_offs[d] = fld->_offs.vals[d] - fld->_sw.vals[d];
    fld->_ghost_dims[d] = fld->_dims.vals[d] + 2 * fld->_sw.vals[d];
    fld->_len *= fld->_ghost_dims[d];
  }

  if (!fld->_with_array) {
    fld->_arr = calloc(fld->_len, fld->_size_of_type);
  }
}

static void
_mrc_f3_setup(struct mrc_f3 *f3)
{
  if (f3->_offs.nr_vals == 0) {
    f3->_offs.nr_vals = f3->_dims.nr_vals;
  }
  if (f3->_sw.nr_vals == 0) {
    f3->_sw.nr_vals = f3->_dims.nr_vals;
  }
  assert(f3->_dims.nr_vals == f3->_offs.nr_vals &&
	 f3->_dims.nr_vals == f3->_sw.nr_vals);
  assert(f3->_dims.nr_vals == 4);

  for (int d = 0; d < 3; d++) {
    f3->_ghost_offs[d] = f3->_offs.vals[d] - f3->_sw.vals[d];
    f3->_ghost_dims[d] = f3->_dims.vals[d] + 2 * f3->_sw.vals[d];
  }
  f3->_len = f3->_ghost_dims[0] * f3->_ghost_dims[1] * f3->_ghost_dims[2] * f3->_nr_comp;

  if (!f3->_arr) {
    f3->_arr = calloc(f3->_len, sizeof(float));
    f3->_with_array = false;
  } else {
    f3->_with_array = true;
  }
}

// ----------------------------------------------------------------------
// mrc_fld_set_array

void
mrc_fld_set_array(struct mrc_fld *fld, void *arr)
{
  assert(!fld->_arr);
  fld->_arr = arr;
  fld->_with_array = true;
}

void
mrc_f3_set_array(struct mrc_f3 *f3, float *arr)
{
  assert(!f3->_arr);
  f3->_arr = arr;
}

// ----------------------------------------------------------------------
// mrc_fld_write

static void
_mrc_fld_write(struct mrc_fld *fld, struct mrc_io *io)
{
  mrc_io_write_fld(io, mrc_io_obj_path(io, fld), fld);
}

static void
_mrc_f3_write(struct mrc_f3 *f3, struct mrc_io *io)
{
  mrc_io_write_ref(io, f3, "domain", f3->_domain);
  mrc_io_write_f3(io, mrc_io_obj_path(io, f3), f3, 1.);
}

// ----------------------------------------------------------------------
// mrc_fld_read

static void
_mrc_fld_read(struct mrc_fld *fld, struct mrc_io *io)
{
  // if we're reading back stuff, there's no way that _with_array
  // would work, so we'll allocate our own array instead.
  fld->_with_array = false;
  mrc_fld_setup(fld);
  mrc_io_read_fld(io, mrc_io_obj_path(io, fld), fld);
}

static void
_mrc_f3_read(struct mrc_f3 *f3, struct mrc_io *io)
{
  f3->_domain = mrc_io_read_ref(io, f3, "domain", mrc_domain);

  // rely on domain rather than read params
  // since the domain may be different (# of procs) when
  // we're reading things back
  // basically, we should use mrc_domain_f3_create()
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(f3->_domain, &nr_patches);
  assert(nr_patches == 1);
  int *dims = patches[0].ldims;
  int nr_comps = mrc_f3_nr_comps(f3);
  mrc_f3_set_param_int_array(f3, "dims", 4,
			     (int[4]) { dims[0], dims[1], dims[2], nr_comps });
  mrc_f3_set_param_int_array(f3, "offs", 4,
			     (int[4]) { 0, 0, 0, 0 });
  mrc_f3_setup(f3);
  // FIXME, the whole _comp_name business is screwy here
  f3->_comp_name = calloc(f3->_nr_comp, sizeof(*f3->_comp_name));
  f3->_nr_allocated_comp_name = f3->_nr_comp;
  mrc_io_read_f3(io, mrc_io_obj_path(io, f3), f3);
}

// ======================================================================
// mrc_fld subclasses

#define MAKE_MRC_FLD_TYPE(type, TYPE)			\
							\
  static void						\
  mrc_fld_##type##_create(struct mrc_fld *fld)		\
  {							\
    fld->_data_type = MRC_NT_##TYPE;			\
    fld->_size_of_type = sizeof(type);			\
  }							\
  							\
  static struct mrc_fld_ops mrc_fld_##type##_ops = {	\
    .name                  = #type,			\
    .create                = mrc_fld_##type##_create,	\
  };							\

MAKE_MRC_FLD_TYPE(float, FLOAT)
MAKE_MRC_FLD_TYPE(double, DOUBLE)
MAKE_MRC_FLD_TYPE(int, INT)

// ----------------------------------------------------------------------
// mrc_fld_init

static void
mrc_fld_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_fld, &mrc_fld_float_ops);
  mrc_class_register_subclass(&mrc_class_mrc_fld, &mrc_fld_double_ops);
  mrc_class_register_subclass(&mrc_class_mrc_fld, &mrc_fld_int_ops);
}

// ----------------------------------------------------------------------
// mrc_class_mrc_fld

#define VAR(x) (void *)offsetof(struct mrc_fld, x)
static struct param mrc_fld_descr[] = {
  { "offs"            , VAR(_offs)        , PARAM_INT_ARRAY(0, 0) },
  { "dims"            , VAR(_dims)        , PARAM_INT_ARRAY(0, 0) },
  { "sw"              , VAR(_sw)          , PARAM_INT_ARRAY(0, 0) },

  { "size_of_type"    , VAR(_size_of_type), MRC_VAR_INT           },
  { "len"             , VAR(_len)         , MRC_VAR_INT           },
  { "with_array"      , VAR(_with_array)  , MRC_VAR_BOOL          },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// mrc_fld class description

struct mrc_class_mrc_fld mrc_class_mrc_fld = {
  .name         = "mrc_fld",
  .size         = sizeof(struct mrc_fld),
  .param_descr  = mrc_fld_descr,
  .init         = mrc_fld_init,
  .destroy      = _mrc_fld_destroy,
  .setup        = _mrc_fld_setup,
  .write        = _mrc_fld_write,
  .read         = _mrc_fld_read,
};

// ======================================================================
// mrc_f3

int
mrc_f3_nr_comps(struct mrc_f3 *f3)
{
  return f3->_nr_comp;
}

void
mrc_f3_set_nr_comps(struct mrc_f3 *f3, int nr_comps)
{
  f3->_nr_comp = nr_comps;
}

void
mrc_f3_set_comp_name(struct mrc_f3 *f3, int m, const char *name)
{
  assert(m < f3->_nr_comp);
  if (f3->_nr_comp > f3->_nr_allocated_comp_name) {
    for (int i = 0; i < f3->_nr_allocated_comp_name; i++) {
      free(f3->_comp_name[m]);
    }
    free(f3->_comp_name);
    f3->_comp_name = calloc(f3->_nr_comp, sizeof(*f3->_comp_name));
    f3->_nr_allocated_comp_name = f3->_nr_comp;
  }
  free(f3->_comp_name[m]);
  f3->_comp_name[m] = name ? strdup(name) : NULL;
}

const char *
mrc_f3_comp_name(struct mrc_f3 *f3, int m)
{
  assert(m < f3->_nr_comp && m < f3->_nr_allocated_comp_name);
  return f3->_comp_name[m];
}

const int *
mrc_f3_off(struct mrc_f3 *f3)
{
  return f3->_offs.vals;
}

const int *
mrc_f3_dims(struct mrc_f3 *f3)
{
  return f3->_dims.vals;
}

const int *
mrc_f3_ghost_off(struct mrc_f3 *f3)
{
  return f3->_ghost_offs;
}

const int *
mrc_f3_ghost_dims(struct mrc_f3 *f3)
{
  return f3->_ghost_dims;
}

struct mrc_f3 *
mrc_f3_duplicate(struct mrc_f3 *f3)
{
  struct mrc_f3 *f3_new = mrc_f3_create(mrc_f3_comm(f3));
  mrc_f3_set_param_int3(f3_new, "offs", f3->_offs.vals);
  mrc_f3_set_param_int3(f3_new, "dims", f3->_dims.vals);
  mrc_f3_set_param_int3(f3_new, "sw", f3->_sw.vals);
  mrc_f3_set_nr_comps(f3_new, f3->_nr_comp);
  f3_new->_domain = f3->_domain;
  mrc_f3_setup(f3_new);
  return f3_new;
}

void
mrc_f3_copy(struct mrc_f3 *f3_to, struct mrc_f3 *f3_from)
{
  assert(mrc_f3_same_shape(f3_to, f3_from));

  memcpy(f3_to->_arr, f3_from->_arr, f3_to->_len * sizeof(float));
}

void
mrc_f3_set(struct mrc_f3 *f3, float val)
{
  float *arr = f3->_arr;
  for (int i = 0; i < f3->_len; i++) {
    arr[i] = val;
  }
}

void
mrc_f3_axpy(struct mrc_f3 *y, float alpha, struct mrc_f3 *x)
{
  assert(mrc_f3_same_shape(x, y));

  mrc_f3_foreach(x, ix, iy, iz, 0, 0) {
    for (int m = 0; m < x->_nr_comp; m++) {
      MRC_F3(y,m, ix,iy,iz) += alpha * MRC_F3(x,m, ix,iy,iz);
    }
  } mrc_f3_foreach_end;
}

void
mrc_f3_waxpy(struct mrc_f3 *w, float alpha, struct mrc_f3 *x, struct mrc_f3 *y)
{
  assert(mrc_f3_same_shape(x, y));
  assert(mrc_f3_same_shape(x, w));

  mrc_f3_foreach(x, ix, iy, iz, 0, 0) {
    for (int m = 0; m < x->_nr_comp; m++) {
      MRC_F3(w,m, ix,iy,iz) = alpha * MRC_F3(x,m, ix,iy,iz) + MRC_F3(y,m, ix,iy,iz);
    }
  } mrc_f3_foreach_end;
}

float
mrc_f3_norm(struct mrc_f3 *x)
{
  float res = 0.;
  mrc_f3_foreach(x, ix, iy, iz, 0, 0) {
    for (int m = 0; m < x->_nr_comp; m++) {
      res = fmaxf(res, fabsf(MRC_F3(x,m, ix,iy,iz)));
    }
  } mrc_f3_foreach_end;

  MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_FLOAT, MPI_MAX, mrc_f3_comm(x));
  return res;
}

void
mrc_f3_write_scaled(struct mrc_f3 *f3, struct mrc_io *io, float scale)
{
  mrc_io_write_ref(io, f3, "domain", f3->_domain);
  mrc_io_write_f3(io, mrc_io_obj_path(io, f3), f3, scale);
}

// ----------------------------------------------------------------------
// mrc_f3_write_comps

void
mrc_f3_write_comps(struct mrc_f3 *f3, struct mrc_io *io, int mm[])
{
  for (int i = 0; mm[i] >= 0; i++) {
    struct mrc_f3 *fld1 = mrc_f3_create(mrc_f3_comm(f3));
    mrc_f3_set_param_int3(fld1, "offs", f3->_offs.vals);
    mrc_f3_set_param_int3(fld1, "dims", f3->_dims.vals);
    mrc_f3_set_param_int3(fld1, "sw", f3->_sw.vals);
    int *ib = f3->_ghost_offs;
    mrc_f3_set_array(fld1, &MRC_F3(f3,mm[i], ib[0], ib[1], ib[2]));
    mrc_f3_set_name(fld1, f3->_comp_name[mm[i]]);
    mrc_f3_set_comp_name(fld1, 0, f3->_comp_name[mm[i]]);
    fld1->_domain = f3->_domain;
    mrc_f3_setup(fld1);
    mrc_f3_write(fld1, io);
    mrc_f3_destroy(fld1);
  }
}

static void
mrc_f3_float_create(struct mrc_fld *fld)
{
  fld->_data_type = MRC_NT_FLOAT;
  fld->_size_of_type = sizeof(float);
}

static struct mrc_f3_ops mrc_f3_float_ops = {
  .name                  = "float",
  .create                = mrc_f3_float_create,
};

// ----------------------------------------------------------------------
// mrc_f3_init

static void
mrc_f3_init()
{
  mrc_class_register_subclass(&mrc_class_mrc_f3, &mrc_f3_float_ops);
}

// ----------------------------------------------------------------------
// mrc_class_mrc_f3

#define VAR(x) (void *)offsetof(struct mrc_f3, x)
static struct param mrc_f3_params_descr[] = {
  { "offs"            , VAR(_offs)        , PARAM_INT_ARRAY(4, 0)  },
  { "dims"            , VAR(_dims)        , PARAM_INT_ARRAY(4, 0)  },
  { "sw"              , VAR(_sw)          , PARAM_INT_ARRAY(4, 0)  },
  { "nr_comps"        , VAR(_nr_comp)     , PARAM_INT(1)           },
  {},
};
#undef VAR

static struct mrc_obj_method mrc_f3_methods[] = {
  MRC_OBJ_METHOD("duplicate", mrc_f3_duplicate),
  MRC_OBJ_METHOD("copy"     , mrc_f3_copy),
  MRC_OBJ_METHOD("axpy"     , mrc_f3_axpy),
  MRC_OBJ_METHOD("waxpy"    , mrc_f3_waxpy),
  MRC_OBJ_METHOD("norm"     , mrc_f3_norm),
  {}
};

struct mrc_class_mrc_f3 mrc_class_mrc_f3 = {
  .name         = "mrc_f3",
  .size         = sizeof(struct mrc_f3),
  .param_descr  = mrc_f3_params_descr,
  .methods      = mrc_f3_methods,
  .init         = mrc_f3_init,
  .destroy      = _mrc_f3_destroy,
  .setup        = _mrc_f3_setup,
  .read         = _mrc_f3_read,
  .write        = _mrc_f3_write,
};

