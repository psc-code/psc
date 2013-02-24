
#include <mrc_fld.h>
#include <mrc_io.h>
#include <mrc_params.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// ======================================================================
// mrc_f1

static void
_mrc_f1_destroy(struct mrc_f1 *f1)
{
  if (!f1->with_array) {
    free(f1->arr);
  }
  for (int m = 0; m < f1->nr_comp; m++) {
    free(f1->_comp_name[m]);
  }
  free(f1->_comp_name);
  f1->arr = NULL;
}

static void
_mrc_f1_setup(struct mrc_f1 *f1)
{
  free(f1->_comp_name);

  f1->_comp_name = calloc(f1->nr_comp, sizeof(*f1->_comp_name));
  if (f1->domain) {
    int nr_patches;
    struct mrc_patch *patches = mrc_domain_get_patches(f1->domain, &nr_patches);
    assert(nr_patches == 1);
    f1->_off[0] = 0;
    f1->_dims[0] = patches[0].ldims[f1->dim];
  }
  f1->_ghost_off[0] = f1->_off[0] - f1->_sw;
  f1->_ghost_dims[0] = f1->_dims[0] + 2 * f1->_sw;
  f1->len = f1->_ghost_dims[0] * f1->nr_comp;

  if (!f1->arr) {
    f1->arr = calloc(f1->len, sizeof(*f1->arr));
    f1->with_array = false;
  } else {
    f1->with_array = true;
  }
}

void
mrc_f1_set_array(struct mrc_f1 *f1, float *arr)
{
  assert(!f1->arr);
  f1->arr = arr;
}

static void
_mrc_f1_read(struct mrc_f1 *f1, struct mrc_io *io)
{
  f1->domain = mrc_io_read_ref(io, f1, "domain", mrc_domain);
  
  mrc_f1_setup(f1);
  mrc_io_read_f1(io, mrc_io_obj_path(io, f1), f1);
}

static void
_mrc_f1_write(struct mrc_f1 *f1, struct mrc_io *io)
{
  mrc_io_write_ref(io, f1, "domain", f1->domain);
  mrc_io_write_f1(io, mrc_io_obj_path(io, f1), f1);
}

struct mrc_f1 *
mrc_f1_duplicate(struct mrc_f1 *f1_in)
{
  struct mrc_f1 *f1 = mrc_f1_create(mrc_f1_comm(f1_in));

  mrc_f1_set_param_int(f1, "offx", f1_in->_off[0]);
  mrc_f1_set_param_int(f1, "dimsx", f1_in->_dims[0]);
  mrc_f1_set_param_int(f1, "nr_comps", f1_in->nr_comp);
  mrc_f1_set_param_int(f1, "sw", f1_in->_sw);
  f1->domain = f1_in->domain;
  mrc_f1_setup(f1);

  return f1;
}

void
mrc_f1_set_comp_name(struct mrc_f1 *f1, int m, const char *name)
{
  assert(m < f1->nr_comp);
  free(f1->_comp_name[m]);
  f1->_comp_name[m] = name ? strdup(name) : NULL;
}

const char *
mrc_f1_comp_name(struct mrc_f1 *f1, int m)
{
  assert(m < f1->nr_comp);
  return f1->_comp_name[m];
}

const int *
mrc_f1_dims(struct mrc_f1 *f1)
{
  return f1->_dims;
}

const int *
mrc_f1_off(struct mrc_f1 *f1)
{
  return f1->_off;
}

const int *
mrc_f1_ghost_dims(struct mrc_f1 *f1)
{
  return f1->_ghost_dims;
}

void
mrc_f1_dump(struct mrc_f1 *x, const char *basename, int n)
{
  struct mrc_io *io = mrc_io_create(MPI_COMM_WORLD);
  mrc_io_set_name(io, "mrc_f1_dump");
  mrc_io_set_param_string(io, "basename", basename);
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", n, n);
  mrc_f1_write(x, io);
  mrc_io_close(io);
  mrc_io_destroy(io);
}

void
mrc_f1_zero(struct mrc_f1 *x)
{
  mrc_f1_foreach(x, ix, 0, 0) {
    for (int m = 0; m < x->nr_comp; m++) {
      MRC_F1(x,m, ix) = 0.;
    }
  } mrc_f1_foreach_end;
}

void
mrc_f1_copy(struct mrc_f1 *x, struct mrc_f1 *y)
{
  assert(x->nr_comp == y->nr_comp);
  assert(x->_ghost_dims[0] == y->_ghost_dims[0]);

  mrc_f1_foreach(x, ix, 0, 0) {
    for (int m = 0; m < y->nr_comp; m++) {
      MRC_F1(x,m, ix) = MRC_F1(y,m, ix);
    }
  } mrc_f1_foreach_end;
}

void
mrc_f1_waxpy(struct mrc_f1 *w, float alpha, struct mrc_f1 *x, struct mrc_f1 *y)
{
  assert(w->nr_comp == x->nr_comp);
  assert(w->nr_comp == y->nr_comp);
  assert(w->_ghost_dims[0] == x->_ghost_dims[0]);
  assert(w->_ghost_dims[0] == y->_ghost_dims[0]);

  mrc_f1_foreach(w, ix, 0, 0) {
    for (int m = 0; m < w->nr_comp; m++) {
      MRC_F1(w,m, ix) = alpha * MRC_F1(x,m, ix) + MRC_F1(y,m, ix);
    }
  } mrc_f1_foreach_end;
}

void
mrc_f1_axpy(struct mrc_f1 *y, float alpha, struct mrc_f1 *x)
{
  assert(x->nr_comp == y->nr_comp);
  assert(x->_ghost_dims[0] == y->_ghost_dims[0]);

  mrc_f1_foreach(x, ix, 0, 0) {
    for (int m = 0; m < y->nr_comp; m++) {
      MRC_F1(y,m, ix) += alpha * MRC_F1(x,m, ix);
    }
  } mrc_f1_foreach_end;
}

float
mrc_f1_norm(struct mrc_f1 *x)
{
  float res = 0.;
  mrc_f1_foreach(x, ix, 0, 0) {
    for (int m = 0; m < x->nr_comp; m++) {
      res = fmaxf(res, fabsf(MRC_F1(x,m, ix)));
    }
  } mrc_f1_foreach_end;

  MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_FLOAT, MPI_MAX, mrc_f1_comm(x));
  return res;
}

float
mrc_f1_norm_comp(struct mrc_f1 *x, int m)
{
  float res = 0.;
  mrc_f1_foreach(x, ix, 0, 0) {
    res = fmaxf(res, fabsf(MRC_F1(x,m, ix)));
  } mrc_f1_foreach_end;

  MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_FLOAT, MPI_MAX, mrc_f1_comm(x));
  return res;
}

// ----------------------------------------------------------------------
// mrc_class_mrc_f1

#define VAR(x) (void *)offsetof(struct mrc_f1, x)
static struct param mrc_f1_params_descr[] = {
  { "offx"            , VAR(_off[0])      , PARAM_INT(0)           },
  { "dimsx"           , VAR(_dims[0])     , PARAM_INT(0)           },
  { "nr_comps"        , VAR(nr_comp)      , PARAM_INT(1)           },
  { "sw"              , VAR(_sw)          , PARAM_INT(0)           },
  { "dim"             , VAR(dim)          , PARAM_INT(0)           },
  {},
};
#undef VAR

static struct mrc_obj_method mrc_f1_methods[] = {
  MRC_OBJ_METHOD("duplicate", mrc_f1_duplicate),
  MRC_OBJ_METHOD("copy"     , mrc_f1_copy),
  MRC_OBJ_METHOD("axpy"     , mrc_f1_axpy),
  MRC_OBJ_METHOD("waxpy"    , mrc_f1_waxpy),
  MRC_OBJ_METHOD("norm"     , mrc_f1_norm),
  {}
};

struct mrc_class_mrc_f1 mrc_class_mrc_f1 = {
  .name         = "mrc_f1",
  .size         = sizeof(struct mrc_f1),
  .param_descr  = mrc_f1_params_descr,
  .methods      = mrc_f1_methods,
  .destroy      = _mrc_f1_destroy,
  .setup        = _mrc_f1_setup,
  .read         = _mrc_f1_read,
  .write        = _mrc_f1_write,
};

// ======================================================================
// mrc_f2

void
mrc_f2_alloc(struct mrc_f2 *f2, int ib[2], int im[2], int nr_comp)
{
  memset(f2, 0, sizeof(*f2));
  f2->len = im[0] * im[1] * nr_comp;
  for (int d = 0; d < 2; d++) {
    f2->im[d] = im[d];
  }
  if (ib) {
    for (int d = 0; d < 2; d++) {
      f2->ib[d] = ib[d];
    }
  }
  f2->nr_comp = nr_comp;
  f2->arr = calloc(f2->len, sizeof(*f2->arr));
  f2->with_array = false;

  f2->name = malloc(nr_comp * sizeof(*f2->name));
  memset(f2->name, 0, nr_comp * sizeof(*f2->name));
}

void
mrc_f2_alloc_with_array(struct mrc_f2 *f2, int ib[2], int im[2], int nr_comp, float *arr)
{
  memset(f2, 0, sizeof(*f2));
  f2->len = im[0] * im[1];
  for (int d = 0; d < 2; d++) {
    f2->im[d] = im[d];
  }
  if (ib) {
    for (int d = 0; d < 2; d++) {
      f2->ib[d] = ib[d];
    }
  }
  f2->nr_comp = nr_comp;
  f2->arr = arr;
  f2->with_array = true;

  f2->name = malloc(nr_comp * sizeof(*f2->name));
  memset(f2->name, 0, nr_comp * sizeof(*f2->name));
}

void
mrc_f2_free(struct mrc_f2 *f2)
{
  if (!f2->with_array) {
    free(f2->arr);
  }
  for (int m = 0; m < f2->nr_comp; m++) {
    free(f2->name[m]);
  }
  free(f2->name);

  f2->arr = NULL;
}

// ======================================================================
// mrc_f3

#define to_mrc_f3(o) container_of(o, struct mrc_f3, obj)

static void
_mrc_f3_destroy(struct mrc_f3 *f3)
{
  if (!f3->with_array) {
    free(f3->arr);
  }
  for (int m = 0; m < f3->nr_comp; m++) {
    free(f3->_comp_name[m]);
  }
  free(f3->_comp_name);
  f3->arr = NULL;
}

static void
_mrc_f3_create(struct mrc_f3 *f3)
{
  f3->_comp_name = calloc(f3->nr_comp, sizeof(*f3->_comp_name));
}

static void
_mrc_f3_setup(struct mrc_f3 *f3)
{
  for (int d = 0; d < 3; d++) {
    f3->_ghost_off[d] = f3->_off[d] - f3->_sw;
    f3->_ghost_dims[d] = f3->_dims[d] + 2 * f3->_sw;
  }
  f3->len = f3->_ghost_dims[0] * f3->_ghost_dims[1] * f3->_ghost_dims[2] * f3->nr_comp;

  if (!f3->arr) {
    f3->arr = calloc(f3->len, sizeof(*f3->arr));
    f3->with_array = false;
  } else {
    f3->with_array = true;
  }
}

void
mrc_f3_set_nr_comps(struct mrc_f3 *f3, int nr_comps)
{
  if (nr_comps == f3->nr_comp)
    return;

  for (int m = 0; m < f3->nr_comp; m++) {
    free(f3->_comp_name[m]);
  }
  f3->nr_comp = nr_comps;
  f3->_comp_name = calloc(nr_comps, sizeof(*f3->_comp_name));
}

void
mrc_f3_set_comp_name(struct mrc_f3 *f3, int m, const char *name)
{
  assert(m < f3->nr_comp);
  free(f3->_comp_name[m]);
  f3->_comp_name[m] = name ? strdup(name) : NULL;
}

const char *
mrc_f3_comp_name(struct mrc_f3 *f3, int m)
{
  assert(m < f3->nr_comp);
  return f3->_comp_name[m];
}

void
mrc_f3_set_array(struct mrc_f3 *f3, float *arr)
{
  assert(!f3->arr);
  f3->arr = arr;
}

const int *
mrc_f3_off(struct mrc_f3 *f3)
{
  return f3->_off;
}

const int *
mrc_f3_dims(struct mrc_f3 *f3)
{
  return f3->_dims;
}

const int *
mrc_f3_ghost_off(struct mrc_f3 *f3)
{
  return f3->_ghost_off;
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
  mrc_f3_set_param_int3(f3_new, "off", f3->_off);
  mrc_f3_set_param_int3(f3_new, "dims", f3->_dims);
  mrc_f3_set_nr_comps(f3_new, f3->nr_comp);
  mrc_f3_set_param_int(f3_new, "sw", f3->_sw);
  f3_new->domain = f3->domain;
  mrc_f3_setup(f3_new);
  return f3_new;
}

void
mrc_f3_copy(struct mrc_f3 *f3_to, struct mrc_f3 *f3_from)
{
  assert(mrc_f3_same_shape(f3_to, f3_from));

  memcpy(f3_to->arr, f3_from->arr, f3_to->len * sizeof(float));
}

void
mrc_f3_set(struct mrc_f3 *f3, float val)
{
  for (int i = 0; i < f3->len; i++) {
    f3->arr[i] = val;
  }
}

void
mrc_f3_axpy(struct mrc_f3 *y, float alpha, struct mrc_f3 *x)
{
  assert(mrc_f3_same_shape(x, y));

  mrc_f3_foreach(x, ix, iy, iz, 0, 0) {
    for (int m = 0; m < x->nr_comp; m++) {
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
    for (int m = 0; m < x->nr_comp; m++) {
      MRC_F3(w,m, ix,iy,iz) = alpha * MRC_F3(x,m, ix,iy,iz) + MRC_F3(y,m, ix,iy,iz);
    }
  } mrc_f3_foreach_end;
}

float
mrc_f3_norm(struct mrc_f3 *x)
{
  float res = 0.;
  mrc_f3_foreach(x, ix, iy, iz, 0, 0) {
    for (int m = 0; m < x->nr_comp; m++) {
      res = fmaxf(res, fabsf(MRC_F3(x,m, ix,iy,iz)));
    }
  } mrc_f3_foreach_end;

  MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_FLOAT, MPI_MAX, mrc_f3_comm(x));
  return res;
}

static void
_mrc_f3_read(struct mrc_f3 *f3, struct mrc_io *io)
{
  f3->domain = mrc_io_read_ref(io, f3, "domain", mrc_domain);

  // rely on domain rather than read params
  // since the domain may be different (# of procs) when
  // we're reading things back
  // basically, we should use mrc_domain_f3_create()
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(f3->domain, &nr_patches);
  assert(nr_patches == 1);
  mrc_f3_set_param_int3(f3, "dims", patches[0].ldims);
  mrc_f3_set_param_int3(f3, "off", (int[3]) { 0, 0, 0 });
  
  mrc_f3_setup(f3);
  mrc_io_read_f3(io, mrc_io_obj_path(io, f3), f3);
}

static void
_mrc_f3_write(struct mrc_f3 *f3, struct mrc_io *io)
{
  mrc_io_write_ref(io, f3, "domain", f3->domain);
  mrc_io_write_f3(io, mrc_io_obj_path(io, f3), f3, 1.);
}

void
mrc_f3_write_scaled(struct mrc_f3 *f3, struct mrc_io *io, float scale)
{
  mrc_io_write_ref(io, f3, "domain", f3->domain);
  mrc_io_write_f3(io, mrc_io_obj_path(io, f3), f3, scale);
}

// ----------------------------------------------------------------------
// mrc_f3_write_comps

void
mrc_f3_write_comps(struct mrc_f3 *f3, struct mrc_io *io, int mm[])
{
  for (int i = 0; mm[i] >= 0; i++) {
    struct mrc_f3 *fld1 = mrc_f3_create(mrc_f3_comm(f3));
    mrc_f3_set_param_int3(fld1, "off", f3->_off);
    mrc_f3_set_param_int3(fld1, "dims", f3->_dims);
    mrc_f3_set_param_int(fld1, "sw", f3->_sw);
    int *ib = f3->_ghost_off;
    mrc_f3_set_array(fld1, &MRC_F3(f3,mm[i], ib[0], ib[1], ib[2]));
    mrc_f3_set_name(fld1, f3->_comp_name[mm[i]]);
    mrc_f3_set_comp_name(fld1, 0, f3->_comp_name[mm[i]]);
    fld1->domain = f3->domain;
    mrc_f3_setup(fld1);
    mrc_f3_write(fld1, io);
    mrc_f3_destroy(fld1);
  }
}

// ----------------------------------------------------------------------
// mrc_class_mrc_f3

#define VAR(x) (void *)offsetof(struct mrc_f3, x)
static struct param mrc_f3_params_descr[] = {
  { "off"             , VAR(_off)         , PARAM_INT3(0, 0, 0)    },
  { "dims"            , VAR(_dims)        , PARAM_INT3(0, 0, 0)    },
  { "nr_comps"        , VAR(nr_comp)      , PARAM_INT(1)           },
  { "sw"              , VAR(_sw)          , PARAM_INT(0)           },
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
  .create       = _mrc_f3_create,
  .destroy      = _mrc_f3_destroy,
  .setup        = _mrc_f3_setup,
  .read         = _mrc_f3_read,
  .write        = _mrc_f3_write,
};

// ======================================================================
// mrc_m1

#define to_mrc_m1(o) container_of(o, struct mrc_m1, obj)

static void
_mrc_m1_create(struct mrc_m1 *m1)
{
  m1->_comp_name = calloc(m1->nr_comp, sizeof(*m1->_comp_name));
}

static void
_mrc_m1_destroy(struct mrc_m1 *m1)
{
  for (int p = 0; p < m1->nr_patches; p++) {
    struct mrc_m1_patch *m1p = &m1->patches[p];
    free(m1p->arr);
  }
  free(m1->patches);

  for (int m = 0; m < m1->nr_comp; m++) {
    free(m1->_comp_name[m]);
  }
  free(m1->_comp_name);
}

static void
_mrc_m1_setup(struct mrc_m1 *m1)
{
  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(m1->domain, &nr_patches);

  m1->nr_patches = nr_patches;
  m1->patches = calloc(nr_patches, sizeof(*patches));
  for (int p = 0; p < nr_patches; p++) {
    struct mrc_m1_patch *m1p = &m1->patches[p];
    m1p->ib[0] = -m1->sw;
    m1p->im[0] = patches[p].ldims[m1->dim] + 2 * m1->sw;
    int len = m1p->im[0] * m1->nr_comp;
    m1p->arr = calloc(len, sizeof(*m1p->arr));
  }
}

static void
_mrc_m1_view(struct mrc_m1 *m1)
{
#if 0
  int rank, size;
  MPI_Comm_rank(obj->comm, &rank);
  MPI_Comm_size(obj->comm, &size);

  for (int r = 0; r < size; r++) {
    if (r == rank) {
      mrc_m1_foreach_patch(m1, p) {
	struct mrc_m1_patch *m1p = mrc_m1_patch_get(m1, p);
	mprintf("patch %d: ib = %d im = %d\n", p,
		m1p->ib[0], m1p->im[0]);
      }
    }
    MPI_Barrier(obj->comm);
  }
#endif
}

static void
_mrc_m1_write(struct mrc_m1 *m1, struct mrc_io *io)
{
  mrc_io_write_ref(io, m1, "domain", m1->domain);
  mrc_io_write_m1(io, mrc_io_obj_path(io, m1), m1);
}

static void
_mrc_m1_read(struct mrc_m1 *m1, struct mrc_io *io)
{
  m1->domain = mrc_io_read_ref(io, m1, "domain", mrc_domain);
  
  mrc_m1_setup(m1);
  mrc_io_read_m1(io, mrc_io_obj_path(io, m1), m1);
}

void
mrc_m1_set_comp_name(struct mrc_m1 *m1, int m, const char *name)
{
  assert(m < m1->nr_comp);
  free(m1->_comp_name[m]);
  m1->_comp_name[m] = name ? strdup(name) : NULL;
}

const char *
mrc_m1_comp_name(struct mrc_m1 *m1, int m)
{
  assert(m < m1->nr_comp);
  return m1->_comp_name[m];
}

bool
mrc_m1_same_shape(struct mrc_m1 *m1_1, struct mrc_m1 *m1_2)
{
  if (m1_1->nr_comp != m1_2->nr_comp) return false;
  if (m1_1->nr_patches != m1_2->nr_patches) return false;
  mrc_m1_foreach_patch(m1_1, p) {
    struct mrc_m1_patch *m1p_1 = mrc_m1_patch_get(m1_1, p);
    struct mrc_m1_patch *m1p_2 = mrc_m1_patch_get(m1_2, p);

    if (m1p_1->im[0] != m1p_2->im[0])
      return false;

    mrc_m1_patch_put(m1_1);
    mrc_m1_patch_put(m1_2);
  }
  return true;
}


// ----------------------------------------------------------------------
// mrc_class_mrc_m1

#define VAR(x) (void *)offsetof(struct mrc_m1, x)
static struct param mrc_m1_params_descr[] = {
  { "nr_comps"        , VAR(nr_comp)      , PARAM_INT(1)           },
  { "sw"              , VAR(sw)           , PARAM_INT(0)           },
  { "dim"             , VAR(dim)          , PARAM_INT(0)           },
  {},
};
#undef VAR

struct mrc_class_mrc_m1 mrc_class_mrc_m1 = {
  .name         = "mrc_m1",
  .size         = sizeof(struct mrc_m1),
  .param_descr  = mrc_m1_params_descr,
  .create       = _mrc_m1_create,
  .destroy      = _mrc_m1_destroy,
  .setup        = _mrc_m1_setup,
  .view         = _mrc_m1_view,
  .read         = _mrc_m1_read,
  .write        = _mrc_m1_write,
};

// ======================================================================
// mrc_m3

#define to_mrc_m3(o) container_of(o, struct mrc_m3, obj)

static void
_mrc_m3_destroy(struct mrc_m3 *m3)
{
  for (int p = 0; p < m3->nr_patches; p++) {
    struct mrc_m3_patch *m3p = &m3->patches[p];
    free(m3p->arr);
  }
  free(m3->patches);

  if (m3->name) {
    for (int m = 0; m < m3->nr_comp; m++) {
      free(m3->name[m]);
    }
    free(m3->name);
  }
}

static void
_mrc_m3_setup(struct mrc_m3 *m3)
{
  m3->name = calloc(m3->nr_comp, sizeof(*m3->name));

  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(m3->domain, &nr_patches);

  m3->nr_patches = nr_patches;
  m3->patches = calloc(nr_patches, sizeof(*m3->patches));
  for (int p = 0; p < nr_patches; p++) {
    struct mrc_m3_patch *m3p = &m3->patches[p];
    for (int d = 0; d < 3; d++) {
      m3p->ib[d] = -m3->sw;
      m3p->im[d] = patches[p].ldims[d] + 2 * m3->sw;
    }
    int len = m3p->im[0] * m3p->im[1] * m3p->im[2] * m3->nr_comp;
    m3p->arr = calloc(len, sizeof(*m3p->arr));
  }
}

static void
_mrc_m3_view(struct mrc_m3 *m3)
{
#if 0
  int rank, size;
  MPI_Comm_rank(obj->comm, &rank);
  MPI_Comm_size(obj->comm, &size);

  for (int r = 0; r < size; r++) {
    if (r == rank) {
      mrc_m3_foreach_patch(m3, p) {
	struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);
	mprintf("patch %d: ib = %dx%dx%d im = %dx%dx%d\n", p,
		m3p->ib[0], m3p->ib[1], m3p->ib[2],
		m3p->im[0], m3p->im[1], m3p->im[2]);
      }
    }
    MPI_Barrier(obj->comm);
  }
#endif
}

void
mrc_m3_set_comp_name(struct mrc_m3 *m3, int m, const char *name)
{
  assert(m < m3->nr_comp);
  free(m3->name[m]);
  m3->name[m] = name ? strdup(name) : NULL;
}

const char *
mrc_m3_comp_name(struct mrc_m3 *m3, int m)
{
  assert(m < m3->nr_comp);
  return m3->name[m];
}

static void
_mrc_m3_write(struct mrc_m3 *m3, struct mrc_io *io)
{
  mrc_io_write_ref(io, m3, "domain", m3->domain);
  mrc_io_write_m3(io, mrc_io_obj_path(io, m3), m3);
}

static void
_mrc_m3_read(struct mrc_m3 *m3, struct mrc_io *io)
{
  m3->domain = mrc_io_read_ref(io, m3, "domain", mrc_domain);
  mrc_m3_setup(m3);
  mrc_io_read_m3(io, mrc_io_obj_path(io, m3), m3);
}

bool
mrc_m3_same_shape(struct mrc_m3 *m3_1, struct mrc_m3 *m3_2)
{
  if (m3_1->nr_comp != m3_2->nr_comp)
    return false;

  if (m3_1->sw != m3_2->sw)
    return false;

  int nr_patches_1, nr_patches_2;
  struct mrc_patch *patches_1 = mrc_domain_get_patches(m3_1->domain, &nr_patches_1);
  struct mrc_patch *patches_2 = mrc_domain_get_patches(m3_2->domain, &nr_patches_2);
  if (nr_patches_1 != nr_patches_2)
    return false;

  mrc_m3_foreach_patch(m3_1, p) {
    for (int d = 0; d < 3; d++) {
      if (patches_1[p].ldims[d] != patches_2[p].ldims[d])
	return false;
      if (patches_1[p].off[d] != patches_2[p].off[d])
	return false;
    }
  }
  return true;
}

// ----------------------------------------------------------------------
// mrc_class_mrc_m3

#define VAR(x) (void *)offsetof(struct mrc_m3, x)
static struct param mrc_m3_params_descr[] = {
  { "nr_comps"        , VAR(nr_comp)      , PARAM_INT(1)           },
  { "sw"              , VAR(sw)           , PARAM_INT(0)           },
  {},
};
#undef VAR

struct mrc_class_mrc_m3 mrc_class_mrc_m3 = {
  .name         = "mrc_m3",
  .size         = sizeof(struct mrc_m3),
  .param_descr  = mrc_m3_params_descr,
  .destroy      = _mrc_m3_destroy,
  .setup        = _mrc_m3_setup,
  .view         = _mrc_m3_view,
  .read         = _mrc_m3_read,
  .write        = _mrc_m3_write,
};

