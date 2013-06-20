
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
  f1->arr = NULL;

  if (f1->_comp_name) {
    for (int m = 0; m < f1->nr_comp; m++) {
      free(f1->_comp_name[m]);
    }
    free(f1->_comp_name);
  }
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
  assert(f1->_comp_name);
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
  
  m1->_comp_name = calloc(m1->nr_comp, sizeof(*m1->_comp_name));
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
_mrc_m3_create(struct mrc_m3 *m3)
{
  m3->_data_type = MRC_NT_FLOAT;
  m3->_size_of_type = sizeof(float);
}

static void
_mrc_m3_destroy(struct mrc_m3 *m3)
{
  free(m3->_arr);
  free(m3->patches);

  if (m3->_comp_name) {
    for (int m = 0; m < mrc_m3_nr_comps(m3); m++) {
      free(m3->_comp_name[m]);
    }
    free(m3->_comp_name);
  }
}

static void
_mrc_m3_setup(struct mrc_m3 *m3)
{
  if (m3->_sw.nr_vals == 0) {
    mrc_m3_set_param_int_array(m3, "sw", 5, NULL);
  }

  m3->_comp_name = calloc(m3->_ghost_dims[3], sizeof(*m3->_comp_name));

  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(m3->domain, &nr_patches);

  assert(nr_patches > 0);
  for (int d = 0; d < 3; d++) {
    m3->_ghost_offs[d] = -m3->_sw.vals[d];
    m3->_ghost_dims[d] = patches[0].ldims[d] + 2 * m3->_sw.vals[d];
  }
  m3->_ghost_dims[4] = nr_patches;
  int len = m3->_ghost_dims[0] * m3->_ghost_dims[1] * m3->_ghost_dims[2] * m3->_ghost_dims[3] * m3->_ghost_dims[4];
  m3->_arr = calloc(len, m3->_size_of_type);

  m3->patches = calloc(nr_patches, sizeof(*m3->patches));
  for (int p = 0; p < nr_patches; p++) {
    struct mrc_m3_patch *m3p = &m3->patches[p];
    m3p->_m3 = m3;
    m3p->_p = p;
    for (int d = 0; d < 3; d++) {
      assert(m3->_ghost_dims[d] = patches[p].ldims[d] + 2 * m3->_sw.vals[d]);
    }
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

int
mrc_m3_nr_patches(struct mrc_m3 *m3)
{
  return m3->_ghost_dims[4];
}

void
mrc_m3_set_sw(struct mrc_m3 *m3, int sw)
{
  assert(m3->domain);
  mrc_m3_set_param_int_array(m3, "sw", 5, (int[5]) { sw, sw, sw, 0, 0 });
}

void
mrc_m3_set_nr_comps(struct mrc_m3 *m3, int nr_comps)
{
  m3->_ghost_dims[3] = nr_comps;
}

int
mrc_m3_nr_comps(struct mrc_m3 *m3)
{
  return m3->_ghost_dims[3];
}

void
mrc_m3_set_comp_name(struct mrc_m3 *m3, int m, const char *name)
{
  assert(m < mrc_m3_nr_comps(m3));
  free(m3->_comp_name[m]);
  m3->_comp_name[m] = name ? strdup(name) : NULL;
}

const char *
mrc_m3_comp_name(struct mrc_m3 *m3, int m)
{
  assert(m < mrc_m3_nr_comps(m3));
  return m3->_comp_name[m];
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
  if (mrc_m3_nr_comps(m3_1) != mrc_m3_nr_comps(m3_2))
    return false;

  if (m3_1->_sw.vals[0] != m3_2->_sw.vals[0])
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
  { "sw"              , VAR(_sw)           , PARAM_INT_ARRAY(0, 0)     },
  {},
};
#undef VAR

struct mrc_class_mrc_m3 mrc_class_mrc_m3 = {
  .name         = "mrc_m3",
  .size         = sizeof(struct mrc_m3),
  .param_descr  = mrc_m3_params_descr,
  .create       = _mrc_m3_create,
  .destroy      = _mrc_m3_destroy,
  .setup        = _mrc_m3_setup,
  .view         = _mrc_m3_view,
  .read         = _mrc_m3_read,
  .write        = _mrc_m3_write,
};

