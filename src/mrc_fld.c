
#include <mrc_fld.h>
#include <mrc_io.h>
#include <mrc_params.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>

// ======================================================================
// mrc_f1

#define to_mrc_f1(o) container_of(o, struct mrc_f1, obj)

static void
_mrc_f1_destroy(struct mrc_obj *obj)
{
  struct mrc_f1 *f1 = to_mrc_f1(obj);

  if (!f1->with_array) {
    free(f1->arr);
  }
  for (int m = 0; m < f1->nr_comp; m++) {
    free(f1->name[m]);
  }
  free(f1->name);
  f1->arr = NULL;
}

static void
_mrc_f1_create(struct mrc_obj *obj)
{
  struct mrc_f1 *f1 = to_mrc_f1(obj);

  f1->name = calloc(f1->nr_comp, sizeof(*f1->name));
}

static void
_mrc_f1_setup(struct mrc_obj *obj)
{
  struct mrc_f1 *f1 = to_mrc_f1(obj);
  
  f1->len = f1->im[0] * f1->nr_comp;

  if (!f1->arr) {
    f1->arr = calloc(f1->len, sizeof(*f1->arr));
    f1->with_array = false;
  } else {
    f1->with_array = true;
  }
}

static void
_mrc_f1_read(struct mrc_obj *obj, struct mrc_io *io)
{
  struct mrc_f1 *f1 = to_mrc_f1(obj);

  mrc_f1_setup(f1);
  mrc_io_read_f1(io, mrc_obj_name(obj), f1);
}

static void
_mrc_f1_write(struct mrc_obj *obj, struct mrc_io *io)
{
  struct mrc_f1 *f1 = to_mrc_f1(obj);

  mrc_io_write_f1(io, mrc_obj_name(obj), f1);
}

// ----------------------------------------------------------------------
// mrc_class_mrc_f1

#define VAR(x) (void *)offsetof(struct mrc_f1, x)
static struct param mrc_f1_params_descr[] = {
  { "ibx"             , VAR(ib[0])        , PARAM_INT(0)           },
  { "imx"             , VAR(im[0])        , PARAM_INT(0)           },

  { "nr_comps"        , VAR(nr_comp)      , PARAM_INT(1)           },
  {},
};
#undef VAR

struct mrc_class mrc_class_mrc_f1 = {
  .name         = "mrc_f1",
  .size         = sizeof(struct mrc_f1),
  .param_descr  = mrc_f1_params_descr,
  .create       = _mrc_f1_create,
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
_mrc_f3_destroy(struct mrc_obj *obj)
{
  struct mrc_f3 *f3 = to_mrc_f3(obj);

  if (!f3->with_array) {
    free(f3->arr);
  }
  for (int m = 0; m < f3->nr_comp; m++) {
    free(f3->name[m]);
  }
  free(f3->name);
  f3->arr = NULL;
}

static void
_mrc_f3_create(struct mrc_obj *obj)
{
  struct mrc_f3 *f3 = to_mrc_f3(obj);

  f3->name = calloc(f3->nr_comp, sizeof(*f3->name));
}

struct mrc_f3 *
mrc_f3_alloc(MPI_Comm comm, int ib[3], int im[3])
{
  struct mrc_f3 *f3 = mrc_f3_create(comm);
  for (int d = 0; d < 3; d++) {
    f3->ib[d] = ib ? ib[d] : 0;
    f3->im[d] = im[d];
  }

  return f3;
}

struct mrc_f3 *
mrc_f3_alloc_with_array(MPI_Comm comm, int ib[3], int im[3], float *arr)
{
  struct mrc_f3 *f3 = mrc_f3_create(comm);
  for (int d = 0; d < 3; d++) {
    f3->ib[d] = ib ? ib[d] : 0;
    f3->im[d] = im[d];
  }

  f3->arr = arr;
  return f3;
}

static void
_mrc_f3_setup(struct mrc_obj *obj)
{
  struct mrc_f3 *f3 = to_mrc_f3(obj);
  
  f3->len = f3->im[0] * f3->im[1] * f3->im[2] * f3->nr_comp;

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
    free(f3->name[m]);
  }
  f3->nr_comp = nr_comps;
  f3->name = calloc(nr_comps, sizeof(*f3->name));
}

void
mrc_f3_set_array(struct mrc_f3 *f3, float *arr)
{
  assert(!f3->arr);
  f3->arr = arr;
}

struct mrc_f3 *
mrc_f3_duplicate(struct mrc_f3 *f3)
{
  struct mrc_f3 *f3_new = mrc_f3_alloc(f3->obj.comm, f3->ib, f3->im);
  mrc_f3_set_nr_comps(f3_new, f3->nr_comp);
  mrc_f3_setup(f3_new); // FIXME
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

static void
_mrc_f3_read(struct mrc_obj *obj, struct mrc_io *io)
{
  struct mrc_f3 *f3 = to_mrc_f3(obj);

  mrc_io_read_attr_int(io, mrc_obj_name(obj), "sw", &f3->sw);
  f3->domain = (struct mrc_domain *)
    mrc_io_read_obj_ref(io, mrc_obj_name(obj), "domain", &mrc_class_mrc_domain);
  
  mrc_f3_setup(f3);
  mrc_io_read_f3(io, mrc_obj_name(obj), f3);
}

static void
_mrc_f3_write(struct mrc_obj *obj, struct mrc_io *io)
{
  struct mrc_f3 *f3 = to_mrc_f3(obj);

  mrc_io_write_attr_int(io, mrc_obj_name(obj), "sw", f3->sw);
  mrc_io_write_obj_ref(io, mrc_obj_name(obj), "domain",
		       (struct mrc_obj *) f3->domain);
  mrc_io_write_f3(io, mrc_obj_name(obj), f3, 1.);
}

void
mrc_f3_write_scaled(struct mrc_f3 *f3, struct mrc_io *io, float scale)
{
  mrc_io_write_f3(io, mrc_f3_name(f3), f3, scale);
}

// ----------------------------------------------------------------------
// mrc_f3_write_comps

void
mrc_f3_write_comps(struct mrc_f3 *f3, struct mrc_io *io, int mm[])
{
  for (int i = 0; mm[i] >= 0; i++) {
    int *ib = f3->ib, *im = f3->im;
    struct mrc_f3 *fld1 = mrc_f3_alloc_with_array(f3->obj.comm, ib, im,
						  &MRC_F3(f3,mm[i], ib[0], ib[1], ib[2]));
    mrc_f3_set_name(fld1, f3->name[mm[i]]);
    fld1->name[0] = strdup(f3->name[mm[i]]);
    fld1->domain = f3->domain;
    fld1->sw = f3->sw;
    mrc_f3_setup(fld1);
    mrc_f3_write(fld1, io);
    mrc_f3_destroy(fld1);
  }
}

// ----------------------------------------------------------------------
// mrc_class_mrc_f3

#define VAR(x) (void *)offsetof(struct mrc_f3, x)
static struct param mrc_f3_params_descr[] = {
  { "ibx"             , VAR(ib[0])        , PARAM_INT(0)           },
  { "iby"             , VAR(ib[1])        , PARAM_INT(0)           },
  { "ibz"             , VAR(ib[2])        , PARAM_INT(0)           },
  { "imx"             , VAR(im[0])        , PARAM_INT(0)           },
  { "imy"             , VAR(im[1])        , PARAM_INT(0)           },
  { "imz"             , VAR(im[2])        , PARAM_INT(0)           },

  { "nr_comps"        , VAR(nr_comp)      , PARAM_INT(1)           },
  {},
};
#undef VAR

struct mrc_class mrc_class_mrc_f3 = {
  .name         = "mrc_f3",
  .size         = sizeof(struct mrc_f3),
  .param_descr  = mrc_f3_params_descr,
  .create       = _mrc_f3_create,
  .destroy      = _mrc_f3_destroy,
  .setup        = _mrc_f3_setup,
  .read         = _mrc_f3_read,
  .write        = _mrc_f3_write,
};

