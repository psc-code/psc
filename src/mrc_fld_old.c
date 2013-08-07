
#include <mrc_fld.h>
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrc_vec.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

// ======================================================================
// mrc_m1

#define to_mrc_m1(o) container_of(o, struct mrc_m1, obj)

static void
_mrc_m1_destroy(struct mrc_m1 *m1)
{
  free(m1->_arr);
  free(m1->_patches);

  for (int m = 0; m < m1->_nr_allocated_comp_name; m++) {
    free(m1->_comp_name[m]);
  }
  free(m1->_comp_name);
}

static void
_mrc_m1_setup(struct mrc_m1 *m1)
{
  assert(mrc_m1_nr_comps(m1) > 0);

  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(m1->_domain, &nr_patches);

  m1->_patches = calloc(nr_patches, sizeof(*patches));
  assert(nr_patches > 0);
  assert(m1->_dims.nr_vals >= 1);
  m1->_dims.vals[0] = patches[0].ldims[m1->_dim];
  m1->_dims.vals[2] = nr_patches;

  if (m1->_offs.nr_vals == 0) {
    mrc_m1_set_param_int_array(m1, "offs", m1->_dims.nr_vals, NULL);
  }
  if (m1->_sw.nr_vals == 0) {
    mrc_m1_set_param_int_array(m1, "sw", m1->_dims.nr_vals, NULL);
  }

  m1->_len = 1;
  for (int d = 0; d < MRC_FLD_MAXDIMS; d++) {
    if (d < m1->_dims.nr_vals) {
      m1->_ghost_offs[d] = m1->_offs.vals[d] - m1->_sw.vals[d];
      m1->_ghost_dims[d] = m1->_dims.vals[d] + 2 * m1->_sw.vals[d];
    } else {
      m1->_ghost_dims[d] = 1;
    }
    m1->_len *= m1->_ghost_dims[d];
  }
  m1->_arr = calloc(m1->_len, sizeof(float));

  for (int p = 0; p < nr_patches; p++) {
    assert(patches[p].ldims[m1->_dim] == patches[0].ldims[m1->_dim]);
    struct mrc_fld_patch *m1p = &m1->_patches[p];
    m1p->_p = p;
    m1p->_fld = m1;
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
  mrc_io_write_ref(io, m1, "domain", m1->_domain);
  mrc_io_write_m1(io, mrc_io_obj_path(io, m1), m1);
}

static void
_mrc_m1_read(struct mrc_m1 *m1, struct mrc_io *io)
{
  m1->_domain = mrc_io_read_ref(io, m1, "domain", mrc_domain);
  
  mrc_m1_setup(m1);
  mrc_io_read_m1(io, mrc_io_obj_path(io, m1), m1);
}

void
mrc_m1_set_comp_name(struct mrc_m1 *fld, int m, const char *name)
{
  int nr_comps = mrc_m1_nr_comps(fld);
  assert(m < nr_comps);
  if (nr_comps > fld->_nr_allocated_comp_name) {
    for (int i = 0; i < fld->_nr_allocated_comp_name; i++) {
      free(fld->_comp_name[m]);
    }
    free(fld->_comp_name);
    fld->_comp_name = calloc(nr_comps, sizeof(*fld->_comp_name));
    fld->_nr_allocated_comp_name = nr_comps;
  }
  free(fld->_comp_name[m]);
  fld->_comp_name[m] = name ? strdup(name) : NULL;
}

const char *
mrc_m1_comp_name(struct mrc_m1 *fld, int m)
{
  assert(m < mrc_m1_nr_comps(fld) && m < fld->_nr_allocated_comp_name);
  return fld->_comp_name[m];
}

void
mrc_m1_set_sw(struct mrc_m1 *fld, int sw)
{
  assert(fld->_domain);
  if (fld->_dims.nr_vals == 3) {
    mrc_m1_set_param_int_array(fld, "sw", fld->_dims.nr_vals,
			       (int[3]) { sw, 0, 0 });
  } else {
    assert(0);
  }
}

bool
mrc_m1_same_shape(struct mrc_m1 *m1_1, struct mrc_m1 *m1_2)
{
  if (mrc_m1_nr_comps(m1_1) != mrc_m1_nr_comps(m1_2)) return false;
  if (mrc_m1_nr_patches(m1_1) != mrc_m1_nr_patches(m1_2)) return false;
  if (m1_1->_ghost_dims[0] != m1_2->_ghost_dims[0]) return false;

  return true;
}

const int *
mrc_m1_dims(struct mrc_m1 *x)
{
  return x->_dims.vals;
}

const int *
mrc_m1_ghost_offs(struct mrc_m1 *x)
{
  return x->_ghost_offs;
}

const int *
mrc_m1_ghost_dims(struct mrc_m1 *x)
{
  return x->_ghost_dims;
}

static int
mrc_m1_comp_dim(struct mrc_m1 *fld)
{
  if (fld->_domain) {
    if (fld->_dims.nr_vals == 3) {
      // emulating mrc_f3, mrc_m3
      return 1;
    } else {
      assert(0);
    }
  }
  assert(0);
}

int
mrc_m1_nr_comps(struct mrc_m1 *fld)
{
  int comp_dim = mrc_m1_comp_dim(fld);
  assert(comp_dim < fld->_dims.nr_vals);
  return fld->_dims.vals[comp_dim];
}

void
mrc_m1_set_nr_comps(struct mrc_m1 *fld, int nr_comps)
{
  int comp_dim = mrc_m1_comp_dim(fld);
  fld->_dims.vals[comp_dim] = nr_comps;
}

int
mrc_m1_nr_patches(struct mrc_m1 *fld)
{
  assert(fld->_domain);
  assert(fld->_dims.nr_vals == 3);
  return fld->_ghost_dims[2];
}

// ----------------------------------------------------------------------
// mrc_class_mrc_m1

#define VAR(x) (void *)offsetof(struct mrc_m1, x)
static struct param mrc_m1_params_descr[] = {
  { "offs"            , VAR(_offs)        , PARAM_INT_ARRAY(0, 0) },
  { "dims"            , VAR(_dims)        , PARAM_INT_ARRAY(0, 0) },
  { "sw"              , VAR(_sw)          , PARAM_INT_ARRAY(0, 0) },

  { "dim"             , VAR(_dim)         , PARAM_INT(0)           },
  {},
};
#undef VAR

struct mrc_class_mrc_m1 mrc_class_mrc_m1 = {
  .name         = "mrc_m1",
  .size         = sizeof(struct mrc_m1),
  .param_descr  = mrc_m1_params_descr,
  .destroy      = _mrc_m1_destroy,
  .setup        = _mrc_m1_setup,
  .view         = _mrc_m1_view,
  .read         = _mrc_m1_read,
  .write        = _mrc_m1_write,
};


