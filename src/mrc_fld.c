
#include <mrc_fld.h>

#include <mrc_vec.h>
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrc_profile.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <execinfo.h>

// ======================================================================
// mrc_fld

#define mrc_fld_ops(fld) ((struct mrc_fld_ops *) (fld)->obj.ops)

// ----------------------------------------------------------------------
// mrc_fld_data_type

int
mrc_fld_data_type(struct mrc_fld *fld)
{
  return fld->_nd->data_type;
}

// ----------------------------------------------------------------------
// mrc_fld_len

int
mrc_fld_len(struct mrc_fld *fld)
{
  return fld->_nd->len;
}

// ----------------------------------------------------------------------
// mrc_fld_create

static void
_mrc_fld_create(struct mrc_fld *fld)
{
  fld->_nd = mrc_ndarray_create(mrc_fld_comm(fld));
}

// ----------------------------------------------------------------------
// mrc_fld_destroy

static void
_mrc_fld_destroy(struct mrc_fld *fld)
{
  for (int m = 0; m < fld->_nr_allocated_comp_name; m++) {
    free(fld->_comp_name[m]);
  }

  if (fld->_comp_name) {
    free(fld->_comp_name);
    fld->_comp_name = NULL;
  }

  if (fld->_patches) {
    free(fld->_patches);
    fld->_patches = NULL;
  }

  mrc_ndarray_destroy(fld->_nd);
}

// ----------------------------------------------------------------------
// setup_ghost_dims_offs

static void
setup_ghost_dims_offs(struct mrc_fld *fld)
{
  int n_dims = fld->_dims.nr_vals;
  assert(n_dims == fld->_offs.nr_vals);
  assert(n_dims == fld->_sw.nr_vals);
  assert(n_dims <= MRC_FLD_MAXDIMS);

  for (int d = 0; d < MRC_FLD_MAXDIMS; d++) {
    if (d < n_dims) {
      fld->_ghost_offs[d] = fld->_offs.vals[d] - fld->_sw.vals[d];
      fld->_ghost_dims[d] = fld->_dims.vals[d] + 2 * fld->_sw.vals[d];
    } else {
      fld->_ghost_dims[d] = 1;
    }
  }
}

// ----------------------------------------------------------------------
// setup_ndarray

static void
setup_ndarray(struct mrc_fld *fld)
{
  setup_ghost_dims_offs(fld);
  
  assert(MRC_FLD_MAXDIMS == 5);
  int *perm;
  if (!fld->_view_base) {
    // This is a field with its own storage

    if (!fld->_aos && !fld->_c_order) {
      // original order
      perm = (int[5]) { 0, 1, 2, 3, 4 };
    } else if (fld->_aos && !fld->_c_order) {
      // x,y,z,m,p to m,x,y,z,p
      perm = (int[5]) { 3, 0, 1, 2, 4 };
    } else if (fld->_aos && fld->_c_order) {
      // x,y,z,m,p to m,z,y,x,p
      perm = (int[5]) { 3, 2, 1, 0, 4 };
    } else if (!fld->_aos && fld->_c_order) {
      // x,y,z,m,p to z,y,x,m,p
      perm = (int[5]) { 2, 1, 0, 3, 4 };
    }
  } else {
    perm = (int[5]) { 0, 1, 2, 3, 4 };
  }
      
  struct mrc_ndarray *nd = fld->_nd;
  int n_dims = fld->_dims.nr_vals;
  mrc_ndarray_set_type(nd, mrc_fld_ops(fld)->vec_type);
  mrc_ndarray_set_param_int_array(nd, "dims", n_dims, fld->_ghost_dims);
  mrc_ndarray_set_param_int_array(nd, "offs", n_dims, fld->_ghost_offs);
  mrc_ndarray_set_param_int_array(nd, "perm", n_dims, perm);

  if (fld->_view_base) {
    // In this case, this field is just a view and has not
    // allocated its own storage
    struct mrc_fld *view_base = fld->_view_base;
    struct mrc_ndarray *nd_base = view_base->_nd;
    assert(nd_base->n_dims == n_dims);
    int view_ghost_offs[MRC_FLD_MAXDIMS];
    for (int d = 0; d < n_dims; d++) {
      view_ghost_offs[d] = fld->_view_offs[d] - fld->_sw.vals[d];
    }
    mrc_ndarray_set_param_obj(nd, "view_base", view_base->_nd);
    mrc_ndarray_set_param_int_array(nd, "view_offs", n_dims, view_ghost_offs);
  }

  mrc_ndarray_setup(nd);
  fld->nd_acc = *mrc_ndarray_access(nd);
}


// ----------------------------------------------------------------------
// mrc_fld_setup

static void
setup_3d_from_domain(struct mrc_fld *fld, int *ldims, int nr_patches, struct mrc_patch *patches)
{
  assert(fld->_dim < 0);
  fld->_patches = calloc(nr_patches, sizeof(*fld->_patches));
  for (int p = 0; p < nr_patches; p++) {
    struct mrc_fld_patch *m3p = &fld->_patches[p];
    m3p->_fld = fld;
    m3p->_p = p;
    for (int d = 0; d < 3; d++) {
      assert(patches[p].ldims[d] == ldims[d]);
    }
  }
  
  int sw[3];
  for (int d = 0; d < 3; d++) {
    if (ldims[d] > 1) {
      sw[d] = fld->_nr_ghosts;
    } else {
      sw[d] = 0;
    }
  }
  
  mrc_fld_set_param_int_array(fld, "dims", 5,
      (int[5]) { ldims[0], ldims[1], ldims[2], fld->_nr_comps, nr_patches });
  mrc_fld_set_param_int_array(fld, "sw", fld->_dims.nr_vals,
      (int[5]) { sw[0], sw[1], sw[2], 0, 0 });
  
}

static void
setup_1d_from_domain(struct mrc_fld *fld, int *ldims, int nr_patches, struct mrc_patch *patches)
{
      assert(fld->_dim >= 0);

      fld->_patches = calloc(nr_patches, sizeof(*fld->_patches));
      for (int p = 0; p < nr_patches; p++) {
	struct mrc_fld_patch *m1p = &fld->_patches[p];
	m1p->_p = p;
	m1p->_fld = fld;
	assert(patches[p].ldims[fld->_dim] == ldims[fld->_dim]);
      }

      int sw = fld->_nr_ghosts;
      mrc_fld_set_param_int_array(fld, "dims", 3,
				  (int[3]) { ldims[fld->_dim], fld->_nr_comps, nr_patches });
      mrc_fld_set_param_int_array(fld, "sw", 3,
				  (int[3]) { sw, 0, 0 });
}

static void
_mrc_fld_setup(struct mrc_fld *fld)
{
  if (fld->_domain) {
    // if we have a domain, use that to set _dims, _sw
    
    // if "domain" is set, can't set "dims", too, which will be set automatically
    // based on "domain"
    assert(fld->_dims.nr_vals == 0);
    assert(fld->_nr_comps > 0);
    int nr_patches;
    struct mrc_patch *patches = mrc_domain_get_patches(fld->_domain, &nr_patches);
    assert(nr_patches > 0);
    int *ldims = patches[0].ldims;

    assert(fld->_nr_spatial_dims == 1 || fld->_nr_spatial_dims == 3);

    if (fld->_nr_spatial_dims == 3) {
      setup_3d_from_domain(fld, ldims, nr_patches, patches);
      
    } else if (fld->_nr_spatial_dims == 1) {
      setup_1d_from_domain(fld, ldims, nr_patches, patches);
    }
  } else {
    assert(fld->_nr_spatial_dims < 0);
    assert(fld->_nr_ghosts == 0);
  }

  if (fld->_offs.nr_vals == 0) {
    mrc_fld_set_param_int_array(fld, "offs", fld->_dims.nr_vals, NULL);
  }
  if (fld->_sw.nr_vals == 0) {
    mrc_fld_set_param_int_array(fld, "sw", fld->_dims.nr_vals, NULL);
  }

  setup_ndarray(fld);
}

// ----------------------------------------------------------------------
// mrc_fld_set_array

void
mrc_fld_set_array(struct mrc_fld *fld, void *arr)
{
  mrc_ndarray_set_array(fld->_nd, arr);
}

// ----------------------------------------------------------------------
// mrc_fld_replace_array

void
mrc_fld_replace_array(struct mrc_fld *fld, void *arr)
{
  mrc_ndarray_replace_array(fld->_nd, arr);
}

// ----------------------------------------------------------------------
// mrc_fld_write

static void
_mrc_fld_write(struct mrc_fld *fld, struct mrc_io *io)
{
  mrc_io_write_fld(io, mrc_io_obj_path(io, fld), fld);
}

// ----------------------------------------------------------------------
// mrc_fld_read

static void
_mrc_fld_read(struct mrc_fld *fld, struct mrc_io *io)
{
  // If the fld has a domain and we don't run a minimal setup
  // the metadata can't be trusted.
  // This will overwrite the obsolete read-in metadata

  // FIXME, this looks rather similar to _setup(), but also
  // there's the question why we shouldn't trust what we're reading?
  // (I guess if were to read on a different number of procs, but that's
  //  hard work)
  if (fld->_domain) {
    if (fld->_patches) free(fld->_patches);
    int nr_patches;
    struct mrc_patch *patches = mrc_domain_get_patches(fld->_domain, &nr_patches);
    assert(nr_patches > 0);
    int *ldims = patches[0].ldims;

    if (fld->_nr_spatial_dims == 3) {
      setup_3d_from_domain(fld, ldims, nr_patches, patches);
      
    } else if (fld->_nr_spatial_dims == 1) {
      setup_1d_from_domain(fld, ldims, nr_patches, patches);
    }
  }

  if (fld->_offs.nr_vals == 0) {
    mrc_fld_set_param_int_array(fld, "offs", fld->_dims.nr_vals, NULL);
  }
  if (fld->_sw.nr_vals == 0) {
    mrc_fld_set_param_int_array(fld, "sw", fld->_dims.nr_vals, NULL);
  }
  // we're not reading back ndarray, rather we're creating a new one
  fld->_nd = mrc_ndarray_create(mrc_fld_comm(fld));
  setup_ndarray(fld);
  fld->nd_acc = *mrc_ndarray_access(fld->_nd);

  setup_ghost_dims_offs(fld);
  // FIXME: Hacky, but we are basically set up now, so we should advertise it
  fld->obj.is_setup = true;

  // Now we can actually read the data
  mrc_io_read_fld(io, mrc_io_obj_path(io, fld), fld);
}

// ----------------------------------------------------------------------
// mrc_fld_set_comp_name

void
mrc_fld_set_comp_name(struct mrc_fld *fld, int m, const char *name)
{
  int nr_comps = mrc_fld_nr_comps(fld);
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

// ----------------------------------------------------------------------
// mrc_fld_comp_name

const char *
mrc_fld_comp_name(struct mrc_fld *fld, int m)
{
  assert(m < mrc_fld_nr_comps(fld));
  if (m < fld->_nr_allocated_comp_name) {
    return fld->_comp_name[m];
  } else {
    return NULL;
  }
}

// ----------------------------------------------------------------------
// mrc_fld_nr_comps

int
mrc_fld_nr_comps(struct mrc_fld *fld)
{
  return fld->_nr_comps;
}

// ----------------------------------------------------------------------
// mrc_fld_set_comp_names
//
// sets all component names from one :-separated string

void
mrc_fld_set_comp_names(struct mrc_fld *fld, const char *comps)
{
  assert(comps);
  char *s1, *s = strdup(comps), *s_save = s;
  // count nr of components first
  int nr_comps = 0;
  while (strsep(&s, ",:")) {
    nr_comps++;
  }
  fld->_nr_comps = nr_comps;

  s = s_save;
  strcpy(s, comps);
  // then parse the names
  for (int m = 0; (s1 = strsep(&s, ",:")); m++) {
    mrc_fld_set_comp_name(fld, m, s1);
  }
  free(s_save);
}

// ----------------------------------------------------------------------
// mrc_fld_nr_patches
//
// returns the number of patches on the this processor that comprise the
// local domain (only works on former "mrc_m1/m3" for now)

int
mrc_fld_nr_patches(struct mrc_fld *fld)
{
  assert(fld->_domain);
  assert(fld->_dims.nr_vals == 3 || fld->_dims.nr_vals == 5);
  return fld->_ghost_dims[fld->_dims.nr_vals - 1];
}

// ----------------------------------------------------------------------
// mrc_fld_offs

const int *
mrc_fld_offs(struct mrc_fld *fld)
{
  return fld->_offs.vals;
}

// ----------------------------------------------------------------------
// mrc_fld_dims

const int *
mrc_fld_dims(struct mrc_fld *fld)
{
  return fld->_dims.vals;
}

// ----------------------------------------------------------------------
// mrc_fld_sw

const int *
mrc_fld_sw(struct mrc_fld *fld)
{
  return fld->_sw.vals;
}

// ----------------------------------------------------------------------
// mrc_fld_ghost_offs

const int *
mrc_fld_ghost_offs(struct mrc_fld *fld)
{
  return fld->_ghost_offs;
}

// ----------------------------------------------------------------------
// mrc_fld_ghost_dims

const int *
mrc_fld_ghost_dims(struct mrc_fld *fld)
{
  return fld->_ghost_dims;
}

// ----------------------------------------------------------------------
// mrc_fld_duplicate

struct mrc_fld *
mrc_fld_duplicate(struct mrc_fld *fld)
{
  struct mrc_fld *fld_new = mrc_fld_create(mrc_fld_comm(fld));
  mrc_fld_set_type(fld_new, mrc_fld_type(fld));
  if (fld->_domain) {
    // if we're based on a domain, dims/offs/sw will be set by setup()
    mrc_fld_set_param_obj(fld_new, "domain", fld->_domain);
    mrc_fld_set_param_int(fld_new, "nr_spatial_dims", fld->_nr_spatial_dims);
    mrc_fld_set_param_int(fld_new, "nr_comps", fld->_nr_comps);
    mrc_fld_set_param_int(fld_new, "nr_ghosts", fld->_nr_ghosts);
    mrc_fld_set_param_int(fld_new, "dim", fld->_dim);
  } else {
    mrc_fld_set_param_int_array(fld_new, "dims", fld->_dims.nr_vals, fld->_dims.vals);
    mrc_fld_set_param_int_array(fld_new, "offs", fld->_offs.nr_vals, fld->_offs.vals);
    mrc_fld_set_param_int_array(fld_new, "sw", fld->_sw.nr_vals, fld->_sw.vals);
  }
  mrc_fld_set_param_bool(fld_new, "aos", fld->_aos);
  mrc_fld_set_param_bool(fld_new, "c_order", fld->_c_order);
  mrc_fld_setup(fld_new);
  return fld_new;
}

// ----------------------------------------------------------------------
// mrc_fld_copy

void
mrc_fld_copy(struct mrc_fld *fld_to, struct mrc_fld *fld_from)
{
  mrc_ndarray_copy(fld_to->_nd, fld_from->_nd);
}


// ======================================================================
// fld math operations
// FIXME: I don't really understand why these math operations are methods...

// ----------------------------------------------------------------------
// mrc_fld_axpy

void
mrc_fld_axpy(struct mrc_fld *y, float alpha, struct mrc_fld *x)
{
  assert(mrc_fld_same_shape(x, y));
  mrc_vec_axpy(y->_nd->vec, (double) alpha, x->_nd->vec);
}

// ----------------------------------------------------------------------
// mrc_fld_waxpy

// FIXME, should take double alpha
void
mrc_fld_waxpy(struct mrc_fld *w, float alpha, struct mrc_fld *x, struct mrc_fld *y)
{
  assert(mrc_fld_same_shape(x, y));
  assert(mrc_fld_same_shape(x, w));
  mrc_vec_waxpy(w->_nd->vec, (double) alpha, x->_nd->vec, y->_nd->vec);
}

// ----------------------------------------------------------------------
// mrc_fld_axpby
//
// FIXME, does interior only, should have some way to specify
// should use mrc_vec*

void
mrc_fld_axpby(struct mrc_fld *y, double a, struct mrc_fld *x, double b)
{
  if (a == 0. && b == 1.) { // nothing to do
    return;
  }

  assert(mrc_fld_same_shape(x, y));
  mrc_vec_axpby(y->_nd->vec, a, x->_nd->vec, b);
}

// ----------------------------------------------------------------------
// mrc_fld_norm

// FIXME, should return double
float
mrc_fld_norm(struct mrc_fld *fld)
{
  double norm = mrc_ndarray_norm(fld->_nd);

  MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_FLOAT, MPI_MAX, mrc_fld_comm(fld));
  return norm;
}

// ----------------------------------------------------------------------
// mrc_fld_norm_comp

// FIXME, should go away, use view instead
float
mrc_fld_norm_comp(struct mrc_fld *x, int m)
{
  assert(mrc_fld_data_type(x) == MRC_NT_FLOAT);
  float res = 0.;
  if (x->_dims.nr_vals == 3) {
    mrc_m1_foreach_patch(x, p) {
      mrc_m1_foreach(x, ix, 0, 0) {
	res = fmaxf(res, fabsf(MRC_M1(x,m, ix, p)));
      } mrc_m1_foreach_end;
    }
  } else {
    assert(0);
  }
  MPI_Allreduce(MPI_IN_PLACE, &res, 1, MPI_FLOAT, MPI_MAX, mrc_fld_comm(x));
  return res;
}


// ----------------------------------------------------------------------
// mrc_fld_set

void
mrc_fld_set(struct mrc_fld *fld, float val)
{
  mrc_ndarray_set(fld->_nd, val);
}

// ----------------------------------------------------------------------
// mrc_fld_write_comps

void
mrc_fld_write_comps(struct mrc_fld *fld, struct mrc_io *io, int mm[])
{
  for (int i = 0; mm[i] >= 0; i++) {
    struct mrc_fld *fld1 = mrc_fld_create(mrc_fld_comm(fld));
    mrc_fld_set_param_int_array(fld1, "offs", fld->_offs.nr_vals, fld->_offs.vals);
    int *dims = fld->_dims.vals;
    mrc_fld_set_param_int_array(fld1, "dims", 5, (int[5]) { dims[0], dims[1], dims[2], 1, 1 });
    mrc_fld_set_param_int_array(fld1, "sw", fld->_sw.nr_vals, fld->_sw.vals);
    int *ib = fld->_ghost_offs;
    mrc_fld_set_array(fld1, &MRC_F3(fld,mm[i], ib[0], ib[1], ib[2]));
    mrc_fld_set_name(fld1, fld->_comp_name[mm[i]]);
    mrc_fld_set_comp_name(fld1, 0, fld->_comp_name[mm[i]]);
    mrc_fld_set_param_obj(fld1, "domain", fld->_domain);
    mrc_fld_setup(fld1);
    mrc_fld_write(fld1, io);
    mrc_fld_destroy(fld1);
  }
}

// ----------------------------------------------------------------------
// mrc_fld_dump

void
mrc_fld_dump(struct mrc_fld *x, const char *basename, int n)
{
  struct mrc_io *io = mrc_io_create(mrc_fld_comm(x));
  mrc_io_set_name(io, "mrc_fld_dump");
  mrc_io_set_param_string(io, "basename", basename);
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", n, n);
  for (int m = 0; m < mrc_fld_nr_comps(x); m++) {
    char s[10];
    sprintf(s, "m%d", m);
    mrc_fld_set_comp_name(x, m, s);
  }
  mrc_fld_write(x, io);
  mrc_io_close(io);
  mrc_io_destroy(io);
}

// ----------------------------------------------------------------------
// mrc_fld_make_view
//
// FIXME having two different versions of making views in the case of
// domain / no-domain is rather confusing, and this one's not very flexible,
// either...

struct mrc_fld *
mrc_fld_make_view(struct mrc_fld *fld, int mb, int me)
{
  assert(mrc_fld_is_setup(fld));
  assert(fld->_domain);
  assert(fld->_dim == -1);

  struct mrc_fld *fld_new = mrc_fld_create(mrc_fld_comm(fld));
  mrc_fld_set_type(fld_new, mrc_fld_type(fld));
  mrc_fld_set_param_obj(fld_new, "domain", fld->_domain);
  mrc_fld_set_param_int(fld_new, "nr_spatial_dims", fld->_nr_spatial_dims);
  mrc_fld_set_param_int(fld_new, "nr_comps", me - mb);
  mrc_fld_set_param_int(fld_new, "nr_ghosts", fld->_nr_ghosts);
  mrc_fld_set_param_int(fld_new, "dim", fld->_dim);
  mrc_fld_set_param_bool(fld_new, "aos", fld->_aos);
  mrc_fld_set_param_bool(fld_new, "c_order", fld->_c_order);

  fld_new->_view_base = fld;
  int offs[MRC_FLD_MAXDIMS] = {};
  offs[3] = mb;
  fld_new->_view_offs = offs;
  mrc_fld_setup(fld_new);
  fld_new->_view_offs = NULL; // so we don't keep a dangling pointer around

  return fld_new;
}

// ----------------------------------------------------------------------
// mrc_fld_get_as
//
// convert fld_base to mrc_fld of type "type"

struct mrc_fld *
mrc_fld_get_as(struct mrc_fld *fld_base, const char *type)
{
  const char *type_base = mrc_fld_type(fld_base);
  // if we're already the subtype, nothing to be done
  if (strcmp(type_base, type) == 0)
    return fld_base;

#if 0
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    mprintf("XXXXXXX mrc_fld_get_as %s -> %s\n", type_base, type);
    void* callstack[128];
    int frames = backtrace(callstack, 128);
    char** strs = backtrace_symbols(callstack, frames);
    for (int i = 0; i < frames; i++) {
      mprintf("%s\n", strs[i]);
    }
    free(strs);
  }
#endif

  static int pr;
  if (!pr) {
    pr = prof_register("mrc_fld_get_as", 1., 0, 0);
  }
  prof_start(pr);

  struct mrc_fld *fld = mrc_fld_create(mrc_fld_comm(fld_base));
  mrc_fld_set_type(fld, type);
  mrc_fld_set_param_bool(fld, "aos", fld_base->_aos);
  mrc_fld_set_param_bool(fld, "c_order", fld_base->_c_order);
  if (fld_base->_domain) {
    // if we're based on a domain, dims/offs/sw will be set by setup()
    mrc_fld_set_param_obj(fld, "domain", fld_base->_domain);
    mrc_fld_set_param_int(fld, "nr_spatial_dims", fld_base->_nr_spatial_dims);
    mrc_fld_set_param_int(fld, "nr_comps", fld_base->_nr_comps);
    mrc_fld_set_param_int(fld, "nr_ghosts", fld_base->_nr_ghosts);
  } else {
    mrc_fld_set_param_int_array(fld, "dims", fld_base->_dims.nr_vals, fld_base->_dims.vals);
    mrc_fld_set_param_int_array(fld, "offs", fld_base->_offs.nr_vals, fld_base->_offs.vals);
    mrc_fld_set_param_int_array(fld, "sw", fld_base->_sw.nr_vals, fld_base->_sw.vals);
  }
  for (int m = 0; m < fld_base->_nr_comps; m++) {
    mrc_fld_set_comp_name(fld, m, mrc_fld_comp_name(fld_base, m));
  }
  mrc_fld_setup(fld);

  char s[strlen(type) + 12]; sprintf(s, "copy_to_%s", type);
  mrc_fld_copy_to_func_t copy_to = (mrc_fld_copy_to_func_t)
    mrc_fld_get_method(fld_base, s);
  if (copy_to) {
    copy_to(fld_base, fld);
  } else {
    sprintf(s, "copy_from_%s", type_base);
    mrc_fld_copy_from_func_t copy_from = (mrc_fld_copy_from_func_t)
      mrc_fld_get_method(fld, s);
    if (copy_from) {
      copy_from(fld, fld_base);
    } else {
      fprintf(stderr, "ERROR: no 'copy_to_%s' in mrc_fld '%s' and "
	      "no 'copy_from_%s' in '%s'!\n",
	      type, mrc_fld_type(fld_base), type_base, mrc_fld_type(fld));
      assert(0);
    }
  }

  prof_stop(pr);
  return fld;
}

// ----------------------------------------------------------------------
// mrc_fld_put_as
//
// after being done with the fields gotten from get_as(), need to put them
// back using this routine, which will copy the contents from fld back
// to fld_base

void
mrc_fld_put_as(struct mrc_fld *fld, struct mrc_fld *fld_base)
{
  const char *type_base = mrc_fld_type(fld_base);
  const char *type = mrc_fld_type(fld);
  // If we're already the subtype, nothing to be done
  if (strcmp(type_base, type) == 0)
    return;

#if 0
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0 && strcmp(type, "float") == 0) {
    mprintf("XXXXXXX mrc_fld_put_as %s <- %s\n", type_base, type);
    void* callstack[128];
    int frames = backtrace(callstack, 128);
    char** strs = backtrace_symbols(callstack, frames);
    for (int i = 0; i < frames; i++) {
      mprintf("%s\n", strs[i]);
    }
    free(strs);
  }
#endif

  static int pr;
  if (!pr) {
    pr = prof_register("mrc_fld_put_as", 1., 0, 0);
  }
  prof_start(pr);

  // FIXME, this is openggcm specific and shouldn't be handled here
  if (strcmp(type, "float") == 0 && strcmp(type_base, "fortran") == 0) {
    // special case: float from fortran, no need to copy back
    mrc_fld_destroy(fld);
    prof_stop(pr);
    return;
  }
  if (strcmp(type, "fortran") == 0) {
    // special case: convert something else to Fortran
    assert(0); // can't happen
  }

  char s[strlen(type) + 12]; sprintf(s, "copy_from_%s", type);
  mrc_fld_copy_from_func_t copy_from = (mrc_fld_copy_from_func_t)
    mrc_fld_get_method(fld_base, s);
  if (copy_from) {
    copy_from(fld_base, fld);
  } else {
    sprintf(s, "copy_to_%s", type_base);
    mrc_fld_copy_to_func_t copy_to = (mrc_fld_copy_to_func_t)
      mrc_fld_get_method(fld, s);
    if (copy_to) {
      copy_to(fld, fld_base);
    } else {
      fprintf(stderr, "ERROR: no 'copy_from_%s' in mrc_fld '%s' and "
	      "no 'copy_to_%s' in '%s'!\n",
	      type, mrc_fld_type(fld_base), type_base, mrc_fld_type(fld));
      assert(0);
    }
  }

  mrc_fld_destroy(fld);

  prof_stop(pr);
}

#ifdef HAVE_PETSC
// Need to surface petsc_vecs somehow. Interface here is less than ideal, 
// and it's not really clear to me how we should handle this.
#include <petscvec.h>
// ----------------------------------------------------------------------
// mrc_fld_set_petsc_vec
//

void
mrc_fld_set_petsc_vec(struct mrc_fld *fld, Vec petsc_vec)
{
  if (mrc_fld_ops(fld)->vec_type) {
    dispatch_vec_type(fld);
  } else {
    fprintf(stderr, "Cannot set petsc vec for fld type %s (This field doesn't have a vec_type)\n", 
	    mrc_fld_type(fld));
    assert(0);
  }
  if ( strcmp(mrc_fld_ops(fld)->vec_type, "petsc")!=0 ){
    fprintf(stderr, "Cannot set petsc vec for fld type %s", mrc_fld_type(fld));
    fprintf(stderr, " (petsc element size: %zu B, mrc_fld element size: %d B)\n",
	    sizeof(PetscScalar), fld->_nd->size_of_type);
    assert(0);
  }
  void (*vec_set_petsc)( struct mrc_vec *, Vec); 
  vec_set_petsc = (void (*)( struct mrc_vec *, Vec)) mrc_vec_get_method(fld->_vec, "set_petsc_vec");
  assert(vec_set_petsc);
  vec_set_petsc(fld->_vec, petsc_vec);
}

// ----------------------------------------------------------------------
// mrc_fld_get_petsc_vec
// We're not really supposed to surface vectors at all, but petsc
// vecs are sort of a special case, since they need to be passed into petsc
// functions.
Vec 
mrc_fld_get_petsc_vec(struct mrc_fld *fld)
{
  assert(strcmp(mrc_fld_ops(fld)->vec_type, "petsc")==0);
  Vec rv;
  void (*vec_get_petsc)( struct mrc_vec *, Vec *); 
  vec_get_petsc = (void (*)(struct mrc_vec *, Vec *)) mrc_vec_get_method(fld->_vec, "get_petsc_vec");
  vec_get_petsc(fld->_vec, &rv);
  return rv;
}

// ----------------------------------------------------------------------
// mrc_fld_put_petsc_vec
// If you get a petsc vec, you damn well better put it back
void
mrc_fld_put_petsc_vec(struct mrc_fld *fld, Vec *invec)
{
  assert(strcmp(mrc_fld_ops(fld)->vec_type, "petsc")==0);
  void (*vec_put_petsc)(struct mrc_vec * , Vec *);
  vec_put_petsc = (void (*)(struct mrc_vec * , Vec *)) mrc_vec_get_method(fld->_vec, "put_petsc_vec");
  vec_put_petsc(fld->_vec, invec);
}

#endif

void
mrc_fld_ddc_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *buf,
			void *_fld)
{
  struct mrc_fld *fld = _fld;
  assert(mrc_fld_ops(fld)->ddc_copy_to_buf);
  mrc_fld_ops(fld)->ddc_copy_to_buf(fld, mb, me, p, ilo, ihi, buf);
}

void
mrc_fld_ddc_copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *buf,
			  void *_fld)
{
  struct mrc_fld *fld = _fld;
  assert(mrc_fld_ops(fld)->ddc_copy_from_buf);
  mrc_fld_ops(fld)->ddc_copy_from_buf(fld, mb, me, p, ilo, ihi, buf);
}

void
mrc_fld_ddc_add_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *buf,
			 void *_fld)
{
  struct mrc_fld *fld = _fld;
  assert(mrc_fld_ops(fld)->ddc_add_from_buf);
  mrc_fld_ops(fld)->ddc_add_from_buf(fld, mb, me, p, ilo, ihi, buf);
}

// ======================================================================
// mrc_fld subclasses

// ----------------------------------------------------------------------
// mrc_fld_float_copy_from_double

static void
mrc_fld_float_copy_from_double(struct mrc_fld *to, struct mrc_fld *from)
{
  assert(mrc_fld_same_shape(to, from));
  assert(mrc_fld_data_type(to) == MRC_NT_FLOAT);
  assert(mrc_fld_data_type(from) == MRC_NT_DOUBLE);

  if (to->_dims.nr_vals == 5) {
    for (int i4 = to->_ghost_offs[4]; i4 < to->_ghost_offs[4] + to->_ghost_dims[4]; i4++) {
      for (int i3 = to->_ghost_offs[3]; i3 < to->_ghost_offs[3] + to->_ghost_dims[3]; i3++) {
	for (int i2 = to->_ghost_offs[2]; i2 < to->_ghost_offs[2] + to->_ghost_dims[2]; i2++) {
	  for (int i1 = to->_ghost_offs[1]; i1 < to->_ghost_offs[1] + to->_ghost_dims[1]; i1++) {
	    for (int i0 = to->_ghost_offs[0]; i0 < to->_ghost_offs[0] + to->_ghost_dims[0]; i0++) {
	      MRC_S5(to, i0,i1,i2,i3,i4) = MRC_D5(from, i0,i1,i2,i3,i4);
	    }
	  }
	}
      }
    }
  } else if (to->_dims.nr_vals == 4) {
    for (int i3 = to->_ghost_offs[3]; i3 < to->_ghost_offs[3] + to->_ghost_dims[3]; i3++) {
      for (int i2 = to->_ghost_offs[2]; i2 < to->_ghost_offs[2] + to->_ghost_dims[2]; i2++) {
	for (int i1 = to->_ghost_offs[1]; i1 < to->_ghost_offs[1] + to->_ghost_dims[1]; i1++) {
	  for (int i0 = to->_ghost_offs[0]; i0 < to->_ghost_offs[0] + to->_ghost_dims[0]; i0++) {
	    MRC_S4(to, i0,i1,i2,i3) = MRC_D4(from, i0,i1,i2,i3);
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mrc_fld_float_copy_to_double

static void
mrc_fld_float_copy_to_double(struct mrc_fld *from, struct mrc_fld *to)
{
  assert(mrc_fld_same_shape(to, from));
  assert(mrc_fld_data_type(to) == MRC_NT_DOUBLE);
  assert(mrc_fld_data_type(from) == MRC_NT_FLOAT);

  if (to->_dims.nr_vals == 5) {
    for (int i4 = to->_ghost_offs[4]; i4 < to->_ghost_offs[4] + to->_ghost_dims[4]; i4++) {
      for (int i3 = to->_ghost_offs[3]; i3 < to->_ghost_offs[3] + to->_ghost_dims[3]; i3++) {
	for (int i2 = to->_ghost_offs[2]; i2 < to->_ghost_offs[2] + to->_ghost_dims[2]; i2++) {
	  for (int i1 = to->_ghost_offs[1]; i1 < to->_ghost_offs[1] + to->_ghost_dims[1]; i1++) {
	    for (int i0 = to->_ghost_offs[0]; i0 < to->_ghost_offs[0] + to->_ghost_dims[0]; i0++) {
	      MRC_D5(to, i0,i1,i2,i3,i4) = MRC_S5(from, i0,i1,i2,i3,i4);
	    }
	  }
	}
      }
    }
  } else if (to->_dims.nr_vals == 4) {
    for (int i3 = to->_ghost_offs[3]; i3 < to->_ghost_offs[3] + to->_ghost_dims[3]; i3++) {
      for (int i2 = to->_ghost_offs[2]; i2 < to->_ghost_offs[2] + to->_ghost_dims[2]; i2++) {
	for (int i1 = to->_ghost_offs[1]; i1 < to->_ghost_offs[1] + to->_ghost_dims[1]; i1++) {
	  for (int i0 = to->_ghost_offs[0]; i0 < to->_ghost_offs[0] + to->_ghost_dims[0]; i0++) {
	    MRC_D4(to, i0,i1,i2,i3) = MRC_S4(from, i0,i1,i2,i3);
	  }
	}
      }
    }
  } else {
    assert(0);
  }
}


// ----------------------------------------------------------------------
// mrc_fld "float" methods

static struct mrc_obj_method mrc_fld_float_methods[] = {
  MRC_OBJ_METHOD("copy_to_double",   mrc_fld_float_copy_to_double),
  MRC_OBJ_METHOD("copy_from_double", mrc_fld_float_copy_from_double),
  {}
};

// ----------------------------------------------------------------------
// mrc_fld_*_methods

static struct mrc_obj_method mrc_fld_double_methods[] = {
  {}
};

static struct mrc_obj_method mrc_fld_int_methods[] = {
  {}
};

// no ddc support for mrc_fld "int" (though it could be added)

void
mrc_fld_int_ddc_copy_from_buf(struct mrc_fld *fld, int mb, int me, int p,
			      int ilo[3], int ihi[3], void *buf)
{
  assert(0);
}

void
mrc_fld_int_ddc_copy_to_buf(struct mrc_fld *fld, int mb, int me, int p,
			    int ilo[3], int ihi[3], void *buf)
{
  assert(0);
}


// ----------------------------------------------------------------------
// create float, double, int subclasses

#define MAKE_MRC_FLD_TYPE(NAME, type, TYPE)				\
									\
  void mrc_fld_##NAME##_ddc_copy_from_buf(struct mrc_fld *, int, int,	\
					  int, int[3], int[3], void *); \
  void mrc_fld_##NAME##_ddc_copy_to_buf(struct mrc_fld *, int, int,	\
					int, int[3], int[3], void *);	\
  									\
  static struct mrc_fld_ops mrc_fld_##NAME##_ops = {			\
    .name                  = #NAME,					\
    .methods               = mrc_fld_##NAME##_methods,			\
    .ddc_copy_from_buf	   = mrc_fld_##NAME##_ddc_copy_from_buf,	\
    .ddc_copy_to_buf	   = mrc_fld_##NAME##_ddc_copy_to_buf,		\
    .vec_type              = #type,					\
  };									\

MAKE_MRC_FLD_TYPE(float, float, FLOAT)
MAKE_MRC_FLD_TYPE(double, double, DOUBLE)
MAKE_MRC_FLD_TYPE(int, int, INT)

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
  { "offs"            , VAR(_offs)           , PARAM_INT_ARRAY(0, 0) },
  { "dims"            , VAR(_dims)           , PARAM_INT_ARRAY(0, 0) },
  { "sw"              , VAR(_sw)             , PARAM_INT_ARRAY(0, 0) },

  { "domain"          , VAR(_domain)         , PARAM_OBJ(mrc_domain) },
  { "nr_spatial_dims" , VAR(_nr_spatial_dims), PARAM_INT(-1)         },
  { "dim"             , VAR(_dim)            , PARAM_INT(-1)         },
  { "nr_comps"        , VAR(_nr_comps)       , PARAM_INT(1)          },
  { "nr_ghosts"       , VAR(_nr_ghosts)      , PARAM_INT(0)          },

  // FIXME, if _nd is listed as member obj, we can't prevent it being written,
  // but usually, we want to write the data ourselves, and don't write nd
  //  { "nd"              , VAR(_nd)             , MRC_VAR_OBJ(mrc_ndarray) },
  { "aos"             , VAR(_aos)            , PARAM_BOOL(false)     },
  { "c_order"         , VAR(_c_order)        , PARAM_BOOL(false)     },
  {},
};
#undef VAR

static struct mrc_obj_method mrc_fld_methods[] = {
  MRC_OBJ_METHOD("duplicate", mrc_fld_duplicate),
  MRC_OBJ_METHOD("copy"     , mrc_fld_copy),
  MRC_OBJ_METHOD("axpy"     , mrc_fld_axpy),
  MRC_OBJ_METHOD("waxpy"    , mrc_fld_waxpy),
  MRC_OBJ_METHOD("norm"     , mrc_fld_norm),
  MRC_OBJ_METHOD("set"      , mrc_fld_set),
#ifdef HAVE_PETSC
  MRC_OBJ_METHOD("set_petsc_vec", mrc_fld_set_petsc_vec),
  MRC_OBJ_METHOD("get_petsc_vec", mrc_fld_get_petsc_vec),
  MRC_OBJ_METHOD("put_petsc_vec", mrc_fld_put_petsc_vec),
#endif
  {}
};

// ----------------------------------------------------------------------
// mrc_fld class description

struct mrc_class_mrc_fld mrc_class_mrc_fld = {
  .name         = "mrc_fld",
  .size         = sizeof(struct mrc_fld),
  .param_descr  = mrc_fld_descr,
  .methods      = mrc_fld_methods,
  .init         = mrc_fld_init,
  .create       = _mrc_fld_create,
  .destroy      = _mrc_fld_destroy,
  .setup        = _mrc_fld_setup,
  .read         = _mrc_fld_read,
  .write        = _mrc_fld_write,
};



