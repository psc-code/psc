
#include "mrc_io_private.h"
#include "mrc_params.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#ifdef HAVE_HDF5_H
#include <hdf5.h>
#endif

// ======================================================================
// mrc_io

static inline struct mrc_io_ops *
mrc_io_ops(struct mrc_io *io)
{
  return (struct mrc_io_ops *) io->obj.ops;
}

// ----------------------------------------------------------------------

struct mrc_obj_entry {
  struct mrc_obj *obj;
  char *path;
  list_t entry;
};

struct mrc_obj *
mrc_io_find_obj(struct mrc_io *io, const char *path)
{
  struct mrc_obj_entry *p;
  __list_for_each_entry(p, &io->obj_list, entry, struct mrc_obj_entry) {
    if (strcmp(p->path, path) == 0) {
      return p->obj;
    }
  }
  return NULL;
}

const char *
__mrc_io_obj_path(struct mrc_io *io, struct mrc_obj *obj)
{
  struct mrc_obj_entry *p;
  __list_for_each_entry(p, &io->obj_list, entry, struct mrc_obj_entry) {
    if (p->obj == obj) {
      return p->path;
    }
  }
  return NULL;
}

int
mrc_io_add_obj(struct mrc_io *io, struct mrc_obj *obj, const char *_path)
{
  char *path = strdup(_path);
  struct mrc_obj *obj2 = mrc_io_find_obj(io, path);
  if (obj2) { // exists
    if (obj != obj2) {
      mprintf("!!! obj  %p '%s' (%s)\n", obj, path, obj->cls->name);
      mprintf("!!! obj2 %p '%s' (%s)\n", obj2, path, obj2->cls->name);
      assert(0);
    }
    free(path);
    return 1;
  }
  struct mrc_obj_entry *p = malloc(sizeof(*p));
  p->obj = mrc_obj_get(obj);
  p->path = path;
  list_add_tail(&p->entry, &io->obj_list);
  return 0;
}

// ----------------------------------------------------------------------
// mrc_io_setup

static void
_mrc_io_setup(struct mrc_io *io)
{
  MPI_Comm_rank(io->obj.comm, &io->rank);
  MPI_Comm_size(io->obj.comm, &io->size);
}

// ----------------------------------------------------------------------
// mrc_io_open

void
mrc_io_open(struct mrc_io *io, const char *mode, int step, float time)
{
  assert(mrc_io_is_setup(io));

  struct mrc_io_ops *ops = mrc_io_ops(io);
  io->step = step;
  io->time = time;
  INIT_LIST_HEAD(&io->obj_list);
  assert(ops->open);
  ops->open(io, mode);
}

// ----------------------------------------------------------------------
// mrc_io_close

void
mrc_io_close(struct mrc_io *io)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  while (!list_empty(&io->obj_list)) {
    struct mrc_obj_entry *p = list_entry(io->obj_list.next, struct mrc_obj_entry,
					 entry);
    list_del(&p->entry);
    mrc_obj_put(p->obj);
    free(p->path);
    free(p);
  }
  assert(ops->close);
  ops->close(io);
  io->step = -1;
  io->time = -1.;
}

// ----------------------------------------------------------------------
// mrc_io_read_f3

void
mrc_io_read_f3(struct mrc_io *io, const char *path, struct mrc_fld *fld)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->read_f3) {
    ops->read_f3(io, path, fld);
  } else {
    assert(fld->_domain);
    struct mrc_fld *m3 = mrc_domain_m3_create(fld->_domain);
    mrc_fld_set_param_int(m3, "nr_ghosts", fld->_sw.vals[0]);
    mrc_fld_set_param_int(m3, "nr_comps", mrc_fld_nr_comps(fld));
    mrc_fld_setup(m3);
    mrc_io_read_m3(io, path, m3);

    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, 0);
    for (int m = 0; m < mrc_fld_nr_comps(m3); m++) {
      mrc_m3_foreach_bnd(m3p, ix,iy,iz) {
	MRC_F3(fld, m, ix,iy,iz) = MRC_M3(m3p, m, ix,iy,iz);
      } mrc_m3_foreach_end;
    }
    mrc_fld_patch_put(m3);
    mrc_fld_destroy(m3);
  }
}

// ----------------------------------------------------------------------
// mrc_io_read_fld

void
mrc_io_read_fld(struct mrc_io *io, const char *path, struct mrc_fld *fld)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->read_fld) {
    ops->read_fld(io, path, fld);
  } else if (fld->_dims.nr_vals == 3) {
    if (strcmp(mrc_fld_type(fld), "float") == 0) {
      mrc_io_read_m1(io, path, fld);
    }
  } else if (fld->_dims.nr_vals == 5) {
    mrc_io_read_m3(io, path, fld);
  } else {
    assert(0);
  }
  // FIXME: To save/recover comp names, just write them as attributes
  for ( int m=0; m < mrc_fld_nr_comps(fld); m++) {
    char comp_label[100];
    sprintf(comp_label, "comp_name_%d", m);
    char *comp_name = NULL;
    // FIXME: Is the string memory returned from this leaked?
    mrc_io_read_attr_string(io, path, (const char *) comp_label, &comp_name);
    if (comp_name) {
      mrc_fld_set_comp_name(fld, m, (const char *) comp_name);
      free(comp_name); // I think this free is needed, otherwise we leak the returned memory
    }
  }
}

// ----------------------------------------------------------------------
// mrc_io_write_fld

void
mrc_io_write_fld(struct mrc_io *io, const char *path, struct mrc_fld *fld)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_fld) {
    ops->write_fld(io, path, fld);
  } else if (fld->_dims.nr_vals == 3) {
    if (strcmp(mrc_fld_type(fld), "float") == 0) {
      mrc_io_write_m1(io, path, fld);
    }
  } else if (fld->_dims.nr_vals == 5) {
    mrc_io_write_m3(io, path, fld);
  } else if (fld->_dims.nr_vals == 1) {
    MHERE;
  } else {
    assert(0);
  }
  // FIXME: To save/recover comp names, just write them as attributes
  for ( int m=0; m < mrc_fld_nr_comps(fld); m++) {
    char comp_label[100];
    sprintf(comp_label, "comp_name_%d", m);
    mrc_io_write_attr_string(io, path, (const char *) comp_label, mrc_fld_comp_name(fld, m));
  }
}

// ----------------------------------------------------------------------
// mrc_io_write_m1

void
mrc_io_write_m1(struct mrc_io *io, const char *path, struct mrc_fld *fld)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_m1) {
    ops->write_m1(io, path, fld);
  } else {
    MHERE; // FIXME
  }
}

// ----------------------------------------------------------------------
// mrc_io_read_m1

void
mrc_io_read_m1(struct mrc_io *io, const char *path, struct mrc_fld *fld)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->read_m1) {
    ops->read_m1(io, path, fld);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mrc_io_write_m3

void
mrc_io_write_m3(struct mrc_io *io, const char *path, struct mrc_fld *fld)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_m3) {
    ops->write_m3(io, path, fld);
  } else if (ops->write_field) {
    assert(mrc_fld_nr_patches(fld) == 1);
    int nr_comps = mrc_fld_nr_comps(fld);
    for (int m = 0; m < nr_comps; m++) {
      assert(mrc_fld_comp_name(fld, m));
      ops->write_field(io, path, 1., fld, m);
    }
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mrc_io_read_m3

void
mrc_io_read_m3(struct mrc_io *io, const char *path, struct mrc_fld *fld)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->read_m3) {
    ops->read_m3(io, path, fld);
  } else if (ops->read_f3) {
    assert(mrc_fld_nr_patches(fld) == 1);
    mrc_io_read_f3(io, path, fld);
  } else {
    assert(0);
  }
}

// ----------------------------------------------------------------------
// mrc_io_write_field2d

void
mrc_io_write_field2d(struct mrc_io *io, float scale, struct mrc_fld *fld,
		     int outtype, float sheet)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(mrc_fld_comp_name(fld, 0));
  assert(ops->write_field2d);
  ops->write_field2d(io, scale, fld, outtype, sheet);
}

// ----------------------------------------------------------------------
// mrc_io_write_field_slice

void
mrc_io_write_field_slice(struct mrc_io *io, float scale, struct mrc_fld *fld,
			 int outtype, float sheet)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);

  //dig out a sheet of constant x/y/z

  assert(outtype >= DIAG_TYPE_2D_X && outtype <= DIAG_TYPE_2D_Z);
  int dim = outtype - DIAG_TYPE_2D_X;

  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(fld->_domain, &nr_patches);
  assert(nr_patches == 1);
  int *dims = patches[0].ldims;

  struct mrc_fld *f2 = mrc_fld_create(MPI_COMM_SELF);
  //check for existence on local proc.
  //0,1, nnx-2, nnx-1 are the ghostpoints
  struct mrc_crds *crds = mrc_domain_get_crds(fld->_domain);
  if (MRC_CRD(crds, dim, 1) < sheet && sheet <= MRC_CRD(crds, dim, dims[dim]+1)) { 
    int ii;
    for (ii = 1; ii < dims[dim]; ii++) {
      if (sheet <= MRC_CRD(crds, dim, ii+1))
	break;
    }

    float s2 = (sheet - MRC_CRD(crds, dim, ii)) / 
      (MRC_CRD(crds, dim, ii+1) - MRC_CRD(crds, dim, ii));
    float s1 = 1. - s2;
    s1 *= scale;
    s2 *= scale;

    struct mrc_fld *f = mrc_fld_get_as(fld, "float");
    switch (dim) {
    case 0: {
      mrc_fld_set_param_int_array(f2, "dims", 5, (int[5]) { dims[1], dims[2], 1, 1, 1 });
      mrc_fld_setup(f2);
      struct mrc_fld *_f2 = mrc_fld_get_as(f2, "float");
      for(int iz = 0; iz < dims[2]; iz++) {
	for(int iy = 0; iy < dims[1]; iy++) {
	  MRC_F2(_f2,0, iy,iz) = (s1 * MRC_F3(f,0, ii-2,iy,iz) +
				  s2 * MRC_F3(f,0, ii-1,iy,iz));
	}
      }
      mrc_fld_put_as(_f2, f2);
    }
      break;
    case 1: {
      mrc_fld_set_param_int_array(f2, "dims", 5, (int[5]) { dims[0], dims[2], 1, 1, 1 });
      mrc_fld_setup(f2);
      struct mrc_fld *_f2 = mrc_fld_get_as(f2, "float");
      for(int iz = 0; iz < dims[2]; iz++) {
	for(int ix = 0; ix < dims[0]; ix++) {
	  MRC_F2(_f2,0, ix,iz) = (s1 * MRC_F3(f,0, ix,ii-2,iz) +
				  s2 * MRC_F3(f,0, ix,ii-1,iz));
	}
      }
      mrc_fld_put_as(_f2, f2);
      break;
    }
    case 2: {
      mrc_fld_set_param_int_array(f2, "dims", 5, (int[5]) { dims[0], dims[1], 1, 1, 1 });
      mrc_fld_setup(f2);
      struct mrc_fld *_f2 = mrc_fld_get_as(f2, "float");
      for(int iy = 0; iy < dims[1]; iy++) {
	for(int ix = 0; ix < dims[0]; ix++) {
	  MRC_F2(_f2,0, ix,iy) = (s1 * MRC_F3(f,0, ix,iy,ii-2) +
				  s2 * MRC_F3(f,0, ix,iy,ii-1));
	}
      }
      mrc_fld_put_as(_f2, f2);
      break;
    }
    }
    mrc_fld_put_as(f, fld);
  } else {
    // not on local proc
    mrc_fld_set_param_int_array(f2, "dims", 5, (int[5]) { 0, 0, 0, 1, 1 });
    mrc_fld_setup(f2);
  }
  f2->_domain = fld->_domain; // FIXME, 2D field isn't directly based on domain, so shouldn't link to it (?)
  mrc_fld_set_comp_name(f2, 0, mrc_fld_comp_name(fld, 0));
  ops->write_field2d(io, 1., f2, outtype, sheet);
  mrc_fld_destroy(f2);
}

// ----------------------------------------------------------------------
// mrc_io_write_ndarray

void
mrc_io_write_ndarray(struct mrc_io *io, const char *path, struct mrc_ndarray *nd)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops && ops->write_ndarray);

  ops->write_ndarray(io, path, nd);
}

// ----------------------------------------------------------------------
// mrc_io_read_ndarray

void
mrc_io_read_ndarray(struct mrc_io *io, const char *path, struct mrc_ndarray *nd)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops && ops->read_ndarray);

  ops->read_ndarray(io, path, nd);
}

// ----------------------------------------------------------------------
// mrc_io_read_attr

void
mrc_io_read_attr(struct mrc_io *io, const char *path, int type, const char *name,
		 union param_u *pv)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->read_attr);
  ops->read_attr(io, path, type, name, pv);
}

void
mrc_io_read_attr_int(struct mrc_io *io, const char *path, const char *name,
		     int *val)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->read_attr);
  union param_u u;
  ops->read_attr(io, path, PT_INT, name, &u);
  *val = u.u_int;
}

void
mrc_io_read_attr_int3(struct mrc_io *io, const char *path, const char *name,
		      int (*val)[3])
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->read_attr);
  union param_u u;
  ops->read_attr(io, path, PT_INT3, name, &u);
  for (int d = 0; d < 3; d++) {
    (*val)[d] = u.u_int3[d];
  }
}

void
mrc_io_read_attr_float(struct mrc_io *io, const char *path, const char *name,
		     float *val)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->read_attr);
  union param_u u;
  ops->read_attr(io, path, PT_FLOAT, name, &u);
  *val = u.u_float;
}

void
mrc_io_read_attr_float3(struct mrc_io *io, const char *path, const char *name,
		      float (*val)[3])
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->read_attr);
  union param_u u;
  ops->read_attr(io, path, PT_FLOAT3, name, &u);
  for (int d = 0; d < 3; d++) {
    (*val)[d] = u.u_float3[d];
  }
}

void
mrc_io_read_attr_double(struct mrc_io *io, const char *path, const char *name,
			double *val)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->read_attr);
  union param_u u;
  ops->read_attr(io, path, PT_DOUBLE, name, &u);
  *val = u.u_double;
}

void
mrc_io_read_attr_double3(struct mrc_io *io, const char *path, const char *name,
		      double (*val)[3])
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->read_attr);
  union param_u u;
  ops->read_attr(io, path, PT_DOUBLE3, name, &u);
  for (int d = 0; d < 3; d++) {
    (*val)[d] = u.u_double3[d];
  }
}


void
mrc_io_read_attr_string(struct mrc_io *io, const char *path, const char *name,
			char **val)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->read_attr);
  union param_u u;
  ops->read_attr(io, path, PT_STRING, name, &u);
  *val = (char *) u.u_string;
}

void
mrc_io_read_attr_bool(struct mrc_io *io, const char *path, const char *name,
		      bool *val)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->read_attr);
  union param_u u;
  ops->read_attr(io, path, PT_BOOL, name, &u);
  *val = u.u_bool;
}


// ----------------------------------------------------------------------
// mrc_io_write_attr

void
mrc_io_write_attr(struct mrc_io *io, const char *path, int type, const char *name,
		  union param_u *pv)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_attr) {
    ops->write_attr(io, path, type, name, pv);
  }
}

void
mrc_io_write_attr_int(struct mrc_io *io, const char *path, const char *name,
		      int val)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_attr) {
    union param_u u = { .u_int = val };
    ops->write_attr(io, path, PT_INT, name, &u);
  }
}

void
mrc_io_write_attr_int3(struct mrc_io *io, const char *path, const char *name,
		       int val[3])
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_attr) {
    union param_u u = { .u_int3 = { val[0], val[1], val[2] } };
    ops->write_attr(io, path, PT_INT3, name, &u);
  }
}



void
mrc_io_write_attr_float(struct mrc_io *io, const char *path, const char *name,
		      float val)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_attr) {
    union param_u u = { .u_float = val };
    ops->write_attr(io, path, PT_FLOAT, name, &u);
  }
}

void
mrc_io_write_attr_float3(struct mrc_io *io, const char *path, const char *name,
		       float val[3])
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_attr) {
    union param_u u = { .u_float3 = { val[0], val[1], val[2] } };
    ops->write_attr(io, path, PT_FLOAT3, name, &u);
  }
}



void
mrc_io_write_attr_string(struct mrc_io *io, const char *path, const char *name,
			 const char *val)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_attr) {
    union param_u u = { .u_string = val };
    ops->write_attr(io, path, PT_STRING, name, &u);
  }
}

void
mrc_io_write_attr_double(struct mrc_io *io, const char *path, const char *name,
			 double val)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_attr) {
    union param_u u = { .u_double = val };
    ops->write_attr(io, path, PT_DOUBLE, name, &u);
  }
}

void
mrc_io_write_attr_double3(struct mrc_io *io, const char *path, const char *name,
		       double val[3])
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_attr) {
    union param_u u = { .u_double3 = { val[0], val[1], val[2] } };
    ops->write_attr(io, path, PT_DOUBLE3, name, &u);
  }
}

void
mrc_io_write_attr_bool(struct mrc_io *io, const char *path, const char *name,
		       bool val)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_attr) {
    union param_u u = { .u_bool = val };
    ops->write_attr(io, path, PT_BOOL, name, &u);
  }
}

// ----------------------------------------------------------------------
// mrc_io_get_h5_file

void
mrc_io_get_h5_file(struct mrc_io *io, long *h5_file)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->get_h5_file);
  ops->get_h5_file(io, h5_file);
}

// ----------------------------------------------------------------------
// __mrc_io_read_path

struct mrc_obj *
__mrc_io_read_path(struct mrc_io *io, const char *path, const char *name,
		   struct mrc_class *class)
{
  char *s;
  mrc_io_read_attr_string(io, path, name, &s);
  struct mrc_obj *obj = mrc_obj_read(io, s, class);
  free(s);
  return obj;
}

// ----------------------------------------------------------------------
// __mrc_io_read_ref

struct mrc_obj *
__mrc_io_read_ref(struct mrc_io *io, struct mrc_obj *obj_parent, const char *name,
		  struct mrc_class *class)
{
  char *s;
  mrc_io_read_attr_string(io, mrc_io_obj_path(io, obj_parent), name, &s);
  if (!s) {
    return NULL;
  }
  struct mrc_obj *obj = mrc_obj_read(io, s, class);
  free(s);
  return obj;
}

// ----------------------------------------------------------------------
// __mrc_io_read_ref_comm

struct mrc_obj *
__mrc_io_read_ref_comm(struct mrc_io *io, struct mrc_obj *obj_parent, const char *name,
		       struct mrc_class *class, MPI_Comm comm)
{
  char *s;
  mrc_io_read_attr_string(io, mrc_io_obj_path(io, obj_parent), name, &s);
  if (!s) {
    return NULL;
  }
  struct mrc_obj *obj = mrc_obj_read_comm(io, s, class, comm);
  free(s);
  return obj;
}

// ----------------------------------------------------------------------
// __mrc_io_write_path

void
__mrc_io_write_path(struct mrc_io *io, const char *path, const char *name,
		    struct mrc_obj *obj)
{
  mrc_obj_write(obj, io);
  mrc_io_write_attr_string(io, path, name, mrc_io_obj_path(io, obj));
}

// ----------------------------------------------------------------------
// __mrc_io_write_ref

void
__mrc_io_write_ref(struct mrc_io *io, struct mrc_obj *obj_parent, const char *name,
		   struct mrc_obj *obj)
{
  if (obj) {
    mrc_obj_write(obj, io);
    mrc_io_write_attr_string(io, mrc_io_obj_path(io, obj_parent), name,
			     mrc_io_obj_path(io, obj));
  } else {
    mrc_io_write_attr_string(io, mrc_io_obj_path(io, obj_parent), name,
			     "(NULL)");
  }
}

// ======================================================================
// mrc_io_init

static void
mrc_io_init()
{
#ifdef HAVE_HDF5
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf_collective_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf2_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf_serial_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf_to_one_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_hdf5_serial_ops);
#ifdef H5_HAVE_PARALLEL
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf_parallel_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_hdf5_parallel_ops);
#endif
#endif
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_ascii_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_vpic_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_combined_ops);
  // ========================================
  // Deprecated / eternally broken io types (FIXME)
  // ====================
  // mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf2_parallel_ops);

}

// ======================================================================
// mrc_io class

#define VAR(x) (void *)offsetof(struct mrc_io_params, x)
static struct param mrc_io_params_descr[] = {
  { "outdir"          , VAR(outdir)       , PARAM_STRING(".")      },
  { "basename"        , VAR(basename)     , PARAM_STRING("run")    },
  {},
};
#undef VAR

struct mrc_class_mrc_io mrc_class_mrc_io = {
  .name         = "mrc_io",
  .size         = sizeof(struct mrc_io),
  .param_descr  = mrc_io_params_descr,
  .param_offset = offsetof(struct mrc_io, par),
  .init         = mrc_io_init,
  .setup        = _mrc_io_setup,
};

