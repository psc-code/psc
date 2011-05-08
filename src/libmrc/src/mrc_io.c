
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

#define check_is_setup(io) do { assert(io->is_setup); } while (0)

static inline struct mrc_io_ops *
mrc_io_ops(struct mrc_io *io)
{
  return (struct mrc_io_ops *) io->obj.ops;
}

// ----------------------------------------------------------------------

struct mrc_obj_entry {
  struct mrc_obj *obj;
  list_t entry;
};

struct mrc_obj *
mrc_io_find_obj(struct mrc_io *io, const char *name)
{
  struct mrc_obj_entry *p;
  list_for_each_entry(p, &io->obj_list, entry) {
    if (strcmp(mrc_obj_name(p->obj), name) == 0) {
      return p->obj;
    }
  }
  return NULL;
}

int
mrc_io_add_obj(struct mrc_io *io, struct mrc_obj *obj)
{
  struct mrc_obj *obj2 = mrc_io_find_obj(io, mrc_obj_name(obj));
  if (obj2) { // exists
    if (obj->class == obj2->class && obj != obj2) { // FIXME
      mprintf("!!! obj  %p '%s' (%s)\n", obj, mrc_obj_name(obj), obj->class->name);
      mprintf("!!! obj2 %p '%s' (%s)\n", obj2, mrc_obj_name(obj2), obj2->class->name);
      assert(0);
    }
    return 1;
  }
  struct mrc_obj_entry *p = malloc(sizeof(*p));
  p->obj = mrc_obj_get(obj);
  //  mprintf("add obj %p %s\n", obj, mrc_obj_name(obj));
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
  if (mrc_io_ops(io)->setup) {
    mrc_io_ops(io)->setup(io);
  }
  io->is_setup = true;
}

// ----------------------------------------------------------------------
// mrc_io_open

void
mrc_io_open(struct mrc_io *io, const char *mode, int step, float time)
{
  check_is_setup(io);
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
    free(p);
  }
  assert(ops->close);
  ops->close(io);
  io->step = -1;
  io->time = -1.;
}

// ----------------------------------------------------------------------
// mrc_io_read_f1

void
mrc_io_read_f1(struct mrc_io *io, const char *path, struct mrc_f1 *fld)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->read_f1);
  ops->read_f1(io, path, fld);
}

// ----------------------------------------------------------------------
// mrc_io_read_f3

void
mrc_io_read_f3(struct mrc_io *io, const char *path, struct mrc_f3 *fld)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->read_f3);
  ops->read_f3(io, path, fld);
}

// ----------------------------------------------------------------------
// mrc_io_write_f1

void
mrc_io_write_f1(struct mrc_io *io, const char *path, struct mrc_f1 *fld)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_f1) { // FIXME
    ops->write_f1(io, path, fld);
  }
}

// ----------------------------------------------------------------------
// mrc_io_write_f3

void
mrc_io_write_f3(struct mrc_io *io, const char *path,
		struct mrc_f3 *fld, float scale)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_f3) {
    ops->write_f3(io, path, fld, scale);
  } else {
    for (int m = 0; m < fld->nr_comp; m++) {
      assert(fld->name[m]);
      ops->write_field(io, path, scale, fld, m);
    }
  }
}

// ----------------------------------------------------------------------
// mrc_io_write_m3

void
mrc_io_write_m3(struct mrc_io *io, const char *path, struct mrc_m3 *fld)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->write_m3);
  ops->write_m3(io, path, fld);
}

// ----------------------------------------------------------------------
// mrc_io_write_field2d

void
mrc_io_write_field2d(struct mrc_io *io, float scale, struct mrc_f2 *fld,
		     int outtype, float sheet)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (!fld->name[0]) {
    char s[10];
    sprintf(s, "m%d", 0);
    fld->name[0] = strdup(s);
  }
  ops->write_field2d(io, scale, fld, outtype, sheet);
}

// ----------------------------------------------------------------------
// mrc_io_write_field_slice

void
mrc_io_write_field_slice(struct mrc_io *io, float scale, struct mrc_f3 *fld,
			 int outtype, float sheet)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);

  //dig out a sheet of constant x/y/z

  assert(outtype >= DIAG_TYPE_2D_X && outtype <= DIAG_TYPE_2D_Z);
  int dim = outtype - DIAG_TYPE_2D_X;

  int nr_patches;
  struct mrc_patch *patches = mrc_domain_get_patches(fld->domain, &nr_patches);
  assert(nr_patches == 1);
  int *dims = patches[0].ldims;

  //check for existence on local proc.
  //0,1, nnx-2, nnx-1 are the ghostpoints
  struct mrc_crds *crds = mrc_domain_get_crds(fld->domain);
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

    struct mrc_f2 f2;
    switch (dim) {
    case 0:
      mrc_f2_alloc(&f2, NULL, (int [2]) { dims[1], dims[2] }, 1);
      for(int iz = 0; iz < dims[2]; iz++) {
	for(int iy = 0; iy < dims[1]; iy++) {
	  MRC_F2(&f2,0, iy,iz) = (s1 * MRC_F3(fld,0, ii-2,iy,iz) +
				  s2 * MRC_F3(fld,0, ii-1,iy,iz));
	}
      }
      break;
    case 1:
      mrc_f2_alloc(&f2, NULL, (int [2]) { dims[0], dims[2] }, 1);
      for(int iz = 0; iz < dims[2]; iz++) {
	for(int ix = 0; ix < dims[0]; ix++) {
	  MRC_F2(&f2,0, ix,iz) = (s1 * MRC_F3(fld,0, ix,ii-2,iz) +
				  s2 * MRC_F3(fld,0, ix,ii-1,iz));
	}
      }
      break;
    case 2:
      mrc_f2_alloc(&f2, NULL, (int [2]) { dims[0], dims[1] }, 1);
      for(int iy = 0; iy < dims[1]; iy++) {
	for(int ix = 0; ix < dims[0]; ix++) {
	  MRC_F2(&f2,0, ix,iy) = (s1 * MRC_F3(fld,0, ix,iy,ii-2) +
				  s2 * MRC_F3(fld,0, ix,iy,ii-1));
	}
      }
      break;
    }

    f2.domain = fld->domain;
    f2.name[0] = strdup(fld->name[0]);
    ops->write_field2d(io, 1., &f2, outtype, sheet);
    mrc_f2_free(&f2);
  } else {
    struct mrc_f2 f2 = {};
    f2.domain = fld->domain;
    f2.name = calloc(1, sizeof(*f2.name));
    f2.name[0] = strdup(fld->name[0]);
    ops->write_field2d(io, 1., &f2, outtype, sheet);
    mrc_f2_free(&f2);
  }
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
mrc_io_read_attr_string(struct mrc_io *io, const char *path, const char *name,
			char **val)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  assert(ops->read_attr);
  union param_u u;
  ops->read_attr(io, path, PT_STRING, name, &u);
  *val = (char *) u.u_string;
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
mrc_io_write_attr_string(struct mrc_io *io, const char *path, const char *name,
			 const char *val)
{
  struct mrc_io_ops *ops = mrc_io_ops(io);
  if (ops->write_attr) {
    union param_u u = { .u_string = val };
    ops->write_attr(io, path, PT_STRING, name, &u);
  }
}

// ----------------------------------------------------------------------

struct mrc_obj *
__mrc_io_read_obj_ref(struct mrc_io *io, const char *path, const char *name,
		      struct mrc_class *class)
{
  char *s;
  mrc_io_read_attr_string(io, path, name, &s);
  struct mrc_obj *obj = mrc_obj_read(io, s, class);
  free(s);
  return obj;
}

void
mrc_io_write_obj_ref(struct mrc_io *io, const char *path, const char *name,
		     struct mrc_obj *obj)
{
  mrc_io_write_attr_string(io, path, name, mrc_obj_name(obj));
  mrc_obj_write(obj, io);
}

// ======================================================================
// mrc_io_init

static void
mrc_io_init()
{
#ifdef HAVE_HDF5_H
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf2_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf_serial_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf_to_one_ops);
#ifdef H5_HAVE_PARALLEL
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf_parallel_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf2_parallel_ops);
#endif
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_xdmf2_collective_ops);
#endif
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_ascii_ops);
  mrc_class_register_subclass(&mrc_class_mrc_io, &mrc_io_combined_ops);
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

