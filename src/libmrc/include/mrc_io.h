
#ifndef MRC_IO_H
#define MRC_IO_H

#include <mrc_obj.h>
#include <mrc_domain.h>
#include <mrc_fld.h>

#include <stdbool.h>

BEGIN_C_DECLS

// ----------------------------------------------------------------------
// diagnostics writing

// These constants correspond to the integer tags for the planes in diagwr2.
enum {
  DIAG_TYPE_3D      = 1,
  DIAG_TYPE_2D_X    = 2,
  DIAG_TYPE_2D_Y    = 3,
  DIAG_TYPE_2D_Z    = 4,
  DIAG_TYPE_2D_IONO = 5,
  NR_DIAG_TYPE,
};

union param_u;

MRC_CLASS_DECLARE(mrc_io, struct mrc_io);

void mrc_io_open(struct mrc_io *io, const char *mode, int step, float time);
void mrc_io_close(struct mrc_io *io);
void mrc_io_read_fld(struct mrc_io *io, const char *path, struct mrc_fld *fld);
void mrc_io_read_f3(struct mrc_io *io, const char *path, struct mrc_fld *fld);
void mrc_io_write_fld(struct mrc_io *io, const char *path, struct mrc_fld *fld);
void mrc_io_read_m1(struct mrc_io *io, const char *path, struct mrc_fld *m1);
void mrc_io_write_m1(struct mrc_io *io, const char *path, struct mrc_fld *m1);
void mrc_io_read_m3(struct mrc_io *io, const char *path, struct mrc_m3 *m3);
void mrc_io_write_m3(struct mrc_io *io, const char *path, struct mrc_m3 *m3);
void mrc_io_write_field2d(struct mrc_io *io, float scale, struct mrc_fld *fld,
			  int outtype, float sheet);
void mrc_io_write_field_slice(struct mrc_io *io, float scale, struct mrc_fld *fld,
			      int outtype, float sheet);
void mrc_io_write_ndarray(struct mrc_io *io, const char *path, struct mrc_ndarray *nd);
void mrc_io_read_ndarray(struct mrc_io *io, const char *path, struct mrc_ndarray *nd);

void mrc_io_read_attr(struct mrc_io *io, const char *path, int type, const char *name,
		      union param_u *pv);
void mrc_io_read_attr_int(struct mrc_io *io, const char *path, const char *name, int *val);
void mrc_io_read_attr_int3(struct mrc_io *io, const char *path, const char *name, int (*val)[3]);
void mrc_io_read_attr_float(struct mrc_io *io, const char *path, const char *name, float *val);
void mrc_io_read_attr_float3(struct mrc_io *io, const char *path, const char *name, float (*val)[3]);
void mrc_io_read_attr_double(struct mrc_io *io, const char *path, const char *name, double *val);
void mrc_io_read_attr_double3(struct mrc_io *io, const char *path, const char *name, double (*val)[3]);
void mrc_io_read_attr_bool(struct mrc_io *io, const char *path, const char *name, bool *pv);
void mrc_io_read_attr_string(struct mrc_io *io, const char *path, const char *name,
			  char **pv);

void mrc_io_write_attr(struct mrc_io *io, const char *path, int type, const char *name,
		       union param_u *pv);
void mrc_io_write_attr_int(struct mrc_io *io, const char *path, const char *name, int val);
void mrc_io_write_attr_int3(struct mrc_io *io, const char *path, const char *name, int val[3]);
void mrc_io_write_attr_float(struct mrc_io *io, const char *path, const char *name, float val);
void mrc_io_write_attr_float3(struct mrc_io *io, const char *path, const char *name, float val[3]);
void mrc_io_write_attr_double(struct mrc_io *io, const char *path, const char *name, double val);
void mrc_io_write_attr_double3(struct mrc_io *io, const char *path, const char *name, double val[3]);
void mrc_io_write_attr_bool(struct mrc_io *io, const char *path, const char *name, const bool val);
void mrc_io_write_attr_string(struct mrc_io *io, const char *path, const char *name,
			      const char *val);

#define mrc_io_read_int(io, obj, name, val) \
  mrc_io_read_attr_int(io, mrc_io_obj_path(io, obj), name, val)
#define mrc_io_read_int3(io, obj, name, val) \
  mrc_io_read_attr_int3(io, mrc_io_obj_path(io, obj), name, val)
#define mrc_io_read_double(io, obj, name, val) \
  mrc_io_read_attr_double(io, mrc_io_obj_path(io, obj), name, val)
#define mrc_io_read_string(io, obj, name, val) \
  mrc_io_read_attr_string(io, mrc_io_obj_path(io, obj), name, val)

#define mrc_io_write_int(io, obj, name, val) \
  mrc_io_write_attr_int(io, mrc_io_obj_path(io, obj), name, val)
#define mrc_io_write_int3(io, obj, name, val) \
  mrc_io_write_attr_int3(io, mrc_io_obj_path(io, obj), name, val)
#define mrc_io_write_double(io, obj, name, val) \
  mrc_io_write_attr_double(io, mrc_io_obj_path(io, obj), name, val)
#define mrc_io_write_string(io, obj, name, val) \
  mrc_io_write_attr_string(io, mrc_io_obj_path(io, obj), name, val)

void mrc_io_get_h5_file(struct mrc_io *io, long *h5_file);

struct mrc_obj *__mrc_io_read_path(struct mrc_io *io, const char *path, const char *name,
				   struct mrc_class *cls);
struct mrc_obj *__mrc_io_read_ref(struct mrc_io *io, struct mrc_obj *obj_parent,
				  const char *name, struct mrc_class *cls);
struct mrc_obj *__mrc_io_read_ref_comm(struct mrc_io *io, struct mrc_obj *obj_parent,
				       const char *name, struct mrc_class *cls,
				       MPI_Comm comm);
void __mrc_io_write_path(struct mrc_io *io, const char *path, const char *name,
			 struct mrc_obj *obj);
void __mrc_io_write_ref(struct mrc_io *io, struct mrc_obj *obj_parent, const char *name,
			struct mrc_obj *obj);

// hide the casts
#define mrc_io_read_path(io, path, name, cls)				\
  (struct cls *)							\
  __mrc_io_read_path(io, path, name,					\
		     (struct mrc_class *)(&mrc_class_ ## cls))
#define mrc_io_read_ref(io, obj, name, cls)				\
  (struct cls *)							\
  __mrc_io_read_ref(io, (struct mrc_obj *) obj, name,			\
		    (struct mrc_class *)(&mrc_class_ ## cls))
#define mrc_io_read_ref_comm(io, obj, name, cls, comm)			\
  (struct cls *)							\
  __mrc_io_read_ref_comm(io, (struct mrc_obj *) obj, name,		\
			 (struct mrc_class *)(&mrc_class_ ## cls), comm)
#define mrc_io_write_path(io, path, name, obj)			\
  __mrc_io_write_path(io, path, name, (struct mrc_obj *) obj)
#define mrc_io_write_ref(io, obj_parent, name, obj)		\
  __mrc_io_write_ref(io, (struct mrc_obj *) obj_parent,		\
		     name, (struct mrc_obj *) obj)

// ----------------------------------------------------------------------
// for mrc_obj use

int mrc_io_add_obj(struct mrc_io *io, struct mrc_obj *obj, const char *path);
struct mrc_obj *mrc_io_find_obj(struct mrc_io *io, const char *path);
const char *__mrc_io_obj_path(struct mrc_io *io, struct mrc_obj *obj);

#define mrc_io_obj_path(io, obj) \
  __mrc_io_obj_path(io, (struct mrc_obj *) obj)

// ----------------------------------------------------------------------
// mrc_io_server

void mrc_io_server(const char *format, const char *ds_srv, int nproc_domain);

END_C_DECLS

#endif

