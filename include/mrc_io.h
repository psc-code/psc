
#ifndef MRC_IO_H
#define MRC_IO_H

#include <mrc_obj.h>
#include <mrc_domain.h>
#include <mrc_fld.h>

#include <stdbool.h>

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
void mrc_io_read_f1(struct mrc_io *io, const char *path, struct mrc_f1 *f1);
void mrc_io_read_f3(struct mrc_io *io, const char *path, struct mrc_f3 *f3);
void mrc_io_write_f1(struct mrc_io *io, const char *path, struct mrc_f1 *f1);
void mrc_io_write_f3(struct mrc_io *io, const char *path, struct mrc_f3 *f3, float scale);
void mrc_io_read_m1(struct mrc_io *io, const char *path, struct mrc_m1 *m1);
void mrc_io_write_m1(struct mrc_io *io, const char *path, struct mrc_m1 *m1);
void mrc_io_read_m3(struct mrc_io *io, const char *path, struct mrc_m3 *m3);
void mrc_io_write_m3(struct mrc_io *io, const char *path, struct mrc_m3 *m3);
void mrc_io_write_field2d(struct mrc_io *io, float scale, struct mrc_f2 *f2,
			  int outtype, float sheet);
void mrc_io_write_field_slice(struct mrc_io *io, float scale, struct mrc_f3 *f3,
			      int outtype, float sheet);

void mrc_io_read_attr(struct mrc_io *io, const char *path, int type, const char *name,
		      union param_u *pv);
void mrc_io_read_attr_int(struct mrc_io *io, const char *path, const char *name, int *val);
void mrc_io_read_attr_string(struct mrc_io *io, const char *path, const char *name,
			  char **pv);

void mrc_io_write_attr(struct mrc_io *io, const char *path, int type, const char *name,
		       union param_u *pv);
void mrc_io_write_attr_int(struct mrc_io *io, const char *path, const char *name, int val);
void mrc_io_write_attr_int3(struct mrc_io *io, const char *path, const char *name, int val[3]);
void mrc_io_write_attr_string(struct mrc_io *io, const char *path, const char *name,
			      const char *val);

// automate the cast to (struct mrc_class *)
#define mrc_io_read_obj_ref(io, path, name, class)			\
  __mrc_io_read_obj_ref(io, path, name, (struct mrc_class *)(class))
struct mrc_obj *__mrc_io_read_obj_ref(struct mrc_io *io, const char *path, const char *name,
				      struct mrc_class *class);
void mrc_io_write_obj_ref(struct mrc_io *io, const char *path, const char *name,
			  struct mrc_obj *obj);

// ----------------------------------------------------------------------
// for mrc_obj use

int mrc_io_add_obj(struct mrc_io *io, struct mrc_obj *obj);
struct mrc_obj *mrc_io_find_obj(struct mrc_io *io, const char *name);

// ----------------------------------------------------------------------
// mrc_io_server

void mrc_io_server(const char *format, const char *ds_srv, int nproc_domain);

#endif

