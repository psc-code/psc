
#ifndef MRC_DIAG_PRIVATE_H
#define MRC_DIAG_PRIVATE_H

#include <mrc_io.h>
#include <mrc_list.h>

#include <assert.h>

// ----------------------------------------------------------------------
// mrc_io internals

struct mrc_io_params {
  char *outdir;
  char *basename; //> filename w/o dir, timestep, suffix
};

struct mrc_io {
  struct mrc_obj obj;
  struct mrc_io_params par;
  bool is_setup;
  int step; //> current output's timestep (for use in creating filename)
  float time; //> current output's time

  int rank; //> this proc's rank in comm
  int size; //> nr procs in comm
  bool diagc_domain_info_sent;
  list_t obj_list; //> objects written / read, lives open() -> close()
};

struct mrc_io_ops {
  MRC_OBJ_OPS;
  bool parallel;
  void (*open)(struct mrc_io *, const char *mode);
  void (*close)(struct mrc_io *);
  void (*read_f1)(struct mrc_io *, const char *path, struct mrc_f1 *fld);
  void (*read_f3)(struct mrc_io *, const char *path, struct mrc_f3 *fld);
  void (*write_f1)(struct mrc_io *, const char *path, struct mrc_f1 *fld);
  void (*write_f3)(struct mrc_io *, const char *path,
		   struct mrc_f3 *fld, float scale);
  void (*write_m3)(struct mrc_io *, const char *path, struct mrc_m3 *m3);
  void (*write_field)(struct mrc_io *, const char *path,
		      float scale, struct mrc_f3 *fld, int m);
  void (*write_field2d)(struct mrc_io *, float scale, struct mrc_f2 *fld,
			int outtype, float sheet);
  void (*read_attr)(struct mrc_io *io, const char *path, int type,
		    const char *name, union param_u *pv);
  void (*write_attr)(struct mrc_io *io, const char *path, int type,
		     const char *name, union param_u *pv);
};

// ----------------------------------------------------------------------
// library functions for mrc_io's to use

char *diagc_make_filename(struct mrc_io *io, const char *ext);

// ----------------------------------------------------------------------
// diagsrv_one

extern struct diagsrv_srv_ops ds_srv_ops;
extern struct diagsrv_srv_ops ds_srv_cache_ops;

void libmrc_io_register_ascii();
void libmrc_io_register_xdmf();
void libmrc_io_register_xdmf2();
void libmrc_io_register_combined();

void libmrc_io_register(struct mrc_io_ops *ops);

// ======================================================================

static inline struct mrc_io *
to_mrc_io(struct mrc_obj *obj)
{
  assert(obj->class == &mrc_class_mrc_io);
  return container_of(obj, struct mrc_io, obj);
}


#endif

