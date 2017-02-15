
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
  int step; //> current output's timestep (for use in creating filename)
  float time; //> current output's time

  int rank; //> this proc's rank in comm
  int size; //> nr procs in comm
  bool diagc_domain_info_sent;
  list_t obj_list; //> objects written / read, lives open() -> close()
};

struct mrc_io_ops {
  MRC_SUBCLASS_OPS(struct mrc_io);
  bool parallel;
  void (*open)(struct mrc_io *, const char *mode);
  void (*close)(struct mrc_io *);
  void (*read_fld)(struct mrc_io *, const char *path, struct mrc_fld *fld);
  void (*read_f3)(struct mrc_io *, const char *path, struct mrc_fld *fld);
  void (*write_fld)(struct mrc_io *, const char *path, struct mrc_fld *fld);
  void (*read_m1)(struct mrc_io *, const char *path, struct mrc_fld *m1);
  void (*write_m1)(struct mrc_io *, const char *path, struct mrc_fld *m1);
  void (*read_m3)(struct mrc_io *, const char *path, struct mrc_m3 *m3);
  void (*write_m3)(struct mrc_io *, const char *path, struct mrc_m3 *m3);
  void (*write_field)(struct mrc_io *, const char *path,
		      float scale, struct mrc_fld *fld, int m);
  void (*write_field2d)(struct mrc_io *, float scale, struct mrc_fld *fld,
			int outtype, float sheet);
  void (*write_ndarray)(struct mrc_io *, const char *path, struct mrc_ndarray *nd);
  void (*read_ndarray)(struct mrc_io *, const char *path, struct mrc_ndarray *nd);
  void (*read_attr)(struct mrc_io *io, const char *path, int type,
		    const char *name, union param_u *pv);
  void (*write_attr)(struct mrc_io *io, const char *path, int type,
		     const char *name, union param_u *pv);
  void (*get_h5_file)(struct mrc_io *io, long *h5_file);
};

// ----------------------------------------------------------------------
// library functions for mrc_io's to use

char *diagc_make_filename(struct mrc_io *io, const char *ext);

// ----------------------------------------------------------------------
// diagsrv_one

extern struct diagsrv_srv_ops ds_srv_ops;
extern struct diagsrv_srv_ops ds_srv_cache_ops;

extern struct mrc_io_ops mrc_io_ascii_ops;
extern struct mrc_io_ops mrc_io_vpic_ops;
extern struct mrc_io_ops mrc_io_xdmf_ops;
extern struct mrc_io_ops mrc_io_xdmf_serial_ops;
extern struct mrc_io_ops mrc_io_xdmf_to_one_ops;
extern struct mrc_io_ops mrc_io_xdmf_parallel_ops;
extern struct mrc_io_ops mrc_io_xdmf2_ops;
extern struct mrc_io_ops mrc_io_xdmf2_parallel_ops;
extern struct mrc_io_ops mrc_io_hdf5_parallel_ops;
extern struct mrc_io_ops mrc_io_xdmf_collective_ops;
extern struct mrc_io_ops mrc_io_hdf5_serial_ops;
extern struct mrc_io_ops mrc_io_combined_ops;

#endif

