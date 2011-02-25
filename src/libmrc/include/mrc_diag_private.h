
#ifndef MRC_DIAG_PRIVATE_H
#define MRC_DIAG_PRIVATE_H

#include <mrc_diag.h>
#include <mrc_list.h>

#include <assert.h>

// ----------------------------------------------------------------------
// diag_format internals

struct diag_format {
  struct mrc_obj obj;
  struct diag_info diag_info;
  bool is_setup;
  char *time_str;
  float time;
  int rank; //> this procs rank in comm
  bool diagc_domain_info_sent;
};

struct diag_format_ops {
  MRC_OBJ_OPS;
  bool parallel;
  void (*open)(struct diag_format *, float sheet, int outtype, int step);
  void (*close)(struct diag_format *);
  void (*write_field)(struct diag_format *, float scale, struct mrc_f3 *fld, int m);
  void (*write_field2d)(struct diag_format *, float scale, struct mrc_f2 *fld,
			const char *fld_name, int outtype, float sheet);
};

// ----------------------------------------------------------------------
// library functions for diag_formats to use

char *diagc_make_filename(struct diag_format *format, const char *sfx,
			  float sheet, int outtype, int step);

// ----------------------------------------------------------------------
// diagsrv_one

extern struct diagsrv_srv_ops ds_srv_ops;
extern struct diagsrv_srv_ops ds_srv_cache_ops;

void libmrc_diag_ascii_register();
void libmrc_diag_hdf5_register();
void libmrc_diag_combined_register();

void libmrc_diag_register_format(struct diag_format_ops *ops);

// ======================================================================

static inline struct diag_format *
to_diag_format(struct mrc_obj *obj)
{
  assert(obj->class == &mrc_class_diag_format);
  return container_of(obj, struct diag_format, obj);
}


#endif

