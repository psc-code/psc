
#ifndef MRC_DIAG_H
#define MRC_DIAG_H

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

struct diag_info {
  const char *outdir, *run;
};

extern struct mrc_class mrc_class_diag_format;

MRC_OBJ_DEFINE_STANDARD_METHODS(diag_format, struct diag_format);
void diag_format_open(struct diag_format *format, float sheet, int outtype, int step,
		      float time, const char *time_str);
void diag_format_close(struct diag_format *format);
void diag_format_write_field(struct diag_format *format, float scale,
			     struct mrc_f3 *f3, int m);
void diag_format_write_field2d(struct diag_format *format, float scale, struct mrc_f2 *f2,
			       const char *fld_name, int outtype, float sheet);
void diag_format_write_field_slice(struct diag_format *format, float scale, struct mrc_f3 *f3,
				   const char *fld_name, int outtype, float sheet);
// simple wrapper to write multiple components in one call
void diag_format_write_fields(struct diag_format *format, struct mrc_f3 *f, int mm[]);

void diagsrv_one(const char *format, const char *ds_srv, int nproc_domain);

#endif

