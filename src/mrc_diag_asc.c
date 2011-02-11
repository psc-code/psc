
#include "mrc_diag_private.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

// ----------------------------------------------------------------------
// diag server ascii

static void
ds_ascii_open(struct diag_format *format, float sheet, int outtype, int step) 
{
  char *filename = diagc_make_filename(format, ".asc", sheet, outtype, step);
  format->obj.subctx = fopen(filename, "w");
  free(filename);
}

static void
ds_ascii_close(struct diag_format *format)
{
  FILE *file = format->obj.subctx;
  fclose(file);
  format->obj.subctx = NULL;
}

static void
ds_ascii_write_field(struct diag_format *format, float scale, struct mrc_f3 *fld,
		     int m)
{
  FILE *file = format->obj.subctx;
  fprintf(file, "# %s\n", fld->name[m]);
  fprintf(file, "# %s\n", format->time_str);

  for (int iz = 0; iz < fld->im[2]; iz++) {
    for (int iy = 0; iy < fld->im[1]; iy++) {
      for (int ix = 0; ix < fld->im[0]; ix++) {
	fprintf(file, "%d %d %d %g\n", ix, iy, iz, MRC_F3(fld,m, ix,iy,iz));
      }
      fprintf(file, "\n");
    }
    fprintf(file, "\n");
  }
}

static struct diag_format_ops ds_ascii_ops = {
  .name        = "ascii",
  .open        = ds_ascii_open,
  .close       = ds_ascii_close,
  .write_field = ds_ascii_write_field,
};

void
libmrc_diag_ascii_register()
{
  libmrc_diag_register_format(&ds_ascii_ops);
}
