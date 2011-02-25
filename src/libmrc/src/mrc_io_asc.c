
#include "mrc_io_private.h"
#include <mrc_params.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

// ----------------------------------------------------------------------
// diag server ascii

struct mrc_io_ascii {
  FILE *file;
};

#define to_mrc_io_ascii(io) ((struct mrc_io_ascii *)io->obj.subctx)

static void
ds_ascii_open(struct mrc_io *io, const char *mode) 
{
  assert(strcmp(mode, "w") == 0); // only writing supported

  struct mrc_io_ascii *ascii = to_mrc_io_ascii(io);
  char *filename = diagc_make_filename(io, ".asc");
  ascii->file = fopen(filename, mode);
  free(filename);
}

static void
ds_ascii_close(struct mrc_io *io)
{
  struct mrc_io_ascii *ascii = to_mrc_io_ascii(io);
  fclose(ascii->file);
  ascii->file = NULL;
}

static void
ds_ascii_write_field(struct mrc_io *io, const char *path,
		     float scale, struct mrc_f3 *fld, int m)
{
  struct mrc_io_ascii *ascii = to_mrc_io_ascii(io);
  FILE *file = ascii->file;
  fprintf(file, "# %s\n", fld->name[m]);

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

static void
ds_ascii_write_attr(struct mrc_io *io, const char *path, int type,
		    const char *name, union param_u *pv)
{
  struct mrc_io_ascii *ascii = to_mrc_io_ascii(io);
  FILE *file = ascii->file;

  switch (type) {
  case PT_INT:
  case PT_SELECT:
    fprintf(file, "# %-20s| %d\n", name, pv->u_int);
    break;
  case PT_BOOL:
    fprintf(file, "# %-20s| %s\n", name, pv->u_bool ? "yes" : "no");
    break;
  case PT_FLOAT:
    fprintf(file, "# %-20s| %g\n", name, pv->u_float);
    break;
  case PT_DOUBLE:
    fprintf(file, "# %-20s| %g\n", name, pv->u_double);
    break;
  case PT_STRING:
    fprintf(file, "# %-20s| %s\n", name, pv->u_string);
    break;
  }
}

static struct mrc_io_ops mrc_io_ops_ascii = {
  .name        = "ascii",
  .size        = sizeof(struct mrc_io_ascii),
  .open        = ds_ascii_open,
  .close       = ds_ascii_close,
  .write_field = ds_ascii_write_field,
  .write_attr  = ds_ascii_write_attr,
};

void
libmrc_io_register_ascii()
{
  libmrc_io_register(&mrc_io_ops_ascii);
}
