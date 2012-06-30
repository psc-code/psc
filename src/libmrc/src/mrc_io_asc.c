
#include "mrc_io_private.h"
#include <mrc_a3.h>
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
  fprintf(file, "# %s\n", mrc_f3_comp_name(fld, m));

  const int *dims = mrc_f3_dims(fld);
  for (int iz = 0; iz < dims[2]; iz++) {
    for (int iy = 0; iy < dims[1]; iy++) {
      for (int ix = 0; ix < dims[0]; ix++) {
	fprintf(file, "%d %d %d %g\n", ix, iy, iz, MRC_F3(fld,m, ix,iy,iz));
      }
      fprintf(file, "\n");
    }
    fprintf(file, "\n");
  }
}

static void
ds_ascii_write_m3(struct mrc_io *io, const char *path, struct mrc_m3 *m3)
{
  struct mrc_io_ascii *ascii = to_mrc_io_ascii(io);
  
  struct mrc_io_params *par = &io->par;

  for (int p = 0; p < m3->nr_patches; p++) {
    char filename[strlen(par->outdir) + strlen(par->basename) + 30];
    fprintf(ascii->file, "# see %s\n", filename);
    sprintf(filename, "%s/%s.%06d_p%06d_%s.asc", par->outdir, par->basename,
	    io->step, p, mrc_m3_name(m3));
    FILE *file = fopen(filename, "w");
    fprintf(file, "# ix iy iz");
    for (int m = 0; m < m3->nr_comp; m++) {
      fprintf(file, " %s", mrc_m3_comp_name(m3, m));
    }
    fprintf(file, "\n");

    struct mrc_m3_patch *m3p = &m3->patches[p];
    mrc_m3_foreach(m3p, ix,iy,iz, 0,0) {
      fprintf(file, "%d %d %d", ix, iy, iz);
      for (int m = 0; m < m3->nr_comp; m++) {
	fprintf(file, " %g", MRC_M3(m3p, m, ix,iy,iz));
      }
      fprintf(file, "\n");
    } mrc_m3_foreach_end;

    fclose(file);
  }
}

static void
ds_ascii_write_a3(struct mrc_io *io, const char *path, struct mrc_a3 *a3)
{
  struct mrc_io_ascii *ascii = to_mrc_io_ascii(io);
  
  struct mrc_io_params *par = &io->par;

  for (int p = 0; p < a3->nr_patches; p++) {
    char filename[strlen(par->outdir) + strlen(par->basename) + 30];
    fprintf(ascii->file, "# see %s\n", filename);
    sprintf(filename, "%s/%s.%06d_p%06d_%s.asc", par->outdir, par->basename,
	    io->step, p, mrc_a3_name(a3));
    FILE *file = fopen(filename, "w");
    fprintf(file, "# ix iy iz");
    for (int m = 0; m < a3->nr_comp; m++) {
      fprintf(file, " %s", mrc_a3_comp_name(a3, m));
    }
    fprintf(file, "\n");

    struct mrc_patch_info info;
    mrc_domain_get_local_patch_info(a3->domain, p, &info);
    float xb[3], xe[3], dx[3];
    for (int d = 0; d < 3; d++) {
      xb[d] = (float) info.off[d] / (1 << info.level);
      xe[d] = (float) (info.off[d] + info.ldims[d]) / (1 << info.level);
      dx[d] = (xe[d] - xb[d]) / info.ldims[d];
    }
    struct mrc_a3_patch *a3p = &a3->patches[p];
    mrc_a3_foreach(a3p, ix,iy,iz, 0,0) {
      if (ix == 0) {
	fprintf(file, "\n");
      }
      float xx[3] = { xb[0] + ix * dx[0], xb[1] + iy * dx[1], xb[2] + iz * dx[2] };
      fprintf(file, "%g %g %g", xx[0], xx[1], xx[2]);
      for (int m = 0; m < a3->nr_comp; m++) {
	fprintf(file, " %g", MRC_A3(a3p, m, ix,iy,iz));
      }
      fprintf(file, "\n");
    } mrc_a3_foreach_end;

    fclose(file);
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

struct mrc_io_ops mrc_io_ascii_ops = {
  .name        = "ascii",
  .size        = sizeof(struct mrc_io_ascii),
  .open        = ds_ascii_open,
  .close       = ds_ascii_close,
  .write_field = ds_ascii_write_field,
  .write_m3    = ds_ascii_write_m3,
  .write_a3    = ds_ascii_write_a3,
  .write_attr  = ds_ascii_write_attr,
};
