
#include "psc.h"
#include "output_fields.h"

#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

//////////////////////////////////////////////////////////////////////
/// Stores information private to the VTK writers.

struct vtk_ctx {
  FILE *file;
};

//////////////////////////////////////////////////////////////////////
/// Helper to open the file and write the header.

static void
vtk_open_file(struct vtk_ctx *vtk, const char *pfx, const char *dataset_type,
	      int extra, fields_base_t *fld, struct psc_output_c *out)
{
  char *filename = psc_output_c_filename(out, pfx);
  vtk->file = fopen(filename, "w");
  free(filename);
  fprintf(vtk->file, "# vtk DataFile Version 3.0\n");
  fprintf(vtk->file, "PSC fields timestep=%d dt=%g\n", psc.timestep, psc.dt);
  fprintf(vtk->file, "ASCII\n");

  fprintf(vtk->file, "DATASET %s\n", dataset_type);
  fprintf(vtk->file, "DIMENSIONS %d %d %d\n",
	  fld->im[0] + extra,
	  fld->im[1] + extra,
	  fld->im[2] + extra);
}

//////////////////////////////////////////////////////////////////////
/// Helper to write the coordinates to the VTK file.

static void
vtk_write_coordinates(struct vtk_ctx *vtk, int extra, double offset,
		      fields_base_t *fld)
{
  for (int d = 0; d < 3; d++) {
    fprintf(vtk->file, "%c_COORDINATES %d float", 'X' + d,
	    fld->im[d] + extra);
    for (int i = fld->ib[d]; i < fld->ib[d] + fld->im[d] + extra; i++) {
      fprintf(vtk->file, " %g", (i + offset) * psc.dx[d]);
    }
    fprintf(vtk->file, "\n");
  }
}

//////////////////////////////////////////////////////////////////////
/// Open VTK file for STRUCTURED_POINTS output.

static void
vtk_open(struct psc_output_c *out, struct psc_fields_list *list, const char *pfx,
	 void **pctx)
{
  struct vtk_ctx *vtk = malloc(sizeof(*vtk));

  assert(list->nr_flds > 0);
  fields_base_t *fld = &list->flds[0];
  int *ib = fld->ib;

  vtk_open_file(vtk, pfx, "STRUCTURED_POINTS", 0, fld, out);

  fprintf(vtk->file, "SPACING %g %g %g\n", psc.dx[0], psc.dx[1], psc.dx[2]);
  fprintf(vtk->file, "ORIGIN %g %g %g\n",
	  psc.dx[0] * ib[0],
	  psc.dx[1] * ib[1],
	  psc.dx[2] * ib[2]);
  fprintf(vtk->file, "\nPOINT_DATA %d\n", fields_base_size(fld));

  *pctx = vtk;
}

//////////////////////////////////////////////////////////////////////
/// Open VTK file for RECTILINEAR_GRID output with data values on points.

static void
vtk_points_open(struct psc_output_c *out, struct psc_fields_list *list,
		const char *pfx, void **pctx)
{
  struct vtk_ctx *vtk = malloc(sizeof(*vtk));

  assert(list->nr_flds > 0);
  fields_base_t *fld = &list->flds[0];
  
  vtk_open_file(vtk, pfx, "RECTILINEAR_GRID", 0, fld, out);
  vtk_write_coordinates(vtk, 0, 0., fld);
  fprintf(vtk->file, "\nPOINT_DATA %d\n", fields_base_size(fld));

  *pctx = vtk;
}

//////////////////////////////////////////////////////////////////////
/// Open VTK file for RECTILINEAR_GRID output with data values on cells.

static void
vtk_cells_open(struct psc_output_c *out, struct psc_fields_list *list,
	       const char *pfx, void **pctx)
{
  struct vtk_ctx *vtk = malloc(sizeof(*vtk));

  assert(list->nr_flds > 0);
  fields_base_t *fld = &list->flds[0];
  
  vtk_open_file(vtk, pfx, "RECTILINEAR_GRID", 1, fld, out);
  vtk_write_coordinates(vtk, 1, .5, fld);
  fprintf(vtk->file, "\nCELL_DATA %d\n", fields_base_size(fld));

  *pctx = vtk;
}

//////////////////////////////////////////////////////////////////////
/// Close VTK file.

static void
vtk_close(void *ctx)
{
  struct vtk_ctx *vtk = ctx;

  fclose(vtk->file);
}

//////////////////////////////////////////////////////////////////////
/// Write one field to VTK file.

static void
vtk_write_field(void *ctx, fields_base_t *fld)
{
  struct vtk_ctx *vtk = ctx;
  
  int *ib = fld->ib, *im = fld->im;
  fprintf(vtk->file, "SCALARS %s float\n", fld->name[0]);
  fprintf(vtk->file, "LOOKUP_TABLE default\n");

  for (int iz = ib[2]; iz < ib[2] + im[2]; iz++) {
    for (int iy = ib[1]; iy < ib[1] + im[1]; iy++) {
      for (int ix = ib[0]; ix < ib[0] + im[0]; ix++) {
	float val = XF3_BASE(fld, 0, ix,iy,iz);
	if (fabsf(val) < 1e-37) {
	  fprintf(vtk->file, "%g ", 0.);
	} else {
	  fprintf(vtk->file, "%g ", val);
	}
      }
      fprintf(vtk->file, "\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////
/// VTK output format writing ASCII STRUCTURED_POINTS file.

struct psc_output_format_ops psc_output_format_ops_vtk = {
  .name         = "vtk",
  .ext          = ".vtk",
  .open         = vtk_open,
  .close        = vtk_close,
  .write_field  = vtk_write_field,
};

//////////////////////////////////////////////////////////////////////
/// VTK output format writing ASCII RECTILINEAR_GRID file on points.
///
/// The grid points here are actually the centers of the computational cells.

struct psc_output_format_ops psc_output_format_ops_vtk_points = {
  .name         = "vtk_points",
  .ext          = ".vtk",
  .open         = vtk_points_open,
  .close        = vtk_close,
  .write_field  = vtk_write_field,
};

//////////////////////////////////////////////////////////////////////
/// VTK output format writing ASCII RECTILINEAR_GRID file on cells.
///
/// The grid points here are the true nodes of the computational cells.
/// The cell-centered field data is written as CELL_DATA, i.e. associated
/// with the cells, not the grid points.

struct psc_output_format_ops psc_output_format_ops_vtk_cells = {
  .name         = "vtk_cells",
  .ext          = ".vtk",
  .open         = vtk_cells_open,
  .close        = vtk_close,
  .write_field  = vtk_write_field,
};

