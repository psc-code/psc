
#include "psc.h"
#include "output_fields.h"

#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define OUTPUT_VTK_STRUCTURED_POINTS  1
#define OUTPUT_VTK_RECTILINEAR_POINTS 2
#define OUTPUT_VTK_RECTILINEAR_CELLS  3

#define OUTPUT_VTK OUTPUT_VTK_STRUCTURED_POINTS

// ======================================================================
// vtk_ctx

struct vtk_ctx {
  FILE *file;
};

static void
vtk_open(struct psc_fields_list *list, const char *filename, void **pctx)
{
  struct vtk_ctx *vtk = malloc(sizeof(*vtk));

  assert(list->nr_flds > 0);
  struct psc_field *fld = &list->flds[0];
  int *ilo = fld->ilo, *ihi = fld->ihi;

  vtk->file = fopen(filename, "w");
  fprintf(vtk->file, "# vtk DataFile Version 3.0\n");
  fprintf(vtk->file, "PSC fields timestep=%d dt=%g\n", psc.timestep, psc.dt);
  fprintf(vtk->file, "ASCII\n");
#ifdef OUTPUT_VTK == OUTPUT_VTK_STRUCTURED_POINTS
  fprintf(vtk->file, "DATASET STRUCTURED_POINTS\n");
#else
  int extra;
  double offset;
#if OUTPUT_VTK == OUTPUT_VTK_RECTILINEAR_POINTS
  extra = 0;
  offset = .5;
#else
  extra = 1;
  offset = 0.;
#endif
  fprintf(vtk->file, "DATASET RECTILINEAR_GRID\n");
  fprintf(vtk->file, "DIMENSIONS %d %d %d\n",
	  ihi[0] - ilo[0] + extra,
	  ihi[1] - ilo[1] + extra,
	  ihi[2] - ilo[2] + extra);
#endif

#ifdef OUTPUT_VTK == OUTPUT_VTK_STRUCTURED_POINTS

  fprintf(vtk->file, "SPACING %g %g %g\n", psc.dx[0], psc.dx[1], psc.dx[2]);
  fprintf(vtk->file, "ORIGIN %g %g %g\n",
	  psc.dx[0] * ilo[0],
	  psc.dx[1] * ilo[1],
	  psc.dx[2] * ilo[2]);
  fprintf(vtk->file, "\nPOINT_DATA %d\n", fld->size);

#else

  for (int d = 0; d < 3; d++) {
    fprintf(vtk->file, "%c_COORDINATES %d float", 'X' + d, ihi[d] - ilo[d]);
    for (int i = ilo[d]; i < ihi[d] + extra; i++) {
      fprintf(vtk->file, " %g", (i + offset) * psc.dx[d]);
    }
    fprintf(vtk->file, "\n");
  }
#if OUTPUT_VTK == OUTPUT_VTK_RECTILINEAR_POINTS
  fprintf(vtk->file, "\nPOINT_DATA %d\n", fld->size);
#else
  fprintf(vtk->file, "\nCELL_DATA %d\n", fld->size);
#endif

#endif


  *pctx = vtk;
}

static void
vtk_close(void *ctx)
{
  struct vtk_ctx *vtk = ctx;

  fclose(vtk->file);
}

static void
vtk_write_field(void *ctx, struct psc_field *fld)
{
  struct vtk_ctx *vtk = ctx;
  
  int mm[3];
  for (int d = 0; d < 3; d++) {
    mm[d] = fld->ihi[d] - fld->ilo[d];
  }

  fprintf(vtk->file, "SCALARS %s float\n", fld->name);
  fprintf(vtk->file, "LOOKUP_TABLE default\n");

  for (int iz = fld->ilo[2]; iz < fld->ihi[2]; iz++) {
    for (int iy = fld->ilo[1]; iy < fld->ihi[1]; iy++) {
      for (int ix = fld->ilo[0]; ix < fld->ihi[0]; ix++) {
	float val = fld->data[((iz-fld->ilo[2]) * mm[1] + iy-fld->ilo[1]) * mm[0] + ix-fld->ilo[0]];
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

// ======================================================================
// psc_output_format_ops_vtk

struct psc_output_format_ops psc_output_format_ops_vtk = {
  .name         = "vtk",
  .ext          = ".vtk",
  .open         = vtk_open,
  .close        = vtk_close,
  .write_field  = vtk_write_field,
};

