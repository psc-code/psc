
#include "psc.h"
#include "output_fields.h"

#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

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
  fprintf(vtk->file, "DATASET RECTILINEAR_GRID\n");
  fprintf(vtk->file, "DIMENSIONS %d %d %d\n",
	  ihi[0] - ilo[0],
	  ihi[1] - ilo[1],
	  ihi[2] - ilo[2]);

  for (int d = 0; d < 3; d++) {
    fprintf(vtk->file, "%c_COORDINATES %d float", 'X' + d, ihi[d] - ilo[d]);
    for (int i = ilo[d]; i < ihi[d]; i++) {
      fprintf(vtk->file, " %g", (i + .5) * psc.dx[d]);
    }
    fprintf(vtk->file, "\n");
  }
  fprintf(vtk->file, "\nPOINT_DATA %d\n", fld->size);

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
	fprintf(vtk->file, "%g ",
		fld->data[((iz-fld->ilo[2]) * mm[1] + iy-fld->ilo[1]) * mm[0] + ix-fld->ilo[0]]);
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

