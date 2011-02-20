
#include "psc.h"
#include "output_fields.h"

#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

//////////////////////////////////////////////////////////////////////
/// Helper to open the file and write the header.

static void
vtk_open_file(const char *pfx, const char *dataset_type, int extra,
	      struct psc_output_c *out, FILE **pfile)
{
  char filename[strlen(out->data_dir) + 30];
  sprintf(filename, "%s/%s_%07d.vtk", out->data_dir, pfx, psc.timestep);

  FILE *file = fopen(filename, "w");
  fprintf(file, "# vtk DataFile Version 3.0\n");
  fprintf(file, "PSC fields timestep=%d dt=%g\n", psc.timestep, psc.dt);
  fprintf(file, "ASCII\n");

  fprintf(file, "DATASET %s\n", dataset_type);
  fprintf(file, "DIMENSIONS %d %d %d\n",
	  psc.domain.gdims[0] + extra,
	  psc.domain.gdims[1] + extra,
	  psc.domain.gdims[2] + extra);
  
  *pfile = file;
}

//////////////////////////////////////////////////////////////////////
/// Helper to write the coordinates to the VTK file.

static void
vtk_write_coordinates(FILE *file, int extra, double offset)
{
  for (int d = 0; d < 3; d++) {
    fprintf(file, "%c_COORDINATES %d float", 'X' + d,
	    psc.domain.gdims[d] + extra);
    for (int i = 0; i < psc.domain.gdims[d] + extra; i++) {
      fprintf(file, " %g", (i + offset) * psc.dx[d]);
    }
    fprintf(file, "\n");
  }
}

//////////////////////////////////////////////////////////////////////
/// Helper to write one field to VTK file.

static void
vtk_write_field(void *ctx, fields_base_t *fld)
{
  FILE *file = ctx;

  int *ib = fld->ib, *im = fld->im;
  fprintf(file, "SCALARS %s float\n", fld->name[0]);
  fprintf(file, "LOOKUP_TABLE default\n");

  for (int iz = ib[2]; iz < ib[2] + im[2]; iz++) {
    for (int iy = ib[1]; iy < ib[1] + im[1]; iy++) {
      for (int ix = ib[0]; ix < ib[0] + im[0]; ix++) {
	float val = F3_BASE(fld, 0, ix,iy,iz);
	if (fabsf(val) < 1e-37) {
	  fprintf(file, "%g ", 0.);
	} else {
	  fprintf(file, "%g ", val);
	}
      }
      fprintf(file, "\n");
    }
  }
}

//////////////////////////////////////////////////////////////////////
/// Write VTK file for STRUCTURED_POINTS output.

static void
vtk_write_fields(struct psc_output_c *out, struct psc_fields_list *flds,
		 const char *pfx)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  FILE *file;
  if (rank == 0) {
    vtk_open_file(pfx, "STRUCTURED_POINTS", 0, out, &file);
    fprintf(file, "SPACING %g %g %g\n", psc.dx[0], psc.dx[1], psc.dx[2]);
    fprintf(file, "ORIGIN 0. 0. 0.\n");
    fprintf(file, "\nPOINT_DATA %d\n",
	    psc.domain.gdims[0] * psc.domain.gdims[1] * psc.domain.gdims[2]);
  }

  write_fields_combine(flds, vtk_write_field, file);

  if (rank == 0) {
    fclose(file);
  }
}

//////////////////////////////////////////////////////////////////////
/// Write VTK file for RECTILINEAR_GRID output with data values on points.

static void
vtk_points_write_fields(struct psc_output_c *out, struct psc_fields_list *flds,
			const char *pfx)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  FILE *file;
  if (rank == 0) {
    vtk_open_file(pfx, "RECTILINEAR_GRID", 0, out, &file);
    vtk_write_coordinates(file, 0, 0.);
    fprintf(file, "\nPOINT_DATA %d\n", fields_base_size(&psc.pf));
  }

  write_fields_combine(flds, vtk_write_field, file);

  if (rank == 0) {
    fclose(file);
  }
}

//////////////////////////////////////////////////////////////////////
/// Write VTK file for RECTILINEAR_GRID output with data values on cells.

static void
vtk_cells_write_fields(struct psc_output_c *out, struct psc_fields_list *flds,
		       const char *pfx)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  FILE *file;
  if (rank == 0) {
    vtk_open_file(pfx, "RECTILINEAR_GRID", 1, out, &file);
    vtk_write_coordinates(file, 1, .5);
    fprintf(file, "\nCELL_DATA %d\n", fields_base_size(&psc.pf));
  }
  write_fields_combine(flds, vtk_write_field, file);

  if (rank == 0) {
    fclose(file);
  }
}

//////////////////////////////////////////////////////////////////////
/// VTK output format writing ASCII STRUCTURED_POINTS file.

struct psc_output_format_ops psc_output_format_ops_vtk = {
  .name         = "vtk",
  .write_fields = vtk_write_fields,
};

//////////////////////////////////////////////////////////////////////
/// VTK output format writing ASCII RECTILINEAR_GRID file on points.
///
/// The grid points here are actually the centers of the computational cells.

struct psc_output_format_ops psc_output_format_ops_vtk_points = {
  .name         = "vtk_points",
  .write_fields = vtk_points_write_fields,
};

//////////////////////////////////////////////////////////////////////
/// VTK output format writing ASCII RECTILINEAR_GRID file on cells.
///
/// The grid points here are the true nodes of the computational cells.
/// The cell-centered field data is written as CELL_DATA, i.e. associated
/// with the cells, not the grid points.

struct psc_output_format_ops psc_output_format_ops_vtk_cells = {
  .name         = "vtk_cells",
  .write_fields = vtk_cells_write_fields,
};

