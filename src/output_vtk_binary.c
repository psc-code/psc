
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
vtk_open_file_binary(const char *pfx, const char *dataset_type, int extra,
	      struct psc_output_c *out, FILE **pfile)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);	
	
  char filename[strlen(out->data_dir) + 30];
  sprintf(filename, "%s/%07d_%s_%07d.vtk", out->data_dir, rank, pfx, psc.timestep);

  FILE *file = fopen(filename, "w");
  fprintf(file, "# vtk DataFile Version 3.0\n");
  fprintf(file, "PSC fields timestep=%d dt=%g\n", psc.timestep, psc.dt);
  fprintf(file, "ASCII\n");

  fprintf(file, "DATASET %s\n", dataset_type);
  /*fprintf(file, "DIMENSIONS %d %d %d\n",
	  psc.ihi[0] - psc.ilo[0] + extra,
	  psc.ihi[1] - psc.ilo[1] + extra,
	  psc.ihi[2] - psc.ilo[2] + extra);*/
	
  fprintf(file, "DIMENSIONS %d %d %d\n",
	  out->rx[0] - out->rn[0] + extra,
	  out->rx[1] - out->rn[1] + extra,
	  out->rx[2] - out->rn[2] + extra);
  
  *pfile = file;
}

#if 0
//////////////////////////////////////////////////////////////////////
/// Helper to write the coordinates to the VTK file.

static void
vtk_write_coordinates_binary(FILE *file, int extra, double offset)
{
  for (int d = 0; d < 3; d++) {
    fprintf(file, "%c_COORDINATES %d float", 'X' + d,
	    psc.domain.ihi[d] + extra);
    for (int i = 0; i < psc.domain.ihi[d] + extra; i++) {
      fprintf(file, " %g", (i + offset) * psc.dx[d]);
    }
    fprintf(file, "\n");
  }

}
#endif

//////////////////////////////////////////////////////////////////////
/// Helper to write one field to VTK file.

static void
vtk_write_field_binary(void *ctx, struct psc_mfields *flds, struct psc_output_c *out)
{
  FILE *file = ctx;

  fields_base_t *fld = &flds->f[0];
  int *ilo = out->rn, *ihi = out->rx;
	
  fprintf(file, "SCALARS %s float\n", fld->name[0]);
  fprintf(file, "LOOKUP_TABLE default\n");
//	fprintf(file, "rnx=%d rny=%d rnz=%d\n", ilo[0], ilo[1], ilo[2]);
//	fprintf(file, "rxx=%d ryx=%d rzx=%d\n", ihi[0], ihi[1], ihi[2]);
   
  for (int iz = ilo[2]; iz < ihi[2]; iz++) {
    for (int iy = ilo[1]; iy < ihi[1]; iy++) {
      for (int ix = ilo[0]; ix < ihi[0]; ix++) {
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
vtk_write_fields_binary(struct psc_output_c *out, struct psc_fields_list *flds,
		 const char *pfx)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  FILE *file;
    vtk_open_file_binary(pfx, "STRUCTURED_POINTS", 0, out, &file);
    fprintf(file, "SPACING %g %g %g\n", psc.dx[0], psc.dx[1], psc.dx[2]);
    fprintf(file, "ORIGIN %g %g %g\n",
	    psc.dx[0] * out->rn[0],
	    psc.dx[1] * out->rn[1],
	    psc.dx[2] * out->rn[2]);
    fprintf(file, "\nPOINT_DATA %d\n",
	    (-out->rn[0] + out->rx[0]) * 
	    (-out->rn[1] + out->rx[1]) *
	    (-out->rn[2] + out->rx[2]));


  
  for (int m = 0; m < flds->nr_flds; m++) 
  {	
	  vtk_write_field_binary(file, &flds->flds[m], out);
  }
  
  fclose(file);

}



//////////////////////////////////////////////////////////////////////
/// VTK output format writing binary STRUCTURED_POINTS files uncombined.

struct psc_output_format_ops psc_output_format_ops_vtk_binary = {
  .name         = "vtk_binary",
  .write_fields = vtk_write_fields_binary,
};

