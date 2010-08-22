
#include "psc.h"
#include "output_fields.h"

#include "util/profile.h"
#include "util/params.h"

#include <mpi.h>
#include <string.h>
#include <hdf5.h>
#include <hdf5_hl.h>

// ======================================================================
// hdf5_ctx

struct hdf5_ctx {
  hid_t file;
  hid_t group;
  hid_t group_fld;
};

static void
hdf5_open(struct psc_output_c *out, struct psc_fields_list *list, const char *pfx,
	  void **pctx)
{
  struct hdf5_ctx *hdf5 = malloc(sizeof(*hdf5));

  char *filename = psc_output_c_filename(out, pfx);
  hdf5->file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  free(filename);

  hdf5->group = H5Gcreate(hdf5->file, "psc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_int(hdf5->group, ".", "timestep", &psc.timestep, 1);
  H5LTset_attribute_double(hdf5->group, ".", "dt", &psc.dt, 1);
  H5LTset_attribute_double(hdf5->group, ".", "dx", psc.dx, 3);
  H5LTset_attribute_int(hdf5->group, ".", "lo", psc.domain.ilo, 3);
  H5LTset_attribute_int(hdf5->group, ".", "hi", psc.domain.ihi, 3);

  hdf5->group_fld = H5Gcreate(hdf5->group, "fields",
			      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  *pctx = hdf5;
}

static void
hdf5_close(void *ctx)
{
  struct hdf5_ctx *hdf5 = ctx;

  H5Gclose(hdf5->group_fld);
  H5Gclose(hdf5->group);
  H5Fclose(hdf5->file);
}

static void
hdf5_write_field(void *ctx, struct psc_field *fld)
{
  struct hdf5_ctx *hdf5 = ctx;

  hsize_t dims[3];
  for (int d = 0; d < 3; d++) {
    // reverse dimensions because of Fortran order
    dims[d] = fld->ihi[2-d] - fld->ilo[2-d];
  }
  
  if (sizeof(*fld->data) == 4) {
    H5LTmake_dataset_float(hdf5->group_fld, fld->name, 3, dims, (float *) fld->data);
  } else {
    H5LTmake_dataset_double(hdf5->group_fld, fld->name, 3, dims, (double *) fld->data);
  }
  H5LTset_attribute_int(hdf5->group_fld, fld->name, "lo", fld->ilo, 3);
  H5LTset_attribute_int(hdf5->group_fld, fld->name, "hi", fld->ihi, 3);
}

// ======================================================================
// psc_output_format_ops_hdf5

struct psc_output_format_ops psc_output_format_ops_hdf5 = {
  .name         = "hdf5",
  .ext          = ".h5",
  .open         = hdf5_open,
  .close        = hdf5_close,
  .write_field  = hdf5_write_field,
};

// ======================================================================

static void
xdmf_open(struct psc_output_c *out, struct psc_fields_list *list, const char *pfx,
	  void **pctx)
{
  assert(!out->output_combine);

  hdf5_open(out, list, pfx, pctx);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    char fname[strlen(out->data_dir) + strlen(pfx) + 20];
    sprintf(fname, "%s/%s_%07d.xdmf", out->data_dir, pfx, psc.timestep);
    FILE *f = fopen(fname, "w");

    struct psc_field *fld = &list->flds[0];
    fprintf(f, "<?xml version=\"1.0\" ?>\n");
    fprintf(f, "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n");
    fprintf(f, "<Domain>\n");
    fprintf(f, "<Grid GridType='Collection' CollectionType='Spatial'>\n");
    fprintf(f, "   <Time Type=\"Single\" Value=\"%g\" />\n", psc.timestep * psc.dt);
    for (int kz = 0; kz < psc.domain.nproc[2]; kz++) {
      for (int ky = 0; ky < psc.domain.nproc[1]; ky++) {
	for (int kx = 0; kx < psc.domain.nproc[0]; kx++) {
	  fprintf(f, "   <Grid Name=\"mesh-%d-%d-%d-%d\" GridType=\"Uniform\">\n",
		  kx, ky, kz, psc.timestep);
	  fprintf(f, "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\"/>\n",
		  fld->ihi[2] - fld->ilo[2] + 1,
		  fld->ihi[1] - fld->ilo[1] + 1,
		  fld->ihi[0] - fld->ilo[0] + 1);
	  fprintf(f, "     <Geometry GeometryType=\"Origin_DxDyDz\">\n");
	  fprintf(f, "     <DataStructure Name=\"Origin\" DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">\n");
	  fprintf(f, "        %g %g %g\n",
		  (psc.domain.ilo[2] + (psc.domain.ihi[2] - psc.domain.ilo[2]) / psc.domain.nproc[2] * kz) * psc.dx[2],
		  (psc.domain.ilo[1] + (psc.domain.ihi[1] - psc.domain.ilo[1]) / psc.domain.nproc[1] * ky) * psc.dx[1],
		  (psc.domain.ilo[0] + (psc.domain.ihi[0] - psc.domain.ilo[0]) / psc.domain.nproc[0] * kx) * psc.dx[0]);
	  fprintf(f, "     </DataStructure>\n");
	  fprintf(f, "     <DataStructure Name=\"Spacing\" DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">\n");
	  fprintf(f, "        %g %g %g\n", psc.dx[2], psc.dx[1], psc.dx[0]);
	  fprintf(f, "     </DataStructure>\n");
	  fprintf(f, "     </Geometry>\n");
	  fprintf(f, "\n");
	  for (int m = 0; m < list->nr_flds; m++) {
	    fld = &list->flds[m];
	    fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
		    fld->name);
	    fprintf(f, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
		    fld->ihi[2] - fld->ilo[2],
		    fld->ihi[1] - fld->ilo[1],
		    fld->ihi[0] - fld->ilo[0]);
	    int proc = (kz * psc.domain.nproc[1] + ky) * psc.domain.nproc[0] + kx;
	    fprintf(f, "        %s_%06d_%07d.h5:/psc/fields/%s\n",
		    pfx, proc, psc.timestep, fld->name);
	    fprintf(f, "       </DataItem>\n");
	    fprintf(f, "     </Attribute>\n");
 	  }
	  fprintf(f, "   </Grid>\n");
	}
      }
    }
    fprintf(f, "</Grid>\n");
    fprintf(f, "</Domain>\n");
    fprintf(f, "</Xdmf>\n");
    fclose(f);


    // It'd be easier to create those files line by line as time goes by.
    // However, then we won't be able to get the timeseries into Paraview
    // until the solution is all finished.
    // So this version rewrites the xdmf file completely every timestep,
    // which needs however some way to figure out what times we've written
    // before.
    sprintf(fname, "%s.xdmf", pfx);
    f = fopen(fname, "w");
    fprintf(f, "<?xml version='1.0' ?>\n");
    fprintf(f, "<Xdmf xmlns:xi='http://www.w3.org/2001/XInclude' Version='2.0'>\n");
    fprintf(f, "<Domain>\n");
    fprintf(f, "  <Grid GridType='Collection' CollectionType='Temporal'>\n");
    int first, step;
    if (strcmp(pfx, "pfd") == 0) {
      first = out->pfield_first;
      step = out->pfield_step;
    } else {
      first = out->tfield_first;
      step = out->tfield_step;
    }
    for (int i = first; i <= psc.timestep; i += step) {
      fprintf(f, "  <xi:include href='%s_%07d.xdmf' xpointer='xpointer(//Xdmf/Domain/Grid)'/>\n", pfx, i);
    }
    fprintf(f, "  </Grid>\n");
    fprintf(f, "  </Domain>\n");
    fprintf(f, "</Xdmf>\n");
    fclose(f);
  }
}

// ======================================================================
// psc_output_format_ops_xdmf

struct psc_output_format_ops psc_output_format_ops_xdmf = {
  .name         = "xdmf",
  .ext          = ".h5",
  .open         = xdmf_open,
  .close        = hdf5_close,
  .write_field  = hdf5_write_field,
};

// ======================================================================
// hdf5_dump_field
//
// dumps the local fields to a file

// FIXME still assuming fortran layout for buf

static void
hdf5_dump_field(int m, const char *fname)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  char *filename = malloc(strlen(fname) + 10);
  sprintf(filename, "%s-p%d.h5", fname, rank);
  mpi_printf(MPI_COMM_WORLD, "hdf5_dump_field: '%s'\n", filename);
  
  hsize_t dims[3] = { psc.ihg[2] - psc.ilg[2],
		      psc.ihg[1] - psc.ilg[1],
		      psc.ihg[0] - psc.ilg[0] };

  hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hid_t group = H5Gcreate(file, "psc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t group_fld = H5Gcreate(group, "fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(sizeof(fields_base_real_t) == sizeof(double));
  H5LTmake_dataset_double(group_fld, fldname[m], 3, dims,
			  &F3_BASE(m, psc.ilg[0], psc.ilg[1], psc.ilg[2]));
  H5Gclose(group_fld);
  H5Gclose(group);
  H5Fclose(file);
}
