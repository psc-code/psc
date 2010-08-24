
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
__hdf5_open(struct psc_output_c *out, struct psc_fields_list *list,
	    const char *filename, void **pctx)
{
  struct hdf5_ctx *hdf5 = malloc(sizeof(*hdf5));

  hdf5->file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

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
__hdf5_close(void *ctx)
{
  struct hdf5_ctx *hdf5 = ctx;

  H5Gclose(hdf5->group_fld);
  H5Gclose(hdf5->group);
  H5Fclose(hdf5->file);
}

static void
__hdf5_write_field(void *ctx, fields_base_t *fld)
{
  struct hdf5_ctx *hdf5 = ctx;

  hid_t mem_type;
  if (sizeof(*fld->flds) == 4) { // FIXME
    mem_type = H5T_NATIVE_FLOAT;
  } else {
    mem_type = H5T_NATIVE_DOUBLE;
  }

  hsize_t file_dims[3], mem_dims[3];
  int ie[3];
  for (int d = 0; d < 3; d++) {
    // reverse dimensions because of Fortran order
    mem_dims[d] = fld->im[2-d];
    ie[d] = fld->ib[d] + fld->im[d];
  }
  hid_t mem_space = H5Screate_simple(3, mem_dims, NULL);
  
  if (fld->im[0] == psc.img[0] &&
      fld->im[1] == psc.img[1] &&
      fld->im[2] == psc.img[2]) {
    // we're writing the local field, let's drop the ghost points
    hsize_t start[3], count[3];
    for (int d = 0; d < 3; d++) {
      // reverse dimensions because of Fortran order
      file_dims[d] = psc.ihi[2-d] - psc.ilo[2-d];
      start[d] = psc.ilo[2-d] - psc.ilg[2-d];
      count[d] = psc.ihi[2-d] - psc.ilo[2-d];
    }
    H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, start, NULL, count, NULL);
  } else {
    for (int d = 0; d < 3; d++) {
      file_dims[d] = mem_dims[d];
    }
  }

  hid_t file_space = H5Screate_simple(3, file_dims, NULL); 
  hid_t dataset = H5Dcreate(hdf5->group_fld, fld->name[0], H5T_NATIVE_FLOAT,
			    file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Sclose(file_space);
  
  H5Dwrite(dataset, mem_type, mem_space, H5S_ALL, H5P_DEFAULT, fld->flds);
  H5Sclose(mem_space);
  
  H5Dclose(dataset);

  H5LTset_attribute_int(hdf5->group_fld, fld->name[0], "lo", fld->ib, 3);
  H5LTset_attribute_int(hdf5->group_fld, fld->name[0], "hi", ie, 3);
}

// ======================================================================

static void
hdf5_open(struct psc_output_c *out, struct psc_fields_list *list, const char *pfx,
	  void **pctx)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    char filename[strlen(out->data_dir) + 30];
    sprintf(filename, "%s/%s_%07d.h5", out->data_dir, pfx, psc.timestep);
    __hdf5_open(out, list, filename, pctx);
  }
}

static void
hdf5_close(void *ctx)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    __hdf5_close(ctx);
  }
}

static void
hdf5_write_fields(struct psc_output_c *out, struct psc_fields_list *flds,
		  const char *pfx)
{
  void *ctx;
  hdf5_open(out, flds, pfx, &ctx);
  write_fields_combine(flds, __hdf5_write_field, ctx);
  hdf5_close(ctx);
}

// ======================================================================
// psc_output_format_ops_hdf5

struct psc_output_format_ops psc_output_format_ops_hdf5 = {
  .name         = "hdf5",
  .ext          = ".h5",
  .write_fields = hdf5_write_fields,
};

// ======================================================================

static void
xdmf_open(struct psc_output_c *out, struct psc_fields_list *list, const char *pfx,
	  void **pctx)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  char filename[strlen(out->data_dir) + 30];
  sprintf(filename, "%s/%s_%06d_%07d.h5", out->data_dir, pfx, rank, psc.timestep);
  
  __hdf5_open(out, list, filename, pctx);

  if (rank == 0) {
    char fname[strlen(out->data_dir) + strlen(pfx) + 20];
    sprintf(fname, "%s/%s_%07d.xdmf", out->data_dir, pfx, psc.timestep);
    FILE *f = fopen(fname, "w");

    fields_base_t *fld = &list->flds[0];
    fprintf(f, "<?xml version=\"1.0\" ?>\n");
    fprintf(f, "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n");
    fprintf(f, "<Domain>\n");
    fprintf(f, "<Grid GridType='Collection' CollectionType='Spatial'>\n");
    fprintf(f, "   <Time Type=\"Single\" Value=\"%g\" />\n", psc.timestep * psc.dt);
    int im[3];
    for (int d = 0; d < 3; d++) {
      im[d] = (psc.domain.ihi[d] - psc.domain.ilo[d]) / psc.domain.nproc[d];
    }
    for (int kz = 0; kz < psc.domain.nproc[2]; kz++) {
      for (int ky = 0; ky < psc.domain.nproc[1]; ky++) {
	for (int kx = 0; kx < psc.domain.nproc[0]; kx++) {
	  fprintf(f, "   <Grid Name=\"mesh-%d-%d-%d-%d\" GridType=\"Uniform\">\n",
		  kx, ky, kz, psc.timestep);
	  fprintf(f, "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"%d %d %d\"/>\n",
		  im[2] + 1, im[1] + 1, im[0] + 1);
	  fprintf(f, "     <Geometry GeometryType=\"Origin_DxDyDz\">\n");
	  fprintf(f, "     <DataStructure Name=\"Origin\" DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">\n");
	  fprintf(f, "        %g %g %g\n",
		  (psc.domain.ilo[2] + im[2] * kz) * psc.dx[2],
		  (psc.domain.ilo[1] + im[1] * ky) * psc.dx[1],
		  (psc.domain.ilo[0] + im[0] * kx) * psc.dx[0]);
	  fprintf(f, "     </DataStructure>\n");
	  fprintf(f, "     <DataStructure Name=\"Spacing\" DataType=\"Float\" Dimensions=\"3\" Format=\"XML\">\n");
	  fprintf(f, "        %g %g %g\n", psc.dx[2], psc.dx[1], psc.dx[0]);
	  fprintf(f, "     </DataStructure>\n");
	  fprintf(f, "     </Geometry>\n");
	  fprintf(f, "\n");
	  for (int m = 0; m < list->nr_flds; m++) {
	    fld = &list->flds[m];
	    fprintf(f, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
		    fld->name[0]);
	    fprintf(f, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
		    im[2], im[1], im[0]);
	    int proc = (kz * psc.domain.nproc[1] + ky) * psc.domain.nproc[0] + kx;
	    fprintf(f, "        %s_%06d_%07d.h5:/psc/fields/%s\n",
		    pfx, proc, psc.timestep, fld->name[0]);
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

static void
xdmf_close(void *ctx)
{
  __hdf5_close(ctx);
}

static void
xdmf_write_fields(struct psc_output_c *out, struct psc_fields_list *flds,
		  const char *pfx)
{
  void *ctx;
  xdmf_open(out, flds, pfx, &ctx);
  for (int m = 0; m < flds->nr_flds; m++) {
    __hdf5_write_field(ctx, &flds->flds[m]);
  }
  xdmf_close(ctx);
}

// ======================================================================
// psc_output_format_ops_xdmf

struct psc_output_format_ops psc_output_format_ops_xdmf = {
  .name         = "xdmf",
  .ext          = ".h5",
  .write_fields = xdmf_write_fields,
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
