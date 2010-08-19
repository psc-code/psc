
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
  
  H5LTmake_dataset_float(hdf5->group_fld, fld->name, 3, dims, fld->data);
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
