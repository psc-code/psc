
#include "psc.h"
#include "util/profile.h"
#include "util/params.h"
#include "output_fields.h"

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
hdf5_open(struct psc_fields_list *list, const char *filename, void **pctx)
{
  struct hdf5_ctx *hdf5 = malloc(sizeof(*hdf5));

  hdf5->file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  hdf5->group = H5Gcreate(hdf5->file, "psc", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_int(hdf5->group, ".", "timestep", &psc.timestep, 1);
  H5LTset_attribute_double(hdf5->group, ".", "dt", &psc.dt, 1);
  H5LTset_attribute_double(hdf5->group, ".", "dx", psc.dx, 3);
  H5LTset_attribute_int(hdf5->group, ".", "lo", psc.glo, 3);
  H5LTset_attribute_int(hdf5->group, ".", "hi", psc.ghi, 3);

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

struct psc_hdf5 {
  int field_next_out;
  int field_step;
};

#define VAR(x) (void *)offsetof(struct psc_hdf5, x)

static struct param psc_hdf5_descr[] = {
  { "field_first_out"    , VAR(field_next_out)       , PARAM_INT(0)        },
  { "field_step_out"     , VAR(field_step)           , PARAM_INT(10)        },
  {},
};

#undef VAR

static struct psc_hdf5 psc_hdf5;

static void hdf5_out_create(void)
{ 
  params_parse_cmdline(&psc_hdf5, psc_hdf5_descr, "PSC HDF5", MPI_COMM_WORLD);
  params_print(&psc_hdf5, psc_hdf5_descr, "PSC HDF5", MPI_COMM_WORLD);
};

static void
hdf5_out_field()
{
  if (psc.timestep < psc_hdf5.field_next_out) {
    return;
  }
  psc_hdf5.field_next_out += psc_hdf5.field_step;

  static int pr;
  if (!pr) {
    pr = prof_register("hdf5_out_field", 1., 0, 0);
  }
  prof_start(pr);

  struct psc_fields_list list = {};
  list.nr_flds = NR_FIELDS;
  for (int m = 0; m < NR_FIELDS; m++) {
    struct psc_field *fld = &list.flds[m];
    fld->name = fldname[m];
    fld->size = 1;
    for (int d = 0; d < 3; d++) {
      fld->ilo[d] = psc.ilo[d];
      fld->ihi[d] = psc.ihi[d];
      fld->size *= (psc.ihi[d] - psc.ilo[d]);
    }
    fld->data = calloc(fld->size, sizeof(*fld->data));
    int my = psc.ihi[1] - psc.ilo[1];
    int mx = psc.ihi[0] - psc.ilo[0];

    // FIXME, this is an extra copy which I'd rather avoid...
    for (int iz = psc.ilo[2]; iz < psc.ihi[2]; iz++) {
      for (int iy = psc.ilo[1]; iy < psc.ihi[1]; iy++) {
	for (int ix = psc.ilo[0]; ix < psc.ihi[0]; ix++) {
	  fld->data[((iz-psc.ilo[2]) * my + iy-psc.ilo[1]) * mx + ix-psc.ilo[0]] =
	    FF3(m, ix, iy, iz);
	}
      }
    }
  }
  write_fields_1proc(&psc_output_format_ops_hdf5, &list, "field");

  for (int m = 0; m < NR_FIELDS; m++) {
    free(list.flds[m].data);
  }

  prof_stop(pr);
}

struct psc_output_ops psc_output_ops_hdf5 = {
  .name      = "hdf5",
  .create    = hdf5_out_create,
  .out_field = hdf5_out_field,
};

// ======================================================================
// psc_output_format_ops_hdf5

struct psc_output_format_ops psc_output_format_ops_hdf5 = {
  .name         = "hdf5",
  .ext          = ".h5",
  .open         = hdf5_open,
  .close        = hdf5_close,
  .write_field  = hdf5_write_field,
};

