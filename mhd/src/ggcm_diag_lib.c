
#include "ggcm_mhd_diag_private.h"

#include <mrc_io.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

static char *
make_basename(const char *run, float sheet, int outtype)
{
  int out_sheet;
  out_sheet=(sheet>=0.0) ? (int)(sheet*10.0001): -(int)(-sheet*10.0001);

  char *filename = malloc(strlen(run) + 20);

  switch (outtype) {
  case DIAG_TYPE_3D:
    sprintf(filename, "%s.3d", run);
    break;
  case DIAG_TYPE_2D_Z:
    sprintf(filename, "%s.pz_%d", run, out_sheet);
    break;
  case DIAG_TYPE_2D_X:
    sprintf(filename, "%s.px_%d", run, out_sheet);
    break;
  case DIAG_TYPE_2D_Y:
    sprintf(filename, "%s.py_%d", run, out_sheet);
    break;
  case DIAG_TYPE_2D_IONO:
    sprintf(filename, "%s.iof", run);
    break;
  default:
    assert(0);
  }

  return filename;
}

// ----------------------------------------------------------------------
// ggcm_diag_lib_create_mrc_io

struct mrc_io *
ggcm_diag_lib_create_mrc_io(MPI_Comm comm, const char *run, const char *outputmode,
			    int outtype, float sheet, int rank_diagsrv)
{
  struct mrc_io *io = mrc_io_create(comm);
  mrc_io_set_type(io, outputmode);
  // FIXME, this is a hack to allow us to set "--mrc_io_type xdmf2" for AMR,
  // but then "--mrc_io_iono_type xdmf_serial" to keep ionosphere output working
  if (outtype == DIAG_TYPE_2D_IONO) {
    mrc_io_set_name(io, "mrc_io_iono");
  }

  const char *outdir = ".";
  mrc_params_get_option_string("outdir", &outdir);
  char *basename = make_basename(run, sheet, outtype);
  mrc_io_set_param_string(io, "basename", basename);
  if (rank_diagsrv >= 0) {
    mrc_io_set_param_int(io, "rank_diagsrv", rank_diagsrv);
  }
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_view(io);
  return io;
}

// ----------------------------------------------------------------------
// ggcm_diag_lib_write_openggcm_attrs

void
ggcm_diag_lib_write_openggcm_attrs(struct mrc_io *io, const char *time_str)
{
  const char *run = "run";
  mrc_params_get_option_string("run", &run);

  struct mrc_obj *obj = mrc_obj_create(mrc_io_comm(io));
  mrc_obj_set_name(obj, "openggcm");
  mrc_obj_dict_add_string(obj, "run", run);
  mrc_obj_dict_add_string(obj, "time_str", time_str);
  mrc_obj_write(obj, io);
  mrc_obj_destroy(obj);
}



