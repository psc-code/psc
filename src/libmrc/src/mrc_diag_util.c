
#include "mrc_diag_private.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

// ======================================================================
// diagc_make_filename

char *
diagc_make_filename(struct diag_format *format, const char *ext,
		    float sheet, int outtype, int step)
{
  struct diag_info *info = &format->diag_info;

  int out_sheet;
  out_sheet=(sheet>=0.0) ? (int)(sheet*10.0001): -(int)(-sheet*10.0001);

  char *filename = malloc(strlen(info->outdir) + strlen(info->run) + 30);

  int size;
  MPI_Comm_size(format->obj.comm, &size);
  char proc[10];
  if (size > 1) {
    sprintf(proc, "_p%06d", format->rank);
  } else {
    strcpy(proc, "");
  }

  switch (outtype) {
  case DIAG_TYPE_3D:
    sprintf(filename, "%s/%s.3df.%06d%s%s", info->outdir, info->run, step,
	    proc, ext);
    break;
  case DIAG_TYPE_2D_Z:
    sprintf(filename, "%s/%s.pz_%d.%06d%s%s", info->outdir, info->run, out_sheet,
	    step, proc, ext);
    break;
  case DIAG_TYPE_2D_X:
    sprintf(filename, "%s/%s.px_%d.%06d%s%s", info->outdir, info->run, out_sheet,
	    step, proc, ext);
    break;
  case DIAG_TYPE_2D_Y:
    sprintf(filename, "%s/%s.py_%d.%06d%s%s", info->outdir, info->run, out_sheet,
	    step, proc, ext);
    break;
  case DIAG_TYPE_2D_IONO:
    sprintf(filename, "%s/%s.iof.%06d%s%s", info->outdir, info->run, step,
	    proc, ext);
    break;
  default:
    assert(0);
  }

  return filename;
}

