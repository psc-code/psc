
#include "mrc_io_private.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

// ======================================================================
// diagc_make_filename

char *
diagc_make_filename(struct mrc_io *io, const char *ext)
{
  struct mrc_io_params *par = &io->par;

  char *filename = malloc(strlen(par->outdir) + strlen(par->basename) + 30);
  sprintf(filename, "%s/%s.%06d%s", par->outdir, par->basename, io->step, ext);
  return filename;
}

