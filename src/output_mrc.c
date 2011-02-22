
#include "psc.h"
#include "output_fields.h"

#include <mrc_common.h>
#include <mrc_io.h>

#include <math.h>

enum {
  IO_TYPE_PFD,
  IO_TYPE_TFD,
  NR_IO_TYPES,
};

static struct mrc_io *ios[NR_IO_TYPES];

// ----------------------------------------------------------------------
// copy_to_mrc_fld

static void
copy_to_mrc_fld(struct mrc_f3 *mrc_fld, fields_base_t *fld)
{
  int sw = mrc_fld->sw;
  mrc_f3_foreach(mrc_fld, ix,iy,iz, 2,2) {
    MRC_F3(mrc_fld,0, ix,iy,iz) = F3_BASE(fld,0, ix-sw,iy-sw,iz-sw);
  } mrc_f3_foreach_end;
}

// ----------------------------------------------------------------------
// mrc_write_fields

static void
mrc_write_fields(struct psc_output_c *out, struct psc_fields_list *list,
		  const char *pfx)
{
  int io_type;
  if (strcmp(pfx, "pfd") == 0) {
    io_type = IO_TYPE_PFD;
  } else if (strcmp(pfx, "tfd") == 0) {
    io_type = IO_TYPE_TFD;
  } else {
    assert(0);
  }
  struct mrc_io *io = ios[io_type];
  if (!io) {
    io = mrc_io_create(MPI_COMM_WORLD);
    mrc_io_set_param_string(io, "basename", pfx);
    mrc_io_set_from_options(io);
    mrc_io_setup(io);
    mrc_io_view(io);
    ios[io_type] = io;
  }

  mrc_io_open(io, "w", psc.timestep, psc.timestep * psc.dt);
  for (int m = 0; m < list->nr_flds; m++) {
    mfields_base_t *flds = &list->flds[m];
    fields_base_t *fld = &flds->f[0];
    assert(fld->nr_comp == 1);

    // FIXME, what if !(ibn[0] == ibn[1] == ibn[2])
    // FIXME, 3 doesn't work -- how about 0?
    struct mrc_f3 *mrc_fld = mrc_domain_f3_create(psc.mrc_domain, 2);
    mrc_f3_set_name(mrc_fld, fld->name[0]);
    mrc_f3_setup(mrc_fld);
    // FIXME, there should be a little function for this
    mrc_fld->name[0] = strdup(fld->name[0]);
    copy_to_mrc_fld(mrc_fld, fld);

    mrc_f3_write(mrc_fld, io);

    mrc_f3_destroy(mrc_fld);
  }
  mrc_io_close(io);
}

// ======================================================================
// psc_output_format_ops_mrc

struct psc_output_format_ops psc_output_format_ops_mrc = {
  .name         = "mrc",
  .write_fields = mrc_write_fields,
};



