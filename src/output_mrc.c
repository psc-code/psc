
#include "psc.h"
#include "psc_output_fields_c.h"

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
copy_to_mrc_fld(struct mrc_m3 *m3, mfields_base_t *flds)
{
  foreach_patch(p) {
    fields_base_t *pf = &flds->f[p];
    struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);
    mrc_m3_foreach(m3p, ix,iy,iz, 0,0) {
      MRC_M3(m3p,0, ix,iy,iz) = F3_BASE(pf,0, ix,iy,iz);
    } mrc_m3_foreach_end;
    mrc_m3_patch_put(m3);
  }
}

// ----------------------------------------------------------------------
// mrc_write_fields

static void
mrc_write_fields(struct psc_output_fields_c *out, struct psc_fields_list *list,
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
    struct mrc_m3 *mrc_fld = mrc_domain_m3_create(psc.mrc_domain);
    mrc_m3_set_name(mrc_fld, fld->name[0]);
    mrc_m3_set_param_int(mrc_fld, "sw", 2);
    mrc_m3_setup(mrc_fld);
    // FIXME, there should be a little function for this
    mrc_fld->name[0] = strdup(fld->name[0]);
    copy_to_mrc_fld(mrc_fld, flds);

    mrc_m3_write(mrc_fld, io);

    mrc_m3_destroy(mrc_fld);
  }
  mrc_io_close(io);
}

// ======================================================================
// psc_output_format_ops_mrc

struct psc_output_format_ops psc_output_format_ops_mrc = {
  .name         = "mrc",
  .write_fields = mrc_write_fields,
};



