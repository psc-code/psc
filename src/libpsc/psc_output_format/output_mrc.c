
#include "psc.h"
#include "psc_output_fields_c.h"
#include "psc_output_format_private.h"

#include <mrc_common.h>
#include <mrc_io.h>

#include <math.h>
#include <string.h>

enum {
  IO_TYPE_PFD,
  IO_TYPE_TFD,
  NR_IO_TYPES,
};

// FIXME, part of subctx
static struct mrc_io *ios[NR_IO_TYPES];

// ----------------------------------------------------------------------
// copy_to_mrc_fld

static void
copy_to_mrc_fld(struct mrc_m3 *m3, mfields_c_t *flds)
{
  psc_foreach_patch(ppsc, p) {
    fields_c_t *pf = psc_mfields_c_get_patch_c(flds, p);
    struct mrc_m3_patch *m3p = mrc_m3_patch_get(m3, p);
    mrc_m3_foreach(m3p, ix,iy,iz, 0,0) {
      MRC_M3(m3p,0, ix,iy,iz) = F3_C(pf,0, ix,iy,iz);
    } mrc_m3_foreach_end;
    mrc_m3_patch_put(m3);
  }
}

// ----------------------------------------------------------------------
// psc_output_format_mrc_destroy

static void
psc_output_format_mrc_destroy(struct psc_output_format *format)
{
  // FIXME, this makes mrc_io persistent across new objects from
  // this class (see that static var above), which is bad, but actually
  // helps getting a proper temporal XDMF file...
  return;
  
  for (int i = 0; i < NR_IO_TYPES; i++) {
    mrc_io_destroy(ios[i]);
    ios[i] = NULL;
  }
}

// ----------------------------------------------------------------------
// psc_output_format_mrc_write_fields

static void
psc_output_format_mrc_write_fields(struct psc_output_format *format,
				   struct psc_output_fields_c *out,
				   struct psc_fields_list *list,
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

  int gdims[3];
  mrc_domain_get_global_dims(ppsc->mrc_domain, gdims);
  int slab_off[3], slab_dims[3];
  for (int d = 0; d < 3; d++) {
    if (out->rx[d] > gdims[d])
      out->rx[d] = gdims[d];
    
    slab_off[d] = out->rn[d];
    slab_dims[d] = out->rx[d] - out->rn[d];
  }

  mrc_io_open(io, "w", ppsc->timestep, ppsc->timestep * ppsc->dt);
  for (int m = 0; m < list->nr_flds; m++) {
    mfields_c_t *flds = list->flds[m];
    fields_c_t *fld = psc_mfields_c_get_patch_c(flds, 0);
    assert(fld->nr_comp == 1);

    // FIXME, what if !(ibn[0] == ibn[1] == ibn[2])
    // FIXME, 3 doesn't work -- how about 0?
    struct mrc_m3 *mrc_fld = mrc_domain_m3_create(ppsc->mrc_domain);
    mrc_m3_set_name(mrc_fld, fld->name[0]);
    mrc_m3_set_param_int(mrc_fld, "sw", 2);
    mrc_m3_setup(mrc_fld);
    // FIXME, there should be a little function for this
    mrc_fld->name[0] = strdup(fld->name[0]);
    copy_to_mrc_fld(mrc_fld, flds);

    if (strcmp(mrc_io_type(io), "xdmf_collective") == 0) {
      mrc_io_set_param_int3(io, "slab_off", slab_off);
      mrc_io_set_param_int3(io, "slab_dims", slab_dims);
    }
    mrc_m3_write(mrc_fld, io);

    mrc_m3_destroy(mrc_fld);
  }
  mrc_io_close(io);
}

// ======================================================================
// psc_output_format: subclass "mrc"

struct psc_output_format_ops psc_output_format_mrc_ops = {
  .name                  = "mrc",
  .write_fields          = psc_output_format_mrc_write_fields,
  .destroy               = psc_output_format_mrc_destroy,
};


