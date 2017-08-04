
#include "psc_output_fields_private.h"

#include <mrc_io.h>

void
psc_output_fields_set_psc(struct psc_output_fields *output_fields, struct psc *psc)
{
  output_fields->psc = psc;
}

// ----------------------------------------------------------------------
// psc_output_fields_write

static void
_psc_output_fields_write(struct psc_output_fields *out, struct mrc_io *io)
{
  mrc_io_write_ref(io, out, "psc", out->psc);
}

// ----------------------------------------------------------------------
// psc_output_fields_read

static void
_psc_output_fields_read(struct psc_output_fields *out, struct mrc_io *io)
{
  out->psc = mrc_io_read_ref(io, out, "psc", psc);

  psc_output_fields_setup(out);
}

// ======================================================================
// forward to subclass

void
psc_output_fields_run(struct psc_output_fields *output_fields,
		      mfields_base_t *flds, struct psc_mparticles *particles)
{
  struct psc_output_fields_ops *ops = psc_output_fields_ops(output_fields);
  assert(ops->run);

  psc_stats_start(st_time_output);
  ops->run(output_fields, flds, particles);
  psc_stats_stop(st_time_output);
}

// ======================================================================
// psc_output_fields_init

static void
psc_output_fields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_output_fields, &psc_output_fields_c_ops);
}

// ======================================================================
// psc_output_fields class

struct mrc_class_psc_output_fields mrc_class_psc_output_fields = {
  .name             = "psc_output_fields",
  .size             = sizeof(struct psc_output_fields),
  .init             = psc_output_fields_init,
  .write            = _psc_output_fields_write,
  .read             = _psc_output_fields_read,
};

