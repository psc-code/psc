
#include "psc_output_fields_private.h"

#include <mrc_io.h>

// ======================================================================
// forward to subclass

void
psc_output_fields_run(struct psc_output_fields *output_fields,
		      MfieldsBase& mflds, MparticlesBase& mprts)
{
  struct psc_output_fields_ops *ops = psc_output_fields_ops(output_fields);
  assert(ops->run);

  psc_stats_start(st_time_output);
  ops->run(output_fields, mflds, mprts);
  psc_stats_stop(st_time_output);
}

// ======================================================================
// psc_output_fields_init

extern struct psc_output_fields_ops psc_output_fields_c_ops;

static void
psc_output_fields_init()
{
  mrc_class_register_subclass(&mrc_class_psc_output_fields, &psc_output_fields_c_ops);
}

// ======================================================================
// psc_output_fields class

struct mrc_class_psc_output_fields_ : mrc_class_psc_output_fields {
  mrc_class_psc_output_fields_() {
    name             = "psc_output_fields";
    size             = sizeof(struct psc_output_fields);
    init             = psc_output_fields_init;
  }
} mrc_class_psc_output_fields;

