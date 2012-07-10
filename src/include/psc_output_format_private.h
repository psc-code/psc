
#ifndef PSC_OUTPUT_FORMAT_PRIVATE_H
#define PSC_OUTPUT_FORMAT_PRIVATE_H

#include <psc_output_format.h>

struct psc_output_format {
  struct mrc_obj obj;
};

struct psc_output_format_ops {
  MRC_SUBCLASS_OPS(struct psc_output_format);
};

// ======================================================================

#define psc_output_format_ops(format) ((struct psc_output_format_ops *)((format)->obj.ops))

#endif
