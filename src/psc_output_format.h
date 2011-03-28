
#include <mrc_obj.h>

#include "psc.h"
#include <psc_output_fields_c.h>

MRC_CLASS_DECLARE(psc_output_format, struct psc_output_format);

void psc_output_format_write_fields(struct psc_output_format *format,
				    struct psc_output_fields_c *out,
				    struct psc_fields_list *list, const char *pfx);
