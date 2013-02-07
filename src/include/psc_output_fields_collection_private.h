
#ifndef PSC_OUTPUT_FIELDS_COLLECTION_PRIVATE_H
#define PSC_OUTPUT_FIELDS_COLLECTION_PRIVATE_H

#include <psc_output_fields_collection.h>

struct psc_output_fields_collection {
  struct mrc_obj obj;
  struct psc *psc;
  char *names;
};

#endif
