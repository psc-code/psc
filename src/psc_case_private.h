
#ifndef PSC_CASE_PRIVATE_H
#define PSC_CASE_PRIVATE_H

#include <psc_case.h>

struct _psc_case {
  struct mrc_obj obj;
  struct psc_case *Case;
};

// ======================================================================

#define _psc_case_ops(c) ((struct _psc_case_ops *)((c)->obj.ops))

#endif
