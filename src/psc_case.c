
#include "psc_case_private.h"

// ======================================================================
// _psc_case_init

static void
_psc_case_init()
{
}

// ======================================================================
// _psc_case class

struct mrc_class__psc_case mrc_class__psc_case = {
  .name             = "_psc_case",
  .size             = sizeof(struct _psc_case),
  .init             = _psc_case_init,
};

