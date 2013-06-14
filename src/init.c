
#include <mrc.h>

unsigned long mrc_flags = 0;

void
mrc_set_flags(unsigned long flags)
{
  mrc_flags |= flags;
}

void
mrc_clear_flags(unsigned long flags)
{
  mrc_flags &= ~flags;
}

