
#include "field_array.h"

#include <mrc_common.h>

void FieldArray::advanceB(double frac)
{
  ::advance_b(this, frac);
}

