
#include "cuda_mfields.h"

#include <cstdio>
#include <cassert>

// ----------------------------------------------------------------------
// cuda_mfields_create

struct cuda_mfields *
cuda_mfields_create()
{
  struct cuda_mfields *cmflds = 
    (struct cuda_mfields *) calloc(1, sizeof(*cmflds));

  return cmflds;
}

// ----------------------------------------------------------------------
// cuda_mfields_destroy

void
cuda_mfields_destroy(struct cuda_mfields *cmflds)
{
  free(cmflds);
}

