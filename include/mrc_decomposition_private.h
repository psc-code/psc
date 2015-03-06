
#ifndef MRC_DECOMPOSITION_PRIVATE_H
#define MRC_DECOMPOSITION_PRIVATE_H

#include "mrc_decomposition.h"

struct mrc_decomposition {
  struct mrc_obj obj;

  int n; // local vector size
  int N; // global vector size
  int off; // offset between local and global offsets (gidx = lidx + off)
};

#endif
