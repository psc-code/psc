
#ifndef VECTOR_PRIVATE_H
#define VECTOR_PRIVATE_H

#include "vector.h"

// ======================================================================
// vector class

struct vector {
  struct mrc_obj obj;
  int nr_elements;
  double *elements;
};

#endif
