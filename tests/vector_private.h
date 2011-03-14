
#ifndef VECTOR_PRIVATE_H
#define VECTOR_PRIVATE_H

#include "vector.h"

// ======================================================================
// vector class

struct vector {
  struct mrc_obj obj;
  int nr_elements;
};

// ======================================================================
// vector subclass

struct vector_ops {
  MRC_SUBCLASS_OPS(struct vector);
  void  (*set_element)(struct vector *vec, int i, double val);
  double (*get_element)(struct vector *vec, int i);
};

extern struct vector_ops vector_double_ops;
extern struct vector_ops vector_float_ops;

#define vector_ops(vec) ((struct vector_ops *)((vec)->obj.ops))

#endif
