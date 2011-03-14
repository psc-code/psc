
#ifndef VECTOR_H
#define VECTOR_H

#include <mrc_obj.h>

MRC_CLASS_DECLARE(vector, struct vector);

void vector_set_element(struct vector *vec, int i, double val);
double vector_get_element(struct vector *vec, int i);

// ======================================================================
// vector class

struct vector {
  struct mrc_obj obj;
  int nr_elements;
  double *elements;
};

#endif
