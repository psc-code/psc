
#ifndef MRC_VEC_PRIVATE_H
#define MRC_VEC_PRIVATE_H

#include "mrc_vec.h"

struct mrc_vec {
  struct mrc_obj obj;
  // parameters
  int len;
  bool with_array;

  // state
  int size_of_type;
  void *arr;
};

struct mrc_vec_ops {
  MRC_SUBCLASS_OPS(struct mrc_vec);
  void (*set_array)(struct mrc_vec *vec, void *arr);
  void (*replace_array)(struct mrc_vec *vec, void *arr);
  void *(*get_array)(struct mrc_vec *vec);
  void (*put_array)(struct mrc_vec *vec, void *arr);
  void (*axpy)(struct mrc_vec *y, double alpha, struct mrc_vec *x);
  void (*waxpy)(struct mrc_vec *w, double alpha, struct mrc_vec *x, struct mrc_vec *y);
  void (*axpby)(struct mrc_vec *y, double alpha, struct mrc_vec *x, double beta);
  void (*set)(struct mrc_vec *x, double val);
  void (*copy)(struct mrc_vec *vec_to, struct mrc_vec *vec_from);
};


extern struct mrc_vec_ops mrc_vec_petsc_ops;
#endif
