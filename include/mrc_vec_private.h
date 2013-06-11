
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
  void *(*get_array)(struct mrc_vec *vec);
  void (*put_array)(struct mrc_vec *vec, void *arr);
};


extern struct mrc_vec_ops mrc_vec_petsc_ops;
#endif
