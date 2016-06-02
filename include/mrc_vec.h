
#ifndef MRC_VEC_H
#define MRC_VEC_H

#include <mrc_common.h>
#include <mrc_obj.h>
#include <string.h>

// ======================================================================
// mrc_vec

MRC_CLASS_DECLARE(mrc_vec, struct mrc_vec);

// instead of allocating memory, use the pointer we provide
void mrc_vec_set_array(struct mrc_vec *vec, void *arr);

// replace our own allocated memory by the pointer provided (and free our mem)
void mrc_vec_replace_array(struct mrc_vec *vec, void *arr);

// gets the pointer to the allocated storage in the vector
void *mrc_vec_get_array(struct mrc_vec *vec);

// indicates that we're done using the data pointer obtained from
// mrc_vec_get_array
void mrc_vec_put_array(struct mrc_vec *vec, void *arr);


// These are also registered as methods for use in the
// time steppers.

// Vector math operations
void mrc_vec_axpy(struct mrc_vec *y, double alpha, struct mrc_vec *x);
void mrc_vec_waxpy(struct mrc_vec *w, double alpha, struct mrc_vec *x, struct mrc_vec *y);
void mrc_vec_axpby(struct mrc_vec *y, double alpha, struct mrc_vec *x, double beta);
// Data management operations
void mrc_vec_set(struct mrc_vec *x, double alpha);
void mrc_vec_copy(struct mrc_vec *vec_to, struct mrc_vec *vec_from);

int mrc_vec_len(struct mrc_vec *x);
int mrc_vec_size_of_type(struct mrc_vec *x);

#endif
