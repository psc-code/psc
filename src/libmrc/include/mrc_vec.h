
#ifndef MRC_VEC_H
#define MRC_VEC_H

#include <mrc_common.h>
#include <mrc_obj.h>

// ======================================================================
// mrc_vec

MRC_CLASS_DECLARE(mrc_vec, struct mrc_vec);

// instead of allocating memory, use the pointer we provide
void mrc_vec_set_array(struct mrc_vec *vec, void *arr);

// gets the pointer to the allocated storage in the vector
void *mrc_vec_get_array(struct mrc_vec *vec);

// indicates that we're done using the data pointer obtained from
// mrc_vec_get_array
void mrc_vec_put_array(struct mrc_vec *vec, void *arr);

#endif
