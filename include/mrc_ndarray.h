
#ifndef MRC_NDARRAY_H
#define MRC_NDARRAY_H

#include <mrc_obj.h>

// ======================================================================
// mrc_ndarray, as one might guess from the name, is somewhat modeled after numpy.
// It's a n_dims array (up to MRC_FLD_MAXDIMS). The underlying data may not be
// contiguous, but it's accessed through a known stride per dimension

#define MRC_NDARRAY_MAXDIMS (5)

enum {
  MRC_NT_FLOAT,
  MRC_NT_DOUBLE,
  MRC_NT_INT,
  MRC_NT_NR,
};

// ----------------------------------------------------------------------
// struct mrc_ndarray_access
//
// This is separated out to allow for replication of the data pointer / strides
// needed to quickly access field values without having to follow a pointer to
// mrc_ndarray. So an mrc_obj like mrc_fld can use mrc_ndarray as a member object
// but maintain performance

struct mrc_ndarray_access {
  void *arr_off; //< same as the data pointer (arr), but shifted by
		 //precalculated offset for faster access
  int stride[MRC_NDARRAY_MAXDIMS];
};

#define __MRC_NDARRAY(nd_acc, type, i0,i1,i2,i3,i4)			\
  (((type *) (nd_acc)->arr_off)[(i4) * (nd_acc)->stride[4] +		\
				(i3) * (nd_acc)->stride[3] +		\
				(i2) * (nd_acc)->stride[2] +		\
				(i1) * (nd_acc)->stride[1] +		\
				(i0) * (nd_acc)->stride[0]])

#ifdef BOUNDS_CHECK

#include <string.h>

#define MRC_NDARRAY(fld, type, i0,i1,i2,i3,i4)				\
  (*({									\
      if (strcmp(#type, "float") == 0) assert(nd->data_type == MRC_NT_FLOAT); \
      if (strcmp(#type, "double") == 0) assert(nd->data_type == MRC_NT_DOUBLE); \
      if (strcmp(#type, "int") == 0) assert(nd->data_type == MRC_NT_INT); \
      assert(i0 >= (nd)->offs.vals[0] && i0 < (nd)->offs.vals[0] + (nd)->dims.vals[0]); \
      assert(i1 >= (nd)->offs.vals[1] && i1 < (nd)->offs.vals[1] + (nd)->dims.vals[1]); \
      assert(i2 >= (nd)->offs.vals[2] && i2 < (nd)->offs.vals[2] + (nd)->dims.vals[2]); \
      assert(i3 >= (nd)->offs.vals[3] && i3 < (nd)->offs.vals[3] + (nd)->dims.vals[3]); \
      assert(i4 >= (nd)->offs.vals[4] && i4 < (nd)->offs.vals[4] + (nd)->dims.vals[4]); \
      assert((nd)->_nd_acc.arr_off);					\
      type *_p = &__MRC_NDARRAY(&nd->acc, type, i0,i1,i2,i3,i4);			\
      _p; }))

#else

#define MRC_NDARRAY(nd, type, i0,i1,i2,i3,i4) __MRC_NDARRAY(&nd->acc, type, i0,i1,i2,i3,i4)

#endif

// ----------------------------------------------------------------------
// struct mrc_ndarray

struct mrc_ndarray {
  struct mrc_obj obj;
  
  // state
  struct mrc_ndarray_access acc;
  void *arr; //< pointer to the actual data
  int start[MRC_NDARRAY_MAXDIMS];
  int n_dims;
  int size_of_type;
  int data_type;
  int len; // number of data values in this ndarray
  struct mrc_vec *vec; //< underlying mrc_vec that manages memory alloc/free (could be petsc)

  // parameters
  struct mrc_param_int_array dims;
  struct mrc_param_int_array offs;
  struct mrc_param_int_array perm;

  // if view_base is set, this mrc_ndarray will be a view
  struct mrc_ndarray *view_base; 
  // the view will be of size offs/dims as usual above, but shifted to start at view_offs in the base nmrc_darray
  struct mrc_param_int_array view_offs;
};

struct mrc_ndarray_ops {
  MRC_SUBCLASS_OPS(struct mrc_ndarray);
};

struct mrc_ndarray_access *mrc_ndarray_access(struct mrc_ndarray *nd);
void mrc_ndarray_set_array(struct mrc_ndarray *nd, void *arr);
void mrc_ndarray_replace_array(struct mrc_ndarray *nd, void *arr);

extern struct mrc_ndarray_ops mrc_ndarray_float_ops;
extern struct mrc_ndarray_ops mrc_ndarray_double_ops;
extern struct mrc_ndarray_ops mrc_ndarray_int_ops;

MRC_CLASS_DECLARE(mrc_ndarray, struct mrc_ndarray);

// ======================================================================
// mrc_ndarray_it

#include <assert.h>

struct mrc_ndarray_it {
  float *ptr;
  int idx[MRC_NDARRAY_MAXDIMS];
  struct mrc_ndarray *nd;
};

static inline void
mrc_ndarray_it_start_all(struct mrc_ndarray_it *it, struct mrc_ndarray *nd)
{
  assert(nd->data_type == MRC_NT_FLOAT);
  
  it->nd = nd;
  for (int d = 0; d < nd->n_dims; d++) {
    it->idx[d] = nd->offs.vals[d];
  }
  it->ptr = (float *) nd->arr;
}

static inline bool
mrc_ndarray_it_done(struct mrc_ndarray_it *it)
{
  return !it->ptr;
}

static inline void
mrc_ndarray_it_next(struct mrc_ndarray_it *it)
{
  assert(!mrc_ndarray_it_done(it));
  
  for (int d = 0; d < it->nd->n_dims; d++) {
    it->idx[d]++;
    it->ptr += it->nd->acc.stride[d];
    if (it->idx[d] >= it->nd->offs.vals[d] + it->nd->dims.vals[d]) {
      it->idx[d] = it->nd->offs.vals[d];
      it->ptr -= it->nd->dims.vals[d] * it->nd->acc.stride[d];
    } else {
      goto done;
    }
  }

  it->ptr = NULL;

 done:
  ;
}

#endif



