
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

#endif



