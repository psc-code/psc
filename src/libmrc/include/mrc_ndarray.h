
#ifndef MRC_NDARRAY_H
#define MRC_NDARRAY_H

#include <mrc_obj.h>

BEGIN_C_DECLS

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
  size_t stride[MRC_NDARRAY_MAXDIMS];
  int beg[MRC_NDARRAY_MAXDIMS];
  int end[MRC_NDARRAY_MAXDIMS];
  int data_type;
};

#define ___MRC_NDARRAY(nd_acc, type, i0,i1,i2,i3,i4)			\
  (((type *) (nd_acc)->arr_off)[(i4) * (nd_acc)->stride[4] +		\
				(i3) * (nd_acc)->stride[3] +		\
				(i2) * (nd_acc)->stride[2] +		\
				(i1) * (nd_acc)->stride[1] +		\
				(i0) * (nd_acc)->stride[0]])

#ifdef BOUNDS_CHECK

#include <string.h>

#define __MRC_NDARRAY(nd_acc, type, i0,i1,i2,i3,i4)			\
  (*({									\
      if (strcmp(#type, "float") == 0) assert((nd_acc)->data_type == MRC_NT_FLOAT); \
      if (strcmp(#type, "double") == 0) assert((nd_acc)->data_type == MRC_NT_DOUBLE); \
      if (strcmp(#type, "int") == 0) assert((nd_acc)->data_type == MRC_NT_INT); \
      assert(i0 >= (nd_acc)->beg[0] && i0 < (nd_acc)->end[0]);		\
      assert(i1 >= (nd_acc)->beg[1] && i1 < (nd_acc)->end[1]);		\
      assert(i2 >= (nd_acc)->beg[2] && i2 < (nd_acc)->end[2]);		\
      assert(i3 >= (nd_acc)->beg[3] && i3 < (nd_acc)->end[3]);		\
      assert(i4 >= (nd_acc)->beg[4] && i4 < (nd_acc)->end[4]);		\
      assert((nd_acc)->arr_off);					\
      type *_p = &___MRC_NDARRAY(nd_acc, type, i0,i1,i2,i3,i4);		\
      _p; }))

#else

#define __MRC_NDARRAY(nd_acc, type, i0,i1,i2,i3,i4) ___MRC_NDARRAY(nd_acc, type, i0,i1,i2,i3,i4)

#endif

#define MRC_NDARRAY(nd, type, i0,i1,i2,i3,i4) __MRC_NDARRAY(&(nd)->nd_acc, type, i0,i1,i2,i3,i4)

#define MRC_S1(nd, i0) MRC_NDARRAY(nd, float, i0,0,0,0,0)
#define MRC_D1(nd, i0) MRC_NDARRAY(nd, double, i0,0,0,0,0)
#define MRC_I1(nd, i0) MRC_NDARRAY(nd, int, i0,0,0,0,0)

#define MRC_S2(nd, i0,i1) MRC_NDARRAY(nd, float, i0,i1,0,0,0)
#define MRC_D2(nd, i0,i1) MRC_NDARRAY(nd, double, i0,i1,0,0,0)
#define MRC_I2(nd, i0,i1) MRC_NDARRAY(nd, int, i0,i1,0,0,0)

#define MRC_S3(nd, i0,i1,i2) MRC_NDARRAY(nd, float, i0,i1,i2,0,0)
#define MRC_D3(nd, i0,i1,i2) MRC_NDARRAY(nd, double, i0,i1,i2,0,0)
#define MRC_I3(nd, i0,i1,i2) MRC_NDARRAY(nd, int, i0,i1,i2,0,0)

#define MRC_S4(nd, i0,i1,i2,i3) MRC_NDARRAY(nd, float, i0,i1,i2,i3,0)
#define MRC_D4(nd, i0,i1,i2,i3) MRC_NDARRAY(nd, double, i0,i1,i2,i3,0)
#define MRC_I4(nd, i0,i1,i2,i3) MRC_NDARRAY(nd, int, i0,i1,i2,i3,0)

#define MRC_S5(nd, i0,i1,i2,i3,i4) MRC_NDARRAY(nd, float, i0,i1,i2,i3,i4)
#define MRC_D5(nd, i0,i1,i2,i3,i4) MRC_NDARRAY(nd, double, i0,i1,i2,i3,i4)
#define MRC_I5(nd, i0,i1,i2,i3,i4) MRC_NDARRAY(nd, int, i0,i1,i2,i3,i4)

// ----------------------------------------------------------------------
// struct mrc_ndarray

struct mrc_ndarray {
  struct mrc_obj obj;
  
  // state
  struct mrc_ndarray_access nd_acc;
  void *arr; //< pointer to the actual data
  int start[MRC_NDARRAY_MAXDIMS];
  int n_dims;
  int size_of_type;
  int data_type;
  size_t len; // number of data values in this ndarray
  struct mrc_vec *vec; //< underlying mrc_vec that manages memory alloc/free (could be petsc)

  // parameters
  struct mrc_param_int_array dims;
  struct mrc_param_int_array offs;
  struct mrc_param_int_array perm;

  // if view_base is set, this mrc_ndarray will be a view
  struct mrc_ndarray *view_base; 
  // the view will be of size offs/dims as usual above, but shifted to start at view_offs in the base mrc_ndarray
  struct mrc_param_int_array view_offs;
};

struct mrc_ndarray_ops {
  MRC_SUBCLASS_OPS(struct mrc_ndarray);
  void (*set)(struct mrc_ndarray *nd, double val);
  void (*copy)(struct mrc_ndarray *to, struct mrc_ndarray *from);
  void (*scale)(struct mrc_ndarray *nd, double val);
  double (*norm)(struct mrc_ndarray *nd);
};

#define mrc_ndarray_ops(nd) ((struct mrc_ndarray_ops *) (nd)->obj.ops)

struct mrc_ndarray_access *mrc_ndarray_access(struct mrc_ndarray *nd);
void mrc_ndarray_set_array(struct mrc_ndarray *nd, void *arr);
void mrc_ndarray_replace_array(struct mrc_ndarray *nd, void *arr);

int mrc_ndarray_n_dims(struct mrc_ndarray *nd);
int *mrc_ndarray_dims(struct mrc_ndarray *nd);
int *mrc_ndarray_offs(struct mrc_ndarray *nd);
int mrc_ndarray_data_type(struct mrc_ndarray *nd);
bool mrc_ndarray_same_shape(struct mrc_ndarray *nd1, struct mrc_ndarray *nd2);
bool mrc_ndarray_f_contiguous(struct mrc_ndarray *nd);

void mrc_ndarray_set(struct mrc_ndarray *nd, double val);
void mrc_ndarray_copy(struct mrc_ndarray *to, struct mrc_ndarray *from);
void mrc_ndarray_scale(struct mrc_ndarray *nd, double val);
double mrc_ndarray_norm(struct mrc_ndarray *nd);

extern struct mrc_ndarray_ops mrc_ndarray_float_ops;
extern struct mrc_ndarray_ops mrc_ndarray_double_ops;
extern struct mrc_ndarray_ops mrc_ndarray_int_ops;

MRC_CLASS_DECLARE(mrc_ndarray, struct mrc_ndarray);

// ======================================================================
// mrc_ndarray_it

#include <assert.h>

struct mrc_ndarray_it {
  char *ptr;
  int idx[MRC_NDARRAY_MAXDIMS];
  int stride[MRC_NDARRAY_MAXDIMS]; // from nd, but in terms of bytes!
  int end[MRC_NDARRAY_MAXDIMS];
  int beg[MRC_NDARRAY_MAXDIMS];
  int n_dims;
};

#define IT_TYPE(it, type) (*(type *) (it)->ptr)
#define IT_S(it) IT_TYPE(it, float)
#define IT_D(it) IT_TYPE(it, double)
#define IT_I(it) IT_TYPE(it, int)

static inline void
mrc_ndarray_it_beg_end(struct mrc_ndarray_it *it, struct mrc_ndarray *nd,
		       int beg[], int end[])
{
  it->n_dims = nd->n_dims;
  it->ptr = (char *) nd->nd_acc.arr_off;
  for (int d = 0; d < it->n_dims; d++) {
    assert(beg[d] >= nd->offs.vals[d]);
    assert(end[d] <= nd->offs.vals[d] + nd->dims.vals[d]);
    it->end[d] = end[d];
    it->beg[d] = beg[d];
    it->idx[d] = beg[d];
    it->stride[d] = nd->nd_acc.stride[d] * nd->size_of_type;
    it->ptr += it->idx[d] * it->stride[d];
  }
}

static inline void
mrc_ndarray_it_all(struct mrc_ndarray_it *it, struct mrc_ndarray *nd)
{
  int beg[MRC_NDARRAY_MAXDIMS], end[MRC_NDARRAY_MAXDIMS];

  for (int d = 0; d < nd->n_dims; d++) {
    end[d] = nd->offs.vals[d] + nd->dims.vals[d];
    beg[d] = nd->offs.vals[d];
  }

  mrc_ndarray_it_beg_end(it, nd, beg, end);
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
  
  for (int d = 0; d < it->n_dims; d++) {
    it->idx[d]++;
    it->ptr += it->stride[d];
    if (it->idx[d] >= it->end[d]) {
      it->idx[d] = it->beg[d];
      it->ptr -= (it->end[d] - it->beg[d]) * it->stride[d];
    } else {
      goto done;
    }
  }

  it->ptr = NULL;

 done:
  ;
}

END_C_DECLS

#endif



