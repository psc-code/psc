
#ifndef MRC_FLD_H
#define MRC_FLD_H

#define mrc_m3 mrc_fld

#include <mrc_common.h>
#include <mrc_obj.h>

#include <stdbool.h>
#include <assert.h>
#include <stdio.h>

//#define BOUNDS_CHECK

#define MRC_F3(f3,m, ix,iy,iz) MRC_S4(f3, ix,iy,iz,m)
#define MRC_F2(f3,m, ix,iy) MRC_S3(f3, ix,iy,m)
#define MRC_F1(f1,m, ix) MRC_S2(f1, ix,m)

struct mrc_io;

// ======================================================================
// mrc_fld

enum {
  MRC_NT_FLOAT,
  MRC_NT_DOUBLE,
  MRC_NT_INT,
  MRC_NT_NR,
};

#define MRC_FLD_MAXDIMS (5)

// for mrc_m3 emulation
struct mrc_fld_patch {
  int _p;
  struct mrc_fld *_fld;
};

struct mrc_fld {
  struct mrc_obj obj;
  // parameters
  struct mrc_param_int_array _dims;
  struct mrc_param_int_array _offs;
  struct mrc_param_int_array _sw;

  // state
  int _ghost_offs[MRC_FLD_MAXDIMS];
  int _ghost_dims[MRC_FLD_MAXDIMS];
  int _data_type;
  int _size_of_type;
  void *_arr;
  int _len;
  struct mrc_vec *_vec; //< underlying mrc_vec that manages memory alloc/free (could be petsc)
  struct mrc_domain *_domain; //< optional, if allocated through mrc_domain
  // for mrc_f3 emulation
  int _nr_allocated_comp_name;
  char **_comp_name;
  // for mrc_m3 emulation (FIXME, should be eliminated eventually (?))
  struct mrc_fld_patch *_patches;
  // for mrc_f1 emulation
  int _dim; //< # along this dim of the domain
};

MRC_CLASS_DECLARE(mrc_fld, struct mrc_fld);

void mrc_fld_set_array(struct mrc_fld *x, void *arr);
void mrc_fld_set_comp_names(struct mrc_fld *fld, const char *comps);
void mrc_fld_set_comp_name(struct mrc_fld *fld, int m, const char *name);
void mrc_fld_set_sw(struct mrc_fld *fld, int sw);
const char *mrc_fld_comp_name(struct mrc_fld *fld, int m);
void mrc_fld_set(struct mrc_fld *x, float val);
void mrc_fld_set_nr_comps(struct mrc_fld *fld, int nr_comps);
int mrc_fld_nr_comps(struct mrc_fld *fld);
const int *mrc_fld_offs(struct mrc_fld *x);
const int *mrc_fld_dims(struct mrc_fld *x);
const int *mrc_fld_ghost_offs(struct mrc_fld *x);
const int *mrc_fld_ghost_dims(struct mrc_fld *x);
struct mrc_fld *mrc_fld_duplicate(struct mrc_fld *fld);
void mrc_fld_copy(struct mrc_fld *fld_to, struct mrc_fld *fld_from);
void mrc_fld_write_comps(struct mrc_fld *fld, struct mrc_io *io, int mm[]);
void mrc_fld_dump(struct mrc_fld *fld, const char *basename, int n);
// for multi-patch mrc_fld only (former mrc_m3)
int mrc_fld_nr_patches(struct mrc_fld *fld);
float mrc_fld_norm_comp(struct mrc_fld *x, int m);

struct mrc_fld *mrc_fld_get_as(struct mrc_fld *fld_base,
			       const char *type);
void mrc_fld_put_as(struct mrc_fld *fld,
		    struct mrc_fld *fld_base);

static inline bool
mrc_fld_same_shape(struct mrc_fld *fld_1, struct mrc_fld *fld_2)
{
  if (fld_1->_dims.nr_vals != fld_2->_dims.nr_vals) {
    return false;
  }

  for (int d = 0; d < fld_1->_dims.nr_vals; d++) {
    if (fld_1->_ghost_dims[d] != fld_2->_ghost_dims[d]) {
      return false;
    }
  }
  return true;
}

#ifdef BOUNDS_CHECK

#include <string.h>

#define MRC_FLD(fld, type, i0,i1,i2,i3,i4)				\
  (*({									\
      if (strcmp(#type, "float") == 0) assert(fld->_data_type == MRC_NT_FLOAT);	\
      if (strcmp(#type, "double") == 0) assert(fld->_data_type == MRC_NT_DOUBLE); \
      if (strcmp(#type, "int") == 0) assert(fld->_data_type == MRC_NT_INT); \
      assert(i0 >= (fld)->_ghost_offs[0] && i0 < (fld)->_ghost_offs[0] + (fld)->_ghost_dims[0]); \
      assert(i1 >= (fld)->_ghost_offs[1] && i1 < (fld)->_ghost_offs[1] + (fld)->_ghost_dims[1]); \
      assert(i2 >= (fld)->_ghost_offs[2] && i2 < (fld)->_ghost_offs[2] + (fld)->_ghost_dims[2]); \
      assert(i3 >= (fld)->_ghost_offs[3] && i3 < (fld)->_ghost_offs[3] + (fld)->_ghost_dims[3]); \
      assert(i4 >= (fld)->_ghost_offs[4] && i4 < (fld)->_ghost_offs[4] + (fld)->_ghost_dims[4]); \
      assert((fld)->_arr);						\
      type *_p  =							\
	&(((type *) (fld)->_arr)[(((((i4) - (fld)->_ghost_offs[4]) *	\
				    (fld)->_ghost_dims[3] + (i3) - (fld)->_ghost_offs[3]) * \
				   (fld)->_ghost_dims[2] + (i2) - (fld)->_ghost_offs[2]) * \
				  (fld)->_ghost_dims[1] + (i1) - (fld)->_ghost_offs[1]) * \
				 (fld)->_ghost_dims[0] + (i0) - (fld)->_ghost_offs[0]]); \
      _p; }))

#else

#define MRC_FLD(fld, type, i0,i1,i2,i3,i4)				\
  (((type *) (fld)->_arr)[(((((i4) - (fld)->_ghost_offs[4]) *		\
			     (fld)->_ghost_dims[3] + (i3) - (fld)->_ghost_offs[3]) * \
			    (fld)->_ghost_dims[2] + (i2) - (fld)->_ghost_offs[2]) * \
			   (fld)->_ghost_dims[1] + (i1) - (fld)->_ghost_offs[1]) * \
			  (fld)->_ghost_dims[0] + (i0) - (fld)->_ghost_offs[0]])

#endif

#define MRC_S2(fld, i0,i1) MRC_FLD(fld, float, i0,i1,0,0,0)
#define MRC_D2(fld, i0,i1) MRC_FLD(fld, double, i0,i1,0,0,0)

#define MRC_S3(fld, i0,i1,i2) MRC_FLD(fld, float, i0,i1,i2,0,0)
#define MRC_D3(fld, i0,i1,i2) MRC_FLD(fld, double, i0,i1,i2,0,0)

#define MRC_S4(fld, i0,i1,i2,i3) MRC_FLD(fld, float, i0,i1,i2,i3,0)
#define MRC_D4(fld, i0,i1,i2,i3) MRC_FLD(fld, double, i0,i1,i2,i3,0)

#define MRC_S5(fld, i0,i1,i2,i3,i4) MRC_FLD(fld, float, i0,i1,i2,i3,i4)
#define MRC_D5(fld, i0,i1,i2,i3,i4) MRC_FLD(fld, double, i0,i1,i2,i3,i4)

#define mrc_fld_foreach(fld, ix,iy,iz, l,r)				\
  for (int iz = (fld)->_offs.vals[2] - (l); iz < (fld)->_offs.vals[2] + (fld)->_dims.vals[2] + (r); iz++) { \
  for (int iy = (fld)->_offs.vals[1] - (l); iy < (fld)->_offs.vals[1] + (fld)->_dims.vals[1] + (r); iy++) { \
  for (int ix = (fld)->_offs.vals[0] - (l); ix < (fld)->_offs.vals[0] + (fld)->_dims.vals[0] + (r); ix++) \

#define mrc_fld_foreach_end			\
  }						\
    } do {} while (0)				\

// ----------------------------------------------------------------------
// mrc_f1

#define mrc_f1_foreach(f1, ix, l,r)					\
  for (int ix = (f1)->_offs.vals[0] - l; ix < (f1)->_offs.vals[0] + (f1)->_dims.vals[0] + r; ix++) \

#define mrc_f1_foreach_end do {} while (0)	\

// ----------------------------------------------------------------------
// mrc_m3 emulation

#define mrc_fld_foreach_patch(m3, p) \
  for (int p = 0; p < mrc_fld_nr_patches(m3); p++)

static inline struct mrc_fld_patch *
mrc_fld_patch_get(struct mrc_fld *fld, int p)
{
  assert(fld->_patches);
  return &fld->_patches[p];
}

static inline void
mrc_fld_patch_put(struct mrc_fld *fld)
{
}

#define MRC_M3(m3p, m, ix,iy,iz) MRC_S5((m3p)->_fld, ix, iy, iz, m, (m3p)->_p)

#define mrc_m3_foreach(m3p, ix,iy,iz, l,r) {			\
  int _l[3] = { -l, -l, -l };					\
  int _r[3] = { m3p->_fld->_ghost_dims[0] + 2 * m3p->_fld->_ghost_offs[0] + r,	\
		m3p->_fld->_ghost_dims[1] + 2 * m3p->_fld->_ghost_offs[1] + r,	\
		m3p->_fld->_ghost_dims[2] + 2 * m3p->_fld->_ghost_offs[2] + r }; \
  for (int iz = _l[2]; iz < _r[2]; iz++) {			\
    for (int iy = _l[1]; iy < _r[1]; iy++) {			\
      for (int ix = _l[0]; ix < _r[0]; ix++)			\

#define mrc_m3_foreach_bnd(m3p, ix,iy,iz) {		\
  int _l[3] = { m3p->_fld->_ghost_offs[0], m3p->_fld->_ghost_offs[1], m3p->_fld->_ghost_offs[2] };	\
  int _r[3] = { m3p->_fld->_ghost_offs[0] + m3p->_fld->_ghost_dims[0],	\
		m3p->_fld->_ghost_offs[1] + m3p->_fld->_ghost_dims[1],	\
		m3p->_fld->_ghost_offs[2] + m3p->_fld->_ghost_dims[2] }; \
  for (int iz = _l[2]; iz < _r[2]; iz++) {				\
    for (int iy = _l[1]; iy < _r[1]; iy++) {				\
      for (int ix = _l[0]; ix < _r[0]; ix++)				\

#define mrc_m3_foreach_end  }}}
  
// ----------------------------------------------------------------------

struct mrc_fld_ops {
  MRC_SUBCLASS_OPS(struct mrc_fld);
  const char *vec_type;
};

typedef void (*mrc_fld_copy_to_func_t)(struct mrc_fld *,
				       struct mrc_fld *);
typedef void (*mrc_fld_copy_from_func_t)(struct mrc_fld *,
					 struct mrc_fld *);

// ======================================================================
// mrc_m1

struct mrc_m1_patch {
  float *arr;
  int im[1];
  int ib[1];
};

struct mrc_m1 {
  struct mrc_obj obj;
  int nr_comp;
  int nr_patches;
  int dim; //< this dimension of the domain
  struct mrc_m1_patch *patches;
  struct mrc_domain *domain;
  int sw;
  char **_comp_name;
};

MRC_CLASS_DECLARE(mrc_m1, struct mrc_m1);

void mrc_m1_set_comp_name(struct mrc_m1 *x, int m, const char *name);
const char *mrc_m1_comp_name(struct mrc_m1 *x, int m);
bool mrc_m1_same_shape(struct mrc_m1 *m1_1, struct mrc_m1 *m1_2);

static inline struct mrc_m1_patch *
mrc_m1_patch_get(struct mrc_m1 *m1, int p)
{
  return &m1->patches[p];
}

static inline void
mrc_m1_patch_put(struct mrc_m1 *m1)
{
}

#define MRC_M1(m1p,m, ix)					\
  ((m1p)->arr[((m) * (m1p)->im[0] + (ix) - (m1p)->ib[0])])

#define mrc_m1_foreach_patch(m1, p) \
  for (int p = 0; p < m1->nr_patches; p++)

#define mrc_m1_foreach(m1p, ix, l,r) {			\
  int _l[1] = { -l };					\
  int _r[1] = { m1p->im[0] + 2 * m1p->ib[0]  + r};	\
  for (int ix = _l[0]; ix < _r[0]; ix++)		\

#define mrc_m1_foreach_bnd(m1p, ix) {			\
  int _l[1] = { m1p->ib[0] };				\
  int _r[1] = { m1p->ib[0] + m1p->im[0]};		\
  for (int ix = _l[0]; ix < _r[0]; ix++)		\

#define mrc_m1_foreach_end  }
  

#endif
