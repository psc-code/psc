
#ifndef MRC_FLD_H
#define MRC_FLD_H

#define mrc_m3 mrc_fld

#include <mrc_common.h>
#include <mrc_obj.h>
#include <mrc_ndarray.h>

#include <stdbool.h>
#include <assert.h>
#include <stdio.h>

//#define BOUNDS_CHECK

BEGIN_C_DECLS

#define MRC_F3(f3,m, ix,iy,iz) MRC_S4(f3, ix,iy,iz,m)
#define MRC_F2(f3,m, ix,iy) MRC_S3(f3, ix,iy,m)
#define MRC_F1(f1,m, ix) MRC_S2(f1, ix,m)

struct mrc_io;

// ======================================================================
// mrc_fld

// mrc_fld is a rather complicated beast. Let's try to break it down:
//
// It is used both as generic multi-d array (up to MRC_FLD_MAXDIMS), as well as a
// 3-d dimensional field living on a given domain (in parallel) (plus 1-d, too), with
// potentially a number of components, and multiple patches per MPI process.
// The local (per process) view of the latter is just a multi-d array, so it uses to
// generic infrastructure -- but the setup happens based on a
// mrc_domain, rather than specifying bounds manually.
//
// In its generic view, the field has up to MRC_FLD_MAXDIMS, each of
// which has a lower and an upper bound. I'll use python notation,
// e.g. bounds 2:5 means 2,3,4 (excl. of 5)
// These bounds are stored as ghost_offs, ghost_dims (where the upper
// bound is hence ghost_offs[d] + ghost_dims[d], FIXME?). To add to
// the confusion, there's also offs and dims -- which are identical to
// the ghost_* values if there are no ghost points guard cells / halos
// / whatever you like to call them. However, the _ghost_* values are really
// the more fundamental ones as far as the field is concerned, though
// in an application, it may be more convenient to think in terms of
// the "real" interior dimensions.
//
// In addition, there is now also _start and _stride, which kinda mirror
// _ghost_offs and _ghost_dims in a regularly allocated contiguous field.
// However, _start and _stride allow to create views of a given field,
// where the points contained inside the view are then not necessarily
// contiguous anymore, but rather just a view into an originally allocated
// contiguous field.

#define MRC_FLD_MAXDIMS MRC_NDARRAY_MAXDIMS

// for mrc_m3 emulation
struct mrc_fld_patch {
  int _p;
  struct mrc_fld *_fld;
};

// ----------------------------------------------------------------------
// struct mrc_fld

struct mrc_fld {
  struct mrc_obj obj;

  // state
  // these are copies from our ::nd member, replicated for fast access
  struct mrc_ndarray_access nd_acc;
  
  // parameters
  struct mrc_param_int_array _dims;
  struct mrc_param_int_array _offs;
  struct mrc_param_int_array _sw;
  // if (optinally) allocated through mrc_domain, _domain will be set, and the 
  // following quantities will be used to set up _dims, _offs, _sw in mrc_fld_setup()
  struct mrc_domain *_domain; //< optional, if allocated through mrc_domain
  int _nr_spatial_dims; //< number of spatial dims to use (1 or 3)
  int _dim; //< if number of spatial dims == 1, field is along this dim of the domain
  int _nr_comps; //< number of components of the field
  int _nr_ghosts; //< number of ghostpoints in non-invariant (dim > 1) directions

  // state
  struct mrc_ndarray *_nd;
  int _ghost_offs[MRC_FLD_MAXDIMS];
  int _ghost_dims[MRC_FLD_MAXDIMS];
  struct mrc_fld *_view_base; //< if this mrc_fld is a view, this is the field it's derived from
  int *_view_offs;
  int _nr_allocated_comp_name;
  char **_comp_name;
  bool _aos; //< indicates whether the layout (w.r.t to domain) is array-of-struct
  bool _c_order; //< indicates whether the layout is C (row-major) order (default false)
  // for mrc_m3 emulation (FIXME, should be eliminated eventually (?))
  struct mrc_fld_patch *_patches;
};

MRC_CLASS_DECLARE(mrc_fld, struct mrc_fld);

void mrc_fld_set_array(struct mrc_fld *x, void *arr);
void mrc_fld_replace_array(struct mrc_fld *x, void *arr);
void mrc_fld_set_comp_names(struct mrc_fld *fld, const char *comps);
void mrc_fld_set_comp_name(struct mrc_fld *fld, int m, const char *name);
const char *mrc_fld_comp_name(struct mrc_fld *fld, int m);
void mrc_fld_set(struct mrc_fld *x, float val);
int mrc_fld_nr_comps(struct mrc_fld *fld);
const int *mrc_fld_offs(struct mrc_fld *x);
const int *mrc_fld_dims(struct mrc_fld *x);
const int *mrc_fld_sw(struct mrc_fld *x);
const int *mrc_fld_ghost_offs(struct mrc_fld *x);
const int *mrc_fld_ghost_dims(struct mrc_fld *x);
int mrc_fld_data_type(struct mrc_fld *fld);
int mrc_fld_len(struct mrc_fld *fld);
struct mrc_fld *mrc_fld_duplicate(struct mrc_fld *fld);
struct mrc_fld *mrc_fld_make_view(struct mrc_fld *fld, int mb, int me);
void mrc_fld_copy(struct mrc_fld *fld_to, struct mrc_fld *fld_from);
void mrc_fld_axpy(struct mrc_fld *y, float alpha, struct mrc_fld *x);
void mrc_fld_axpby(struct mrc_fld *y, double alpha, struct mrc_fld *x, double beta);
float mrc_fld_norm(struct mrc_fld *fld);
void mrc_fld_write_comps(struct mrc_fld *fld, struct mrc_io *io, int mm[]);
void mrc_fld_dump(struct mrc_fld *fld, const char *basename, int n);
// for multi-patch mrc_fld only (former mrc_m3)
int mrc_fld_nr_patches(struct mrc_fld *fld);
float mrc_fld_norm_comp(struct mrc_fld *x, int m);
void mrc_fld_setup_vec(struct mrc_fld *fld);

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

static inline const int *
mrc_fld_spatial_dims(struct mrc_fld *x)
{
  return mrc_fld_dims(x);
}

static inline const int *
mrc_fld_spatial_offs(struct mrc_fld *x)
{
  return mrc_fld_offs(x);
}

static inline const int *
mrc_fld_spatial_sw(struct mrc_fld *x)
{
  return mrc_fld_sw(x);
}

#define MRC_FLD(fld, type, i0,i1,i2,i3,i4) MRC_NDARRAY(fld, type, i0,i1,i2,i3,i4)

#define mrc_fld_foreach(fld, ix,iy,iz, l,r) do {			\
  const int *_offs = mrc_fld_spatial_offs(fld);				\
  const int *_dims = mrc_fld_spatial_dims(fld);				\
  int _l[3] = { _offs[0]            - (_dims[0] > 1 ? (l) : 0),		\
		_offs[1]            - (_dims[1] > 1 ? (l) : 0),		\
		_offs[2]            - (_dims[2] > 1 ? (l) : 0) };	\
  int _r[3] = { _offs[0] + _dims[0] + (_dims[0] > 1 ? (r) : 0),		\
		_offs[1] + _dims[1] + (_dims[1] > 1 ? (r) : 0),		\
		_offs[2] + _dims[2] + (_dims[2] > 1 ? (r) : 0) };	\
  for (int iz = _l[2]; iz < _r[2]; iz++) {				\
  for (int iy = _l[1]; iy < _r[1]; iy++) {				\
  for (int ix = _l[0]; ix < _r[0]; ix++)				\

#define mrc_fld_foreach_end			\
  }						\
  }						\
  } while (0)				\

// FIXME? should use something like mrc_fld_spatial_ghost_offs/dims?
#define mrc_fld_foreach_bnd(fld, ix,iy,iz) do {				\
  const int *_offs = mrc_fld_ghost_offs(fld);				\
  const int *_dims = mrc_fld_ghost_dims(fld);				\
  int _l[3] = { _offs[0], _offs[1], _offs[2] };				\
  int _r[3] = { _offs[0] + _dims[0], _offs[1] + _dims[1], _offs[2] + _dims[2] }; \
  for (int iz = _l[2]; iz < _r[2]; iz++) {				\
  for (int iy = _l[1]; iy < _r[1]; iy++) {				\
  for (int ix = _l[0]; ix < _r[0]; ix++)				\

#define mrc_fld_foreach_bnd_end			\
  }						\
  }						\
  } while (0)				\

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
  void (*ddc_copy_to_buf)(struct mrc_fld *fld, int mb, int me, int p,
			  int ilo[3], int ihi[3], void *buf);
  void (*ddc_copy_from_buf)(struct mrc_fld *fld, int mb, int me, int p,
			    int ilo[3], int ihi[3], void *buf);
  void (*ddc_add_from_buf)(struct mrc_fld *fld, int mb, int me, int p,
			   int ilo[3], int ihi[3], void *buf);
};

// FIXME, fld should be the first argument (and struct mrc_fld *), but we want compatibility with the
// old-style mrc_ddc_funcs until the go away.
void mrc_fld_ddc_copy_to_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *buf,
			     void *fld);
void mrc_fld_ddc_copy_from_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *buf,
			       void *fld);
void mrc_fld_ddc_add_form_buf(int mb, int me, int p, int ilo[3], int ihi[3], void *buf,
			      void *fld);

typedef void (*mrc_fld_copy_to_func_t)(struct mrc_fld *,
				       struct mrc_fld *);
typedef void (*mrc_fld_copy_from_func_t)(struct mrc_fld *,
					 struct mrc_fld *);

// ======================================================================
// mrc_m1

#define MRC_M1(fld, m, ix, p) MRC_S3(fld, ix, m, p)

#define mrc_m1_foreach_patch(m1, p) \
  for (int p = 0; p < mrc_fld_nr_patches(m1); p++)

#define mrc_m1_foreach(m1, ix, l,r) {					\
  int _l[1] = { -l };							\
  int _r[1] = { m1->_ghost_dims[0] + 2 * m1->_ghost_offs[0]  + r};	\
  for (int ix = _l[0]; ix < _r[0]; ix++)				\

#define mrc_m1_foreach_bnd(m1, ix) {					\
  int _l[1] = { m1->_ghost_offs[0] };					\
  int _r[1] = { m1->_ghost_offs[0] + m1->_ghost_dims[0]};		\
  for (int ix = _l[0]; ix < _r[0]; ix++)				\

#define mrc_m1_foreach_end  }
  
END_C_DECLS

#endif
