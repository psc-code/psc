
#ifndef MRC_A3_H
#define MRC_A3_H

#include <mrc_obj.h>

// FIXME, mrc_m3 (and f3) functionality is a subset of mrc_a3

// ======================================================================
// mrc_a3

struct mrc_a3_patch {
  float *arr;
  int im[3];
  int ib[3];
};

struct mrc_a3 {
  struct mrc_obj obj;
  int nr_comp;
  int nr_patches;
  struct mrc_a3_patch *patches;
  struct mrc_domain *domain; //< based on this mrc_domain
  int sw; //< # of ghost points
  char **name;
};

MRC_CLASS_DECLARE(mrc_a3, struct mrc_a3);

void mrc_a3_set_comp_name(struct mrc_a3 *a3, int m, const char *name);
const char *mrc_a3_comp_name(struct mrc_a3 *a3, int m);

static inline struct mrc_a3_patch *
mrc_a3_patch_get(struct mrc_a3 *a3, int p)
{
  return &a3->patches[p];
}

static inline void
mrc_a3_patch_put(struct mrc_a3 *a3)
{
}

#ifdef BOUNDS_CHECK

#define MRC_A3(a3p,m, ix,iy,iz) (*({					\
  assert((ix) >= (a3p)->ib[0] && ix < (a3p)->ib[0] + (a3p)->im[0]);     \
  assert((iy) >= (a3p)->ib[1] && iy < (a3p)->ib[1] + (a3p)->im[1]);     \
  assert((iz) >= (a3p)->ib[2] && iz < (a3p)->ib[2] + (a3p)->im[2]);     \
  float *__p = &((a3p)->arr[(((m) * (a3p)->im[2] + (iz) - (a3p)->ib[2]) * \
			    (a3p)->im[1] + (iy) - (a3p)->ib[1]) *	\
			   (a3p)->im[0] + (ix) - (a3p)->ib[0]]);	\
  __p;									\
      }))

#else

#define MRC_A3(a3p,m, ix,iy,iz)					\
  ((a3p)->arr[(((m) * (a3p)->im[2] + (iz) - (a3p)->ib[2]) *	\
	       (a3p)->im[1] + (iy) - (a3p)->ib[1]) *		\
	      (a3p)->im[0] + (ix) - (a3p)->ib[0]])

#endif

#define mrc_a3_foreach_patch(a3, p) \
  for (int p = 0; p < a3->nr_patches; p++)

#define mrc_a3_foreach(a3p, ix,iy,iz, l,r) {		\
  int _l[3] = { -l, -l, -l };				\
  int _r[3] = { a3p->im[0] + 2 * a3p->ib[0] + r,	\
		a3p->im[1] + 2 * a3p->ib[1] + r,	\
		a3p->im[2] + 2 * a3p->ib[2] + r};	\
  for (int iz = _l[2]; iz < _r[2]; iz++) {		\
    for (int iy = _l[1]; iy < _r[1]; iy++) {		\
      for (int ix = _l[0]; ix < _r[0]; ix++)		\

#define mrc_a3_foreach_bnd(a3p, ix,iy,iz) {		\
  int _l[3] = { a3p->ib[0], a3p->ib[1], a3p->ib[2] };	\
  int _r[3] = { a3p->ib[0] + a3p->im[0],		\
		a3p->ib[1] + a3p->im[1],		\
		a3p->ib[2] + a3p->im[2] };		\
  for (int iz = _l[2]; iz < _r[2]; iz++) {		\
    for (int iy = _l[1]; iy < _r[1]; iy++) {		\
      for (int ix = _l[0]; ix < _r[0]; ix++)		\

#define mrc_a3_foreach_end  }}}
  

#endif
