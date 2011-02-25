
#ifndef MRC_FLD_H
#define MRC_FLD_H

#include <mrc_common.h>
#include <mrc_obj.h>

#include <stdbool.h>

#define MRC_F1(f1,m, ix)					\
  ((f1)->arr[(m) * (f1)->im[0] + (ix) - (f1)->ib[0]])

#define MRC_F2(f2,m, ix,iy)					\
  ((f2)->arr[((m) * (f2)->im[1] + (iy) - (f2)->ib[1]) *		\
	      (f2)->im[0] + (ix) - (f2)->ib[0]])

#define MRC_F3(f3,m, ix,iy,iz)					\
  ((f3)->arr[(((m) * (f3)->im[2] + (iz) - (f3)->ib[2]) *	\
	      (f3)->im[1] + (iy) - (f3)->ib[1]) *		\
	     (f3)->im[0] + (ix) - (f3)->ib[0]])

struct mrc_io;

// ======================================================================
// mrc_f1

extern struct mrc_class mrc_class_mrc_f1;

struct mrc_f1 {
  struct mrc_obj obj;
  float *arr;
  int im[1];
  int ib[1];
  int nr_comp;
  int len;
  bool with_array;
  char **name;
};

MRC_OBJ_DEFINE_STANDARD_METHODS(mrc_f1, struct mrc_f1);

// ======================================================================

struct mrc_f2 {
  float *arr;
  int im[2];
  int ib[2];
  int nr_comp;
  int len;
  bool with_array;
  struct mrc_domain *domain; //< optional, if allocated through mrc_domain
  int sw; //< # of ghost points
  char **name;
};

void mrc_f2_alloc(struct mrc_f2 *f2, int ib[2], int im[2], int nr_comp);
void mrc_f2_alloc_with_array(struct mrc_f2 *f2, int ib[2], int im[2], int nr_comp, float *arr);
void mrc_f2_free(struct mrc_f2 *f2);

extern struct mrc_class mrc_class_mrc_f3;

struct mrc_f3 {
  struct mrc_obj obj;
  float *arr;
  int im[3];
  int ib[3];
  int nr_comp;
  int len;
  bool with_array;
  struct mrc_domain *domain; //< optional, if allocated through mrc_domain
  int sw; //< # of ghost points
  char **name;
};

MRC_OBJ_DEFINE_STANDARD_METHODS(mrc_f3, struct mrc_f3);
struct mrc_f3 *mrc_f3_alloc(MPI_Comm comm, int ib[3], int im[3]);
struct mrc_f3 *mrc_f3_duplicate(struct mrc_f3 *f3);
void mrc_f3_set_nr_comps(struct mrc_f3 *f3, int nr_comps);
void mrc_f3_set_array(struct mrc_f3 *f3, float *arr);
void mrc_f3_copy(struct mrc_f3 *f3_to, struct mrc_f3 *f3_from);
void mrc_f3_set(struct mrc_f3 *f3, float val);
void mrc_f3_write(struct mrc_f3 *f3, struct mrc_io *io);
void mrc_f3_write_scaled(struct mrc_f3 *f3, struct mrc_io *io, float scale);
void mrc_f3_write_comps(struct mrc_f3 *f3, struct mrc_io *io, int mm[]);

static inline bool
mrc_f3_same_shape(struct mrc_f3 *f3_1, struct mrc_f3 *f3_2)
{
  return (f3_1->nr_comp == f3_2->nr_comp &&
	  f3_1->im[0] == f3_2->im[0] &&
	  f3_1->im[1] == f3_2->im[1] &&
	  f3_1->im[2] == f3_2->im[2]);
}

#define mrc_f3_foreach(f3, ix,iy,iz, l,r)				        \
  for (int iz = (f3)->ib[2] + l; iz < (f3)->ib[2] + (f3)->im[2] - r; iz++) {    \
    for (int iy = (f3)->ib[1] + l; iy < (f3)->ib[1] + (f3)->im[1] - r; iy++) {	\
      for (int ix = (f3)->ib[0] + l; ix < (f3)->ib[0] + (f3)->im[0] - r; ix++)	\

#define mrc_f3_foreach_end			\
  }						\
    } do {} while (0)				\

#endif
