
#include "psc.h"

// ----------------------------------------------------------------------
// generic C data structures

struct c_particle {
  real xi, yi, zi;
  real pxi, pyi, pzi;
  real qni;
  real mni;
  real wni;
};

struct psc_genc {
  struct c_particle *part;
  float *flds;
};

void genc_push_part_yz_a();
void genc_push_part_yz_b();

#define F3_OFF(fldnr, jx,jy,jz)						\
  (((((fldnr								\
       *psc.img[2] + ((jz)-psc.ilg[2]))					\
      *psc.img[1] + ((jy)-psc.ilg[1]))					\
     *psc.img[0] + ((jx)-psc.ilg[0]))))

#if 1

#define F3(fldnr, jx,jy,jz)			\
  (genc->flds[F3_OFF(fldnr, jx,jy,jz)])

#else

#define F3(fldnr, jx,jy,jz)						\
  (*({int off = F3_OFF(fldnr, jx,jy,jz);				\
      assert(off >= 0);							\
      assert(off < 6*psc.fld_size);					\
      &(genc->flds[off]);						\
    }))

#endif

static inline int
nint(real x)
{
  return (int)(x + real(10.5)) - 10;
}

