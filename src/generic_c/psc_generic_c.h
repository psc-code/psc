
#include "psc.h"

// switch between double and float in generic_c
// constants need to always be given like 1.5f
// they will propagated to double precision as necessary
// when actually doing double computations

#if 1

typedef double creal;
#define creal_abs(x) fabs(x)
#define creal_sqrt(x) sqrt(x)

#else

typedef float creal;
#define creal_abs(x) fabsf(x)
#define creal_sqrt(x) sqrtf(x)
#endif

// ----------------------------------------------------------------------
// generic C data structures

struct c_particle {
  creal xi, yi, zi;
  creal pxi, pyi, pzi;
  creal qni;
  creal mni;
  creal wni;
};

struct psc_genc {
  struct c_particle *part;
  creal *flds;
};

void genc_push_part_xz();
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
nint(creal x)
{
  return (int)(x + 10.5f) - 10;
}

