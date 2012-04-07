
#ifndef PSC_FIELD_SINGLE_H
#define PSC_FIELD_SINGLE_H

#include "psc.h"

typedef float fields_single_real_t;
#define MPI_FIELDS_SINGLE_REAL MPI_FLOAT

// Lower bounds and dims are intentionally not called ilg, ihg, img,
// to lessen confusion with psc.ilg, psc.ihg, etc.. These bounds may
// be psc.ilo, psc.ihi or psc.ilg, psc.ihg, or something yet different.

typedef struct {
  fields_single_real_t *flds;
  int ib[3], im[3]; //> lower bounds and length per direction
  int nr_comp; //> nr of components
  int first_comp; // first component
} fields_single_t;

#define F3_OFF_S(pf, fldnr, jx,jy,jz)					\
  ((((((fldnr - (pf)->first_comp)					\
       * (pf)->im[2] + ((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))				\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#ifndef BOUNDS_CHECK

#define F3_S(pf, fldnr, jx,jy,jz)		\
  ((pf)->flds[F3_OFF_S(pf, fldnr, jx,jy,jz)])

#else

#define F3_S(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_S(pf, fldnr, jx,jy,jz);				\
      assert(fldnr >= (pf)->first_comp && fldnr < (pf)->first_comp + (pf)->nr_comp); \
      assert(jx >= (pf)->ib[0] && jx < (pf)->ib[0] + (pf)->im[0]);	\
      assert(jy >= (pf)->ib[1] && jy < (pf)->ib[1] + (pf)->im[1]);	\
      assert(jz >= (pf)->ib[2] && jz < (pf)->ib[2] + (pf)->im[2]);	\
      &((pf)->flds[off]);						\
    }))

#endif

void fields_single_alloc(fields_single_t *pf, int ib[3], int ie[3],
			 int nr_comp, int first_comp);
void fields_single_free(fields_single_t *pf);

static inline unsigned int
fields_single_size(fields_single_t *pf)
{
  return pf->im[0] * pf->im[1] * pf->im[2];
}

#endif
