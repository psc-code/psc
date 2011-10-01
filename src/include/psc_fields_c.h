
#ifndef PSC_FIELD_C_H
#define PSC_FIELD_C_H

#include "psc.h"

typedef double fields_c_real_t;
#define MPI_FIELDS_C_REAL MPI_DOUBLE

// Lower bounds and dims are intentionally not called ilg, ihg, img,
// to lessen confusion with psc.ilg, psc.ihg, etc.. These bounds may
// be psc.ilo, psc.ihi or psc.ilg, psc.ihg, or something yet different.

typedef struct {
  fields_c_real_t *flds;
  int ib[3], im[3]; //> lower bounds and length per direction
  int nr_comp; //> nr of components
  char **name; //> name for each component
  bool with_array; //> indicates whether data array was passed in instead of alloc'd
} fields_c_t;

#define F3_OFF_C(pf, fldnr, jx,jy,jz)					\
  ((((((fldnr)								\
       * (pf)->im[2] + ((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))				\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#ifndef BOUNDS_CHECK

#define F3_C(pf, fldnr, jx,jy,jz)		\
  ((pf)->flds[F3_OFF_C(pf, fldnr, jx,jy,jz)])

#else

#define F3_C(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_C(pf, fldnr, jx,jy,jz);				\
      assert(fldnr >= 0 && fldnr < (pf)->nr_comp);			\
      assert(jx >= (pf)->ib[0] && jx < (pf)->ib[0] + (pf)->im[0]);	\
      assert(jy >= (pf)->ib[1] && jy < (pf)->ib[1] + (pf)->im[1]);	\
      assert(jz >= (pf)->ib[2] && jz < (pf)->ib[2] + (pf)->im[2]);	\
      &((pf)->flds[off]);						\
    }))

#endif

void fields_c_alloc(fields_c_t *pf, int ib[3], int ie[3], int nr_comp);
void fields_c_alloc_with_array(fields_c_t *pf, int ib[3], int ie[3], int nr_comp,
			       fields_c_real_t *arr);
void fields_c_free(fields_c_t *pf);
void fields_c_zero(fields_c_t *pf, int m);
void fields_c_set(fields_c_t *pf, int m, fields_c_real_t val);
void fields_c_copy(fields_c_t *pf, int m_to, int m_from);
void fields_c_axpy_all(fields_c_t *y, fields_c_real_t a, fields_c_t *x);
void fields_c_scale_all(fields_c_t *pf, fields_c_real_t s);

static inline unsigned int
fields_c_size(fields_c_t *pf)
{
  return pf->im[0] * pf->im[1] * pf->im[2];
}

#endif
