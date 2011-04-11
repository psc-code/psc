
#ifndef PSC_FIELD_FORTRAN_H
#define PSC_FIELD_FORTRAN_H

#include "psc.h"

typedef double fields_fortran_real_t;
#define MPI_FIELDS_FORTRAN_REAL MPI_DOUBLE

typedef struct {
  fields_fortran_real_t **flds;
  int ib[3], im[3]; //> lower bounds and length per direction
  int nr_comp; //> nr of components
  char **name; //> name for each component
  bool with_array; // array was passed in, not alloc'ed
} fields_fortran_t;

typedef struct {
  fields_fortran_t *f;
  int nr_patches;
  list_t entry;
} mfields_fortran_t;

#define F3_OFF_FORTRAN(pf, jx,jy,jz)			\
  (((((((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))		\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#ifndef BOUNDS_CHECK

#define F3_FORTRAN(pf, fldnr, jx,jy,jz)                \
  ((pf)->flds[fldnr][F3_OFF_FORTRAN(pf, jx,jy,jz)])

#else

#define F3_FORTRAN(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_FORTRAN(pf, jx,jy,jz);				\
      assert(fldnr >= 0 && fldnr < (pf)->nr_comp);			\
      assert(jx >= (pf)->ib[0] && jx < (pf)->ib[0] + (pf)->im[0]);	\
      assert(jy >= (pf)->ib[1] && jy < (pf)->ib[1] + (pf)->im[1]);	\
      assert(jz >= (pf)->ib[2] && jz < (pf)->ib[2] + (pf)->im[2]);	\
      &((pf)->flds[fldnr][off]);					\
    }))

#endif

void fields_fortran_alloc(fields_fortran_t *pf, int ib[3], int ie[3], int nr_comp);
void fields_fortran_alloc_with_array(fields_fortran_t *pf, int ib[3], int ie[3],
				     int nr_comp, fields_fortran_real_t *arr);
void fields_fortran_free(fields_fortran_t *pf);
void fields_fortran_get(mfields_fortran_t *pf, int mb, int me, void *flds_base);
void fields_fortran_get_from(mfields_fortran_t *pf, int mb, int me, void *flds_base,
			     int mb_base);
void fields_fortran_put(mfields_fortran_t *pf, int mb, int me, void *flds_base);
void fields_fortran_put_to(mfields_fortran_t *pf, int mb, int me, void *flds_base,
			   int mb_base);
void fields_fortran_zero(fields_fortran_t *pf, int m);
void fields_fortran_zero_all(fields_fortran_t *pf);
void fields_fortran_zero(fields_fortran_t *pf, int m);
void fields_fortran_set(fields_fortran_t *pf, int m, fields_fortran_real_t val);
void fields_fortran_copy(fields_fortran_t *pf, int m_to, int m_from);
void fields_fortran_axpy_all(fields_fortran_t *y, fields_fortran_real_t a,
			     fields_fortran_t *x);
void fields_fortran_scale_all(fields_fortran_t *pf, fields_fortran_real_t s);

static inline unsigned int
fields_fortran_size(fields_fortran_t *pf)
{
  return pf->im[0] * pf->im[1] * pf->im[2];
}

#endif
