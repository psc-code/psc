
#include <mrc_common.h>

#include "fields3d.hxx"

#include <string.h>
#include <stdlib.h>
#include <assert.h>

#ifndef __cplusplus
#error need C++
#endif

BEGIN_C_DECLS

#define FTYPE_SINGLE          1
#define FTYPE_C               2
#define FTYPE_FORTRAN         3
#define FTYPE_CUDA            4
#define FTYPE_VPIC            5

#if FTYPE == FTYPE_SINGLE

#define fields_FTYPE_real_t fields_single_real_t
#define fields_FTYPE_t fields_single_t
#define fields_FTYPE_t_ctor fields_single_t_ctor
#define fields_FTYPE_t_dtor fields_single_t_dtor
#define fields_FTYPE_t_mflds fields_single_t_mflds
#define fields_FTYPE_t_size fields_single_t_size
#define fields_FTYPE_t_zero_range fields_single_t_zero_range

#elif FTYPE == FTYPE_C

#define fields_FTYPE_real_t fields_c_real_t
#define fields_FTYPE_t fields_c_t
#define fields_FTYPE_t_ctor fields_c_t_ctor
#define fields_FTYPE_t_dtor fields_c_t_dtor
#define fields_FTYPE_t_mflds fields_c_t_mflds
#define fields_FTYPE_t_size fields_c_t_size
#define fields_FTYPE_t_zero_range fields_c_t_zero_range

#elif FTYPE == FTYPE_FORTRAN

#define fields_FTYPE_real_t fields_fortran_real_t
#define fields_FTYPE_t fields_fortran_t
#define fields_FTYPE_t_ctor fields_fortran_t_ctor
#define fields_FTYPE_t_dtor fields_fortran_t_dtor
#define fields_FTYPE_t_mflds fields_fortran_t_mflds
#define fields_FTYPE_t_size fields_fortran_t_size
#define fields_FTYPE_t_zero_range fields_fortran_t_zero_range

#elif FTYPE == FTYPE_CUDA

#define fields_FTYPE_real_t fields_cuda_real_t

#elif FTYPE == FTYPE_VPIC

#define fields_FTYPE_real_t fields_vpic_real_t
#define fields_FTYPE_t fields_vpic_t

#endif

// ----------------------------------------------------------------------
// fields_FYTPE_real_t

#if FTYPE == FTYPE_SINGLE || FTYPE == FTYPE_CUDA || FTYPE == FTYPE_VPIC

typedef float fields_FTYPE_real_t;

#elif FTYPE == FTYPE_C || FTYPE == FTYPE_FORTRAN

typedef double fields_FTYPE_real_t;

#endif

// ----------------------------------------------------------------------
// MPI_FIELDS_FTYPE_REAL

#if FTYPE == FTYPE_SINGLE

#define MPI_FIELDS_SINGLE_REAL MPI_FLOAT

#elif FTYPE == FTYPE_C

#define MPI_FIELDS_C_REAL MPI_DOUBLE

#elif FTYPE == FTYPE_FORTRAN

#define MPI_FIELDS_FORTRAN_REAL MPI_DOUBLE

#elif FTYPE == FTYPE_CUDA

#define MPI_FIELDS_CUDA_REAL MPI_FLOAT

#elif FTYPE == FTYPE_VPIC

#define MPI_FIELDS_VPIC_REAL MPI_FLOAT

#endif

// ----------------------------------------------------------------------
// _F3

// Lower bounds and dims are intentionally not called ilg, ihg, img,
// to lessen confusion with psc.ilg, psc.ihg, etc.. These bounds may
// be psc.ilo, psc.ihi or psc.ilg, psc.ihg, or something yet different.

#define _F3_OFF(flds, m, i,j,k)						\
  (((((((m) - (flds).first_comp)					\
       * (flds).im[2] + ((k)-(flds).ib[2]))				\
      * (flds).im[1] + ((j)-(flds).ib[1]))				\
     * (flds).im[0] + ((i)-(flds).ib[0]))))

#if FTYPE == FTYPE_VPIC

#define _F3_VPIC_OFF(flds, m, i,j,k)					\
  ((((((k)-(flds).ib[2])) * (flds).im[1] +				\
     ((j)-(flds).ib[1])) * (flds).im[0] +				\
    ((i)-(flds).ib[0])) * (flds).nr_comp +				\
   m)

#endif

#ifndef BOUNDS_CHECK // ------------------------------

#if FTYPE == FTYPE_SINGLE

#define _F3_S(flds, m, i,j,k)			\
  ((flds).data[_F3_OFF(flds, m, i,j,k)])

#elif FTYPE == FTYPE_C

#define _F3_C(flds, m, i,j,k)			\
  ((flds).data[_F3_OFF(flds, m, i,j,k)])

#elif FTYPE == FTYPE_FORTRAN

#define _F3_FORTRAN(flds, m, i,j,k)		\
  ((flds).data[_F3_OFF(flds, m, i,j,k)])

#elif FTYPE == FTYPE_VPIC

#define _F3_VPIC(flds, m, i,j,k)		\
  ((flds).data[_F3_VPIC_OFF(flds, m, i,j,k)])

#endif

#else // BOUNDS_CHECK ------------------------------

#if FTYPE == FTYPE_SINGLE

#define _F3_S(flds, m, i,j,k)						\
  (*({assert(m >= 0 && m < (flds).nr_comp);				\
      assert(i >= (flds).ib[0] && i < (flds).ib[0] + (flds).im[0]);	\
      assert(j >= (flds).ib[1] && j < (flds).ib[1] + (flds).im[1]);	\
      assert(k >= (flds).ib[2] && k < (flds).ib[2] + (flds).im[2]);	\
      &((flds).data[_F3_OFF(flds, m, i,j,k)]);				\
    }))

#elif FTYPE == FTYPE_C

#define _F3_C(flds, m, i,j,k)						\
  (*({assert(m >= 0 && m < (flds).nr_comp);				\
      assert(i >= (flds).ib[0] && i < (flds).ib[0] + (flds).im[0]);	\
      assert(j >= (flds).ib[1] && j < (flds).ib[1] + (flds).im[1]);	\
      assert(k >= (flds).ib[2] && k < (flds).ib[2] + (flds).im[2]);	\
      &((flds).data[_F3_OFF(flds, m, i,j,k)]);				\
    }))

#elif FTYPE == FTYPE_FORTRAN

#define _F3_FORTRAN(flds, m, i,j,k)					\
  (*({assert(m >= 0 && m < (flds).nr_comp);				\
      assert(i >= (flds).ib[0] && i < (flds).ib[0] + (flds).im[0]);	\
      assert(j >= (flds).ib[1] && j < (flds).ib[1] + (flds).im[1]);	\
      assert(k >= (flds).ib[2] && k < (flds).ib[2] + (flds).im[2]);	\
      &((flds).data[_F3_OFF(flds, m, i,j,k)]);				\
    }))

#elif FTYPE == FTYPE_VPIC

#define _F3_VPIC(flds, m, i,j,k)					\
  (*({assert(m >= 0 && m < (flds).nr_comp);				\
      assert(i >= (flds).ib[0] && i < (flds).ib[0] + (flds).im[0]);	\
      assert(j >= (flds).ib[1] && j < (flds).ib[1] + (flds).im[1]);	\
      assert(k >= (flds).ib[2] && k < (flds).ib[2] + (flds).im[2]);	\
      &((flds).data[_F3_VPIC_OFF(flds, m, i,j,k)]);			\
    }))

#endif

#endif // BOUNDS_CHECK ------------------------------

#if FTYPE == FTYPE_SINGLE || FTYPE == FTYPE_C || FTYPE == FTYPE_FORTRAN || FTYPE == FTYPE_VPIC

// ======================================================================
// fields_FTYPE_t

struct fields_FTYPE_t : fields3d<fields_FTYPE_real_t>
{
};

// ----------------------------------------------------------------------
// fields_t_ctor

static inline fields_FTYPE_t
fields_FTYPE_t_ctor(int ib[3], int im[3], int n_comps)
{
  fields_FTYPE_t flds;

  unsigned int size = 1;
  for (int d = 0; d < 3; d++) {
    flds.ib[d] = ib[d];
    flds.im[d] = im[d];
    size *= im[d];
  }
  flds.nr_comp = n_comps;
  flds.first_comp = 0;
  flds.data = (fields_FTYPE_real_t *) calloc(size * flds.nr_comp, sizeof(*flds.data));

  return flds;
}

// ----------------------------------------------------------------------
// fields_t_dtor

static inline void
fields_FTYPE_t_dtor(fields_FTYPE_t *flds)
{
  free(flds->data);
  flds->data = NULL;
}

// ----------------------------------------------------------------------
// fields_t_mflds

#if FTYPE == FTYPE_SINGLE

struct psc_mfields;
fields_single_t psc_mfields_single_get_field_t(struct psc_mfields *mflds, int p);

static inline fields_single_t
fields_single_t_mflds(struct psc_mfields *mflds, int p)
{
  return psc_mfields_single_get_field_t(mflds, p);
}

#elif FTYPE == FTYPE_C

struct psc_mfields;
fields_c_t psc_mfields_c_get_field_t(struct psc_mfields *mflds, int p);

static inline fields_c_t
fields_c_t_mflds(struct psc_mfields *mflds, int p)
{
  return psc_mfields_c_get_field_t(mflds, p);
}

#elif FTYPE == FTYPE_VPIC

struct psc_mfields;
fields_vpic_t psc_mfields_vpic_get_field_t(struct psc_mfields *mflds, int p);

static inline fields_vpic_t
fields_vpic_t_mflds(struct psc_mfields *mflds, int p)
{
  return psc_mfields_vpic_get_field_t(mflds, p);
}

#elif FTYPE == FTYPE_FORTRAN

struct psc_mfields;
fields_fortran_t psc_mfields_fortran_get_field_t(struct psc_mfields *mflds, int p);

static inline fields_fortran_t
fields_fortran_t_mflds(struct psc_mfields *mflds, int p)
{
  return psc_mfields_fortran_get_field_t(mflds, p);
}

#endif

// ----------------------------------------------------------------------
// fields_t_size

static inline unsigned int
fields_FTYPE_t_size(fields_FTYPE_t flds)
{
  return flds.im[0] * flds.im[1] * flds.im[2];
}

// ----------------------------------------------------------------------
// fields_t_zero_range

static inline void
fields_FTYPE_t_zero_range(fields_FTYPE_t flds, int mb, int me)
{
  for (int m = mb; m < me; m++) {
#if FTYPE == FTYPE_SINGLE
    memset(&_F3_S(flds, m, flds.ib[0], flds.ib[1], flds.ib[2]), 0,
	   flds.im[0] * flds.im[1] * flds.im[2] * sizeof(fields_FTYPE_real_t));
#elif FTYPE == FTYPE_C
    memset(&_F3_C(flds, m, flds.ib[0], flds.ib[1], flds.ib[2]), 0,
	   flds.im[0] * flds.im[1] * flds.im[2] * sizeof(fields_FTYPE_real_t));
#else
    assert(0);
#endif
  }
}

#endif // FTYPE == FTYPE_SINGLE || FTYPE == FTYPE_C || FTYPE == FTYPE_FORTRAN

// ----------------------------------------------------------------------

#undef fields_FTYPE_real_t
#undef fields_FTYPE_t
#undef fields_FTYPE_t_ctor
#undef fields_FTYPE_t_dtor
#undef fields_FTYPE_t_mflds
#undef fields_FTYPE_t_size
#undef fields_FTYPE_t_zero_range

END_C_DECLS

