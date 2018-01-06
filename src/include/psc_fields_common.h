
#include <mrc_common.h>

#include "fields3d.hxx"

#include <string.h>
#include <stdlib.h>
#include <assert.h>

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
#define mfields_FTYPE_t mfields_single_t

#elif FTYPE == FTYPE_C

#define fields_FTYPE_real_t fields_c_real_t
#define fields_FTYPE_t fields_c_t
#define fields_FTYPE_t_ctor fields_c_t_ctor
#define fields_FTYPE_t_dtor fields_c_t_dtor
#define fields_FTYPE_t_mflds fields_c_t_mflds
#define fields_FTYPE_t_size fields_c_t_size
#define fields_FTYPE_t_zero_range fields_c_t_zero_range
#define mfields_FTYPE_t mfields_c_t

#elif FTYPE == FTYPE_FORTRAN

#define fields_FTYPE_real_t fields_fortran_real_t
#define fields_FTYPE_t fields_fortran_t
#define fields_FTYPE_t_ctor fields_fortran_t_ctor
#define fields_FTYPE_t_dtor fields_fortran_t_dtor
#define fields_FTYPE_t_mflds fields_fortran_t_mflds
#define fields_FTYPE_t_size fields_fortran_t_size
#define fields_FTYPE_t_zero_range fields_fortran_t_zero_range
#define mfields_FTYPE_t mfields_fortran_t

#elif FTYPE == FTYPE_CUDA

#define fields_FTYPE_real_t fields_cuda_real_t
#define fields_FTYPE_t fields_cuda_t
#define mfields_FTYPE_t mfields_cuda_t

#elif FTYPE == FTYPE_VPIC

#define fields_FTYPE_real_t fields_vpic_real_t
#define fields_FTYPE_t fields_vpic_t
#define mfields_FTYPE_t mfields_vpic_t

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

#if FTYPE == FTYPE_SINGLE || FTYPE == FTYPE_C || FTYPE == FTYPE_FORTRAN || FTYPE == FTYPE_VPIC

// ======================================================================
// fields_FTYPE_t

#if FTYPE == VPIC

struct fields_FTYPE_t : fields3d<fields_FTYPE_real_t, LayoutAOS>
{
  using mfields_t = mfields3d<fields_FTYPE_t>;
  
  using fields3d<fields_FTYPE_real_t>::fields3d;

  static fields_FTYPE_t psc_mfields_get_field_t(struct psc_mfields *mflds, int p)
  {
    fields_vpic_t psc_mfields_vpic_get_field_t(struct psc_mfields *mflds, int p);
    return psc_mfields_vpic_get_field_t(mflds, p);
  }
};

#else

struct fields_FTYPE_t : fields3d<fields_FTYPE_real_t>
{
  using mfields_t = mfields3d<fields_FTYPE_t>;

  using fields3d<fields_FTYPE_real_t>::fields3d;

#if FTYPE == FTYPE_SINGLE
  static fields_FTYPE_t psc_mfields_get_field_t(struct psc_mfields *mflds, int p)
  {
    fields_FTYPE_t psc_mfields_single_get_field_t(struct psc_mfields *mflds, int p);
    return psc_mfields_single_get_field_t(mflds, p);
  }
#elif FTYPE == FTYPE_C
  static fields_FTYPE_t psc_mfields_get_field_t(struct psc_mfields *mflds, int p)
  {
    fields_c_t psc_mfields_c_get_field_t(struct psc_mfields *mflds, int p);
    return psc_mfields_c_get_field_t(mflds, p);
  }
#elif FTYPE == FTYPE_FORTRAN
  static fields_FTYPE_t psc_mfields_get_field_t(struct psc_mfields *mflds, int p)
  {
    fields_FTYPE_t psc_mfields_fortran_get_field_t(struct psc_mfields *mflds, int p);
    return psc_mfields_fortran_get_field_t(mflds, p);
  }
#elif FTYPE == FTYPE_VPIC
  static fields_FTYPE_t psc_mfields_get_field_t(struct psc_mfields *mflds, int p)
  {
    fields_FTYPE_t psc_mfields_vpic_get_field_t(struct psc_mfields *mflds, int p);
    return psc_mfields_vpic_get_field_t(mflds, p);
  }
#endif
};

#endif

// ----------------------------------------------------------------------
// fields_t_ctor

static inline fields_FTYPE_t
fields_FTYPE_t_ctor(int ib[3], int im[3], int n_comps)
{
  return fields_FTYPE_t(ib, im, n_comps);
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
// mfields_t

using mfields_FTYPE_t = mfields3d<fields_FTYPE_t>;

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
  flds.zero(mb, me);
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
#undef mfields_FTYPE_t

