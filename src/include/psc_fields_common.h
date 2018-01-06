
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

#define fields_FTYPE_real_t _fields_single_real_t
#define fields_FTYPE_t fields_single_t
#define mfields_FTYPE_t mfields_single_t

#elif FTYPE == FTYPE_C

#define fields_FTYPE_real_t _fields_c_real_t
#define fields_FTYPE_t fields_c_t
#define mfields_FTYPE_t mfields_c_t

#elif FTYPE == FTYPE_FORTRAN

#define fields_FTYPE_real_t _fields_fortran_real_t
#define fields_FTYPE_t fields_fortran_t
#define mfields_FTYPE_t mfields_fortran_t

#elif FTYPE == FTYPE_CUDA

#define fields_FTYPE_real_t _fields_cuda_real_t
#define fields_FTYPE_t fields_cuda_t
#define mfields_FTYPE_t mfields_cuda_t

#elif FTYPE == FTYPE_VPIC

#define fields_FTYPE_real_t _fields_vpic_real_t
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

// ======================================================================
// fields_FTYPE_t

#if FTYPE == VPIC

struct fields_FTYPE_t : fields3d<fields_FTYPE_real_t, LayoutAOS>
{
  using mfields_t = mfields3d<fields_FTYPE_t>;
  
  using fields3d<fields_FTYPE_real_t>::fields3d;

};

#else

struct fields_FTYPE_t : fields3d<fields_FTYPE_real_t>
{
  using mfields_t = mfields3d<fields_FTYPE_t>;

  using fields3d<fields_FTYPE_real_t>::fields3d;

};

#endif // FTYPE == FTYPE_SINGLE || FTYPE == FTYPE_C || FTYPE == FTYPE_FORTRAN

// ----------------------------------------------------------------------

#undef fields_FTYPE_real_t
#undef fields_FTYPE_t
#undef mfields_FTYPE_t

