
#define FTYPE_SINGLE          1
#define FTYPE_C               2
#define FTYPE_FORTRAN         3
#define FTYPE_CUDA            4

#if FTYPE == FTYPE_SINGLE

#define fields_FTYPE_real_t fields_single_real_t

#elif FTYPE == FTYPE_C

#define fields_FTYPE_real_t fields_c_real_t

#elif FTYPE == FTYPE_FORTRAN

#define fields_FTYPE_real_t fields_fortran_real_t

#elif FTYPE == FTYPE_CUDA

#define fields_FTYPE_real_t fields_cuda_real_t

#endif

// ----------------------------------------------------------------------
// fields_FYTPE_real_t

#if FTYPE == FTYPE_SINGLE || FTYPE == FTYPE_CUDA

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

#endif

// ----------------------------------------------------------------------
// F3

// Lower bounds and dims are intentionally not called ilg, ihg, img,
// to lessen confusion with psc.ilg, psc.ihg, etc.. These bounds may
// be psc.ilo, psc.ihi or psc.ilg, psc.ihg, or something yet different.

#if FTYPE == FTYPE_SINGLE

#define F3_OFF_S(pf, fldnr, jx,jy,jz)					\
  ((((((fldnr - (pf)->first_comp)					\
       * (pf)->im[2] + ((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))				\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#elif FTYPE == FTYPE_C

#define F3_OFF_C(pf, fldnr, jx,jy,jz)					\
  ((((((fldnr - (pf)->first_comp)					\
       * (pf)->im[2] + ((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))				\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#elif FTYPE == FTYPE_FORTRAN

#define F3_OFF_FORTRAN(pf, jx,jy,jz)			\
  (((((((jz)-(pf)->ib[2]))				\
      * (pf)->im[1] + ((jy)-(pf)->ib[1]))		\
     * (pf)->im[0] + ((jx)-(pf)->ib[0]))))

#endif

#ifndef BOUNDS_CHECK // ------------------------------

#if FTYPE == FTYPE_SINGLE

#define F3_S(pf, fldnr, jx,jy,jz)		\
  (((fields_single_real_t *) (pf)->data)[F3_OFF_S(pf, fldnr, jx,jy,jz)])

#elif FTYPE == FTYPE_C

#define F3_C(pf, fldnr, jx,jy,jz)		\
  (((fields_c_real_t *) (pf)->data)[F3_OFF_C(pf, fldnr, jx,jy,jz)])

#elif FTYPE == FTYPE_FORTRAN

#define F3_FORTRAN(pf, fldnr, jx,jy,jz)					\
  (((fields_fortran_real_t **) (pf)->data)[fldnr][F3_OFF_FORTRAN(pf, jx,jy,jz)])

#endif

#else // BOUNDS_CHECK ------------------------------

#if FTYPE == FTYPE_SINGLE

#define F3_S(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_S(pf, fldnr, jx,jy,jz);				\
      assert(fldnr >= (pf)->first_comp && fldnr < (pf)->first_comp + (pf)->nr_comp); \
      assert(jx >= (pf)->ib[0] && jx < (pf)->ib[0] + (pf)->im[0]);	\
      assert(jy >= (pf)->ib[1] && jy < (pf)->ib[1] + (pf)->im[1]);	\
      assert(jz >= (pf)->ib[2] && jz < (pf)->ib[2] + (pf)->im[2]);	\
      &(((fields_single_real_t *) (pf)->data)[off]);			\
    }))

#elif FTYPE == FTYPE_C

#define F3_C(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_C(pf, fldnr, jx,jy,jz);				\
      assert(fldnr >= (pf)->first_comp && fldnr < (pf)->first_comp + (pf)->nr_comp); \
      assert(jx >= (pf)->ib[0] && jx < (pf)->ib[0] + (pf)->im[0]);	\
      assert(jy >= (pf)->ib[1] && jy < (pf)->ib[1] + (pf)->im[1]);	\
      assert(jz >= (pf)->ib[2] && jz < (pf)->ib[2] + (pf)->im[2]);	\
      &(((fields_c_real_t *) (pf)->data)[off]);				\
    }))

#elif FTYPE == FTYPE_FORTRAN

#define F3_FORTRAN(pf, fldnr, jx,jy,jz)					\
  (*({int off = F3_OFF_FORTRAN(pf, jx,jy,jz);				\
      assert(fldnr >= 0 && fldnr < (pf)->nr_comp);			\
      assert(jx >= (pf)->ib[0] && jx < (pf)->ib[0] + (pf)->im[0]);	\
      assert(jy >= (pf)->ib[1] && jy < (pf)->ib[1] + (pf)->im[1]);	\
      assert(jz >= (pf)->ib[2] && jz < (pf)->ib[2] + (pf)->im[2]);	\
      &(((fields_fortran_real_t **) (pf)->data)[fldnr][off]);		\
    }))

#endif

#endif // BOUNDS_CHECK ------------------------------

// ----------------------------------------------------------------------

#undef fields_FTYPE_real_t




