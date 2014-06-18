
#ifndef MRC_FLD_AS_FLOAT_H
#define MRC_FLD_AS_FLOAT_H

typedef float mrc_fld_data_t;
#define mrc_fld_abs fabsf
#define mrc_fld_min fminf
#define mrc_fld_max fmaxf
#define mrc_fld_sqrt sqrtf

#define F3(f, m, i,j,k) MRC_S4(f, i,j,k, m)
#define M3(f, m, i,j,k, p) MRC_S5(f, i,j,k, m, p)
#define FLD_TYPE "float"
#define MPI_MRC_FLD_DATA_T MPI_FLOAT

#endif
