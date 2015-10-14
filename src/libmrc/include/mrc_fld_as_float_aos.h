
#ifndef MRC_FLD_AS_FLOAT_AOS_H
#define MRC_FLD_AS_FLOAT_AOS_H

typedef float mrc_fld_data_t;
#define mrc_fld_abs fabsf
#define mrc_fld_max fmaxf

#define F3(f, m, i,j,k) MRC_S4(f, m, i,j,k)
#define M3(f, m, i,j,k, p) MRC_S5(f, m, i,j,k, p)
#define FLD_TYPE "float_aos"

#endif
