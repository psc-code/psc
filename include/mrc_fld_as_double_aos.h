
#ifndef MRC_FLD_AS_DOUBLE_AOS_H
#define MRC_FLD_AS_DOUBLE_AOS_H

typedef double mrc_fld_data_t;

#define F3(f, m, i,j,k) MRC_D4(f, m, i,j,k)
#define M3(f, m, i,j,k, p) MRC_D5(f, m, i,j,k, p)
#define FLD_TYPE "double_aos"
#define mrc_ddc_funcs_fld mrc_ddc_funcs_fld_double_aos

#endif
