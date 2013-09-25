
#ifndef MRC_FLD_AS_DOUBLE_H
#define MRC_FLD_AS_DOUBLE_H

typedef double mrc_fld_data_t;

#define F3(f, m, i,j,k) MRC_D4(f, i,j,k, m)
#define M3(f, m, i,j,k, p) MRC_D5(f, i,j,k, m, p)
#define FLD_TYPE "double"
#define mrc_ddc_funcs_fld mrc_ddc_funcs_fld_double

#endif
