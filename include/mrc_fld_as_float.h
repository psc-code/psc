
#ifndef MRC_FLD_AS_FLOAT_H
#define MRC_FLD_AS_FLOAT_H

typedef float mrc_fld_data_t;

#define F3(f, m, i,j,k) MRC_S4(f, i,j,k, m)
#define M3(f, m, i,j,k, p) MRC_S5(f, i,j,k, m, p)
#define FLD_TYPE "float"
#define mrc_ddc_funcs_fld mrc_ddc_funcs_fld_float

#endif
