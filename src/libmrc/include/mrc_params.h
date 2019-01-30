
#ifndef MRC_PARAMS_H
#define MRC_PARAMS_H

// FIXME, mrc_ prefix and cleanup

#include <mpi.h>
#include <stdbool.h>
#include <stddef.h>

#include <mrc_common.h>

BEGIN_C_DECLS

enum param_type {
  PT_INT,   // 0
  PT_UINT,
  PT_BOOL,
  PT_FLOAT,
  PT_DOUBLE,
  PT_STRING,
  PT_SELECT,
  PT_INT3,
  PT_FLOAT3,
  PT_DOUBLE3,
  PT_PTR,   // 10
  PT_OBJ,
  PT_INT_ARRAY,
  PT_FLOAT_ARRAY,
  MRC_VAR_INT,
  MRC_VAR_BOOL,
  MRC_VAR_FLOAT,
  MRC_VAR_DOUBLE,
  MRC_VAR_OBJ,
  MRC_VAR_DOUBLE3,
};

#define PARAM_INT(x)      PT_INT,    .u = { .ini_int = (x), }
#define PARAM_UINT(x)     PT_UINT,   .u = { .ini_uint = (x), }
#define PARAM_BOOL(x)     PT_BOOL,   .u = { .ini_bool = (x), }
#define PARAM_FLOAT(x)    PT_FLOAT,  .u = { .ini_float = (x), }
#define PARAM_DOUBLE(x)   PT_DOUBLE, .u = { .ini_double = (x), }
#define PARAM_STRING(x)   PT_STRING, .u = { .ini_string = (x), }
#define PARAM_SELECT(x,d) PT_SELECT, .u = { .select = { .default_value = (x), .descr = (d), }, }
#define PARAM_INT3(x,y,z) PT_INT3,   .u = { .ini_int3 = { (x), (y), (z) }, }
#define PARAM_INT_ARRAY(n, x)   PT_INT_ARRAY,   .u = { .int_array = { .default_value = (x), .nr_vals = (n), }, }
#define PARAM_FLOAT_ARRAY(n, x) PT_FLOAT_ARRAY, .u = { .float_array = { .default_value = (x), .nr_vals = (n), }, }
#define PARAM_FLOAT3(x,y,z)     PT_FLOAT3,      .u = { .ini_float3 = { (x), (y), (z) }, }
#define PARAM_DOUBLE3(x,y,z)    PT_DOUBLE3,     .u = { .ini_double3 = { (x), (y), (z) }, }
#define PARAM_PTR(x)            PT_PTR,         .u = { .ini_ptr = (x), }
#define PARAM_OBJ(c)            PT_OBJ,         .u = { .mrc_obj = { .cls = (struct mrc_class *) &mrc_class_ ## c, }, }
#define MRC_VAR_OBJ(c)          MRC_VAR_OBJ,    .u = { .mrc_obj = { .cls = (struct mrc_class *) &mrc_class_ ## c, }, }

struct mrc_param_int_array {
  int nr_vals;
  int *vals;
};

struct mrc_param_float_array {
  int nr_vals;
  float *vals;
};

union param_u {
  int u_int;
  unsigned int u_uint;
  bool u_bool;
  float u_float;
  double u_double;
  const char *u_string;
  int u_select;
  int u_int3[3];
  float u_float3[3];
  double u_double3[3];
  void* u_ptr;
  struct mrc_param_int_array u_int_array;
  struct mrc_param_float_array u_float_array;
  struct mrc_obj *u_obj;
};

struct mrc_param_select {
  const char *str;
  int val;
};

struct param {
  const char *name;
  void *var;
  enum param_type type;
  union {
    int    ini_int;
    int    ini_uint;
    int    ini_bool;
    float  ini_float;
    double ini_double;
    const char *ini_string;
    struct { int default_value; struct mrc_param_select *descr; } select;
    int    ini_int3[3];
    float  ini_float3[3];
    float  ini_double3[3];
    void*  ini_ptr;
    struct { int default_value; int nr_vals; } int_array;
    struct { float default_value; int nr_vals; } float_array;
    struct { struct mrc_class *cls; } mrc_obj;
  } u;
  const char *help;
};

void libmrc_params_init(int argc, char **argv);
void libmrc_params_finalize(void);

void mrc_params_print_all(MPI_Comm comm);
void mrc_params_insert_option(const char *name, const char *val);
int  mrc_params_get_option_string(const char *name, const char **pval);
int  mrc_params_get_option_int(const char *name, int *pval);
int  mrc_params_get_option_uint(const char *name, unsigned int *pval);
int  mrc_params_get_option_float(const char *name, float *pval);
int  mrc_params_get_option_double(const char *name, double *pval);
int  mrc_params_get_option_bool(const char *name, bool *pval);
int  mrc_params_get_option_select(const char *name, struct mrc_param_select *descr, int *pval);
int  mrc_params_get_option_int_array(const char *name, struct mrc_param_int_array *pval);
int  mrc_params_get_option_float_array(const char *name, struct mrc_param_float_array *pval);
int  mrc_params_get_option_string_help(const char *name, const char **pval, const char *help);
int  mrc_params_get_option_int_help(const char *name, int *pval, const char *help);
int  mrc_params_get_option_uint_help(const char *name, unsigned int *pval, const char *help);
int  mrc_params_get_option_float_help(const char *name, float *pval, const char *help);
int  mrc_params_get_option_double_help(const char *name, double *pval, const char *help);
int  mrc_params_get_option_bool_help(const char *name, bool *pval, const char *help);
int  mrc_params_get_option_select_help(const char *name, struct mrc_param_select *descr, int *pval, const char *help);
int  mrc_params_get_option_int_array_help(const char *name, struct mrc_param_int_array *pval, const char *help);
int  mrc_params_get_option_float_array_help(const char *name, struct mrc_param_float_array *pval, const char *help);

// sets defaults from the descriptions
void mrc_params_set_default(void *p, struct param *params);

int  mrc_params_set_type(void *p, struct param *params, const char *name,
			 int type, union param_u *pval);
int  mrc_params_get_type(void *p, struct param *params, const char *name,
			 int type, union param_u *pval);
int  mrc_params_set_int(void *p, struct param *params, const char *name, int val);
int  mrc_params_set_string(void *p, struct param *params, const char *name, const char *val);

// parses the cmd line for the parameters described
// if an option is not provided, use the default value given in
// PARAM_DOUBLE(default_value)
void mrc_params_parse(void *p, struct param *params, const char *title,
		      MPI_Comm comm);

// parses the cmd line for the parameters described
// if an option is not provided, preserve the value unchanged
void mrc_params_parse_nodefault(void *p, struct param *params, const char *title,
				MPI_Comm comm);
// parses the cmd line for the parameters described, prefixed with "title_"
// if an option is not provided, preserve the value unchanged
void mrc_params_parse_pfx(void *p, struct param *params, const char *title,
			  MPI_Comm comm);
void mrc_params_print(void *p, struct param *params, const char *title,
		      MPI_Comm comm);
void mrc_params_print_one(void *p, struct param *prm, MPI_Comm comm);

// FIXME!!!
// This function should go away and users should convert to proper mrc_params infrastructure
int parse_float_array(const char *str, float *arr, int n);

END_C_DECLS

#endif
