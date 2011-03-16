
#ifndef MRC_PARAMS_H
#define MRC_PARAMS_H

// FIXME, mrc_ prefix and cleanup

#include <mpi.h>
#include <stdbool.h>
#include <stddef.h>

enum param_type {
  PT_INT,
  PT_BOOL,
  PT_FLOAT,
  PT_DOUBLE,
  PT_STRING,
  PT_SELECT,
  PT_INT3,
  PT_FLOAT3,
};

#define PARAM_INT(x)      PT_INT,    .u = { .ini_int = (x), }
#define PARAM_BOOL(x)     PT_BOOL,   .u = { .ini_bool = (x), }
#define PARAM_FLOAT(x)    PT_FLOAT,  .u = { .ini_float = (x), }
#define PARAM_DOUBLE(x)   PT_DOUBLE, .u = { .ini_double = (x), }
#define PARAM_STRING(x)   PT_STRING, .u = { .ini_string = (x), }
#define PARAM_SELECT(x,d) PT_SELECT, .u = { .ini_select = (x), }, .descr = d
#define PARAM_INT3(x,y,z) PT_INT3,   .u = { .ini_int3 = { (x), (y), (z) }, }
#define PARAM_FLOAT3(x,y,z) PT_FLOAT3, .u = { .ini_float3 = { (x), (y), (z) }, }

union param_u {
  int u_int;
  bool u_bool;
  float u_float;
  double u_double;
  const char *u_string;
  int u_select;
  int u_int3[3];
  float u_float3[3];
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
    int    ini_bool;
    float  ini_float;
    double ini_double;
    const char *ini_string;
    int    ini_select;
    int    ini_int3[3];
    float  ini_float3[3];
  } u;
  struct mrc_param_select *descr;
};

void libmrc_params_init(int argc, char **argv);

void mrc_params_print_all(MPI_Comm comm);
void mrc_params_insert_option(const char *name, const char *val);
int  mrc_params_get_option_string(const char *name, const char **pval);
int  mrc_params_get_option_int(const char *name, int *pval);

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

struct mrc_io;

void mrc_params_read(void *p, struct param *params, const char *title,
		     struct mrc_io *io);
void mrc_params_write(void *p, struct param *params, const char *title,
		      struct mrc_io *io);

#endif
