
#ifndef PARAMS_H
#define PARAMS_H

#include <mpi.h>
#include <stdbool.h>
#include <stddef.h>

enum param_type {
  PT_INT,
  PT_BOOL,
  PT_DOUBLE,
  PT_STRING,
  PT_SELECT
};

#define PARAM_INT(x)    PT_INT,    .u = { .ini_int = (x), }
#define PARAM_BOOL(x)   PT_BOOL,   .u = { .ini_bool = (x), }
#define PARAM_DOUBLE(x) PT_DOUBLE, .u = { .ini_double = (x), }
#define PARAM_STRING(x) PT_STRING, .u = { .ini_string = (x), }
#define PARAM_SELECT(x, d) PT_SELECT,  .u = { .ini_select = (x), }, .descr = d

union param_u {
  int u_int;
  bool u_bool;
  double u_double;
  const char *u_string;
  int u_select;
};

struct param_select {
  const char *str;
  int val;
};

struct param {
  const char *name;
  void *var;
  enum param_type type;
  union {
    int         ini_int;
    int         ini_bool;
    double      ini_double;
    const char *ini_string;
    int         ini_select;
  } u;
  struct param_select *descr;
};

void params_init(int argc, char **argv);
void params_print_all(void);

// parses the cmd line for the parameters described
// if an option is not provided, use the default value given in
// PARAM_DOUBLE(default_value)
void params_parse_cmdline(void *p, struct param *params, const char *title,
			  MPI_Comm comm);

// parses the cmd line for the parameters described
// if an option is not provided, preserve the value unchanged
void params_parse_cmdline_nodefault(void *p, struct param *params, const char *title,
				  MPI_Comm comm);
void params_print(void *p, struct param *params, const char *title,
		  MPI_Comm comm);

// hack that allows to print a message only on the first proc

#define mpi_printf(comm, args...) do { int __rank; MPI_Comm_rank(comm, &__rank); if (__rank == 0) { printf(args); } } while(0)

#endif
