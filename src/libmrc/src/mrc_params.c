
#include <mrc_common.h>
#include <mrc_params.h>
#include <mrc_list.h>
#include <mrc_io.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

struct option {
  const char *name;
  const char *value;
  list_t entry;
};

static LIST_HEAD(option_list);

void
libmrc_params_init(int argc, char **argv)
{
  for (int i = 1; i < argc; i++) {
    if (strncmp(argv[i], "--", 2) != 0) {
      fprintf(stderr,
	      "error: expected argument '%s' to specify an option like '--something'\n",
	      argv[i]);
      abort();
    }
    struct option *opt = malloc(sizeof(*opt));
    opt->name = strdup(argv[i] + 2);
    opt->value = NULL;

    // another arg left, which doesn't start with --
    if (i < argc - 1 && strncmp(argv[i+1], "--", 2) != 0) {
      opt->value = strdup(argv[i+1]);
      i++;
    }
    list_add_tail(&opt->entry, &option_list);
  }
}

void
mrc_params_insert_option(const char *name, const char *val)
{
  struct option *opt = malloc(sizeof(*opt));
  opt->name = strdup(name);
  opt->value = val ? strdup(val) : NULL;
  list_add_tail(&opt->entry, &option_list);
}

void
mrc_params_print_all(MPI_Comm comm)
{
  mpi_printf(comm, "%-20s| %s\n", "parameter", "value");
  mpi_printf(comm, "--------------------+----------------------------------------\n");
  struct option *p;
  list_for_each_entry(p, &option_list, entry) {
    mpi_printf(comm, "%-20s| %s\n", p->name, p->value);
  }
  mpi_printf(comm, "\n");
}

static struct option *
find_option(const char *name) {
  struct option *p;
  list_for_each_entry(p, &option_list, entry) {
    if (strcmp(p->name, name) == 0)
      return p;
  }
  return NULL;
}

int
mrc_params_get_option_int(const char *name, int *pval)
{
  struct option *p = find_option(name);
  
  if (!p)
    return -1;

  int rv = sscanf(p->value, "%d", pval);
  if (rv != 1) {
    fprintf(stderr, "error: cannot parse integer from '%s'\n", p->value);
    abort();
  }
  return 0;
}

static void
get_option_float(const char *name, float *pval)
{
  struct option *p = find_option(name);
  
  if (!p)
    return;

  int rv = sscanf(p->value, "%g", pval);
  if (rv != 1) {
    fprintf(stderr, "error: cannot parse float from '%s'\n", p->value);
    abort();
  }
}

static void
get_option_double(const char *name, double *pval)
{
  struct option *p = find_option(name);
  
  if (!p)
    return;

  int rv = sscanf(p->value, "%lg", pval);
  if (rv != 1) {
    fprintf(stderr, "error: cannot parse double from '%s'\n", p->value);
    abort();
  }
}

int
mrc_params_get_option_string(const char *name, const char **pval)
{
  struct option *p = find_option(name);
  
  if (!p)
    return -1;

  *pval = p->value;
  return 0;
}

static void
get_option_bool(const char *name, bool *pval)
{
  struct option *p = find_option(name);

  if (!p)
    return;

  if (!p->value) { // just "--something"
    *pval = true;
    return;
  }

  if (strcasecmp(p->value, "yes") == 0 ||
      strcasecmp(p->value, "true") == 0) {
    *pval = true;
  } else if (strcasecmp(p->value, "no") == 0 ||
	     strcasecmp(p->value, "false") == 0) {
    *pval = false;
  } else {
    fprintf(stderr, "error: cannot parse bool from '%s'\n", p->value);
    abort();
  }
}

void
mrc_params_get_option_select(const char *name, struct mrc_param_select *descr,
			     int *pval)
{
  struct option *p = find_option(name);
  if (!p) 
    return;

  if (!p->value) {
    fprintf(stderr, "ERROR: need to specify value for '%s'\n",
            name);
    abort();
  }

  for (int i = 0; descr[i].str; i++) {
    if (strcasecmp(descr[i].str, p->value) == 0) {
      *pval = i;
      return;
    }
  }

  fprintf(stderr, "ERROR: Select value '%s' not found. Valid options are:",
          p->value);
  for (int i = 0; descr[i].str; i++) {
    fprintf(stderr, " '%s'", descr[i].str);
  }
  fprintf(stderr, "\n");
  abort();
}

void
mrc_params_set_default(void *p, struct param *params)
{
  for (int i = 0; params[i].name; i++) {
    union param_u *pv = p + (unsigned long) params[i].var;
    switch (params[i].type) {
    case PT_INT:
      pv->u_int = params[i].u.ini_int;
      break;
    case PT_BOOL:
      pv->u_bool = params[i].u.ini_bool;
      break;
    case PT_FLOAT:
      pv->u_float = params[i].u.ini_float;
      break;
    case PT_DOUBLE:
      pv->u_double = params[i].u.ini_double;
      break;
    case PT_STRING:
      pv->u_string = params[i].u.ini_string;
      break;
    case PT_SELECT:
      pv->u_select = params[i].u.ini_select;
      break;
    }
  }
}

int
mrc_params_set_type(void *p, struct param *params, const char *name,
		    int type, union param_u *pval)
{
  for (int i = 0; params[i].name; i++) {
    if (strcmp(params[i].name, name) != 0)
      continue;

    union param_u *pv = p + (unsigned long) params[i].var;
    // types have to match, except a PT_SELECT can be set by a PT_INT
    if (params[i].type != type &&
	!(params[i].type == PT_SELECT && type == PT_INT)) {
      fprintf(stderr, "ERROR: option '%s' is not of type %d!\n",
	      name, type);
      abort();
    }
    switch (type) {
    case PT_INT:
      pv->u_int = pval->u_int;
      break;
    case PT_FLOAT:
      pv->u_float = pval->u_float;
      break;
    case PT_STRING:
      pv->u_string = pval->u_string;
      break;
    case PT_SELECT:
      pv->u_select = pval->u_select;
      break;
    default:
      assert(0);
    }
    return 0;
  }
  return -1; // not found
}

int
mrc_params_get_type(void *p, struct param *params, const char *name,
		    int type, union param_u *pval)
{
  for (int i = 0; params[i].name; i++) {
    if (strcmp(params[i].name, name) != 0)
      continue;

    union param_u *pv = p + (unsigned long) params[i].var;
    // types have to match, except a PT_SELECT can be set by a PT_INT
    if (params[i].type != type &&
	!(params[i].type == PT_SELECT && type == PT_INT)) {
      fprintf(stderr, "ERROR: option '%s' is not of type %d!\n",
	      name, type);
      abort();
    }
    switch (type) {
    case PT_INT:
      pval->u_int = pv->u_int;
      break;
    case PT_FLOAT:
      pval->u_float = pv->u_float;
      break;
    case PT_STRING:
      pval->u_string = pv->u_string;
      break;
    case PT_SELECT:
      pval->u_select = pv->u_select;
      break;
    default:
      assert(0);
    }
    return 0;
  }
  return -1; // not found
}

int
mrc_params_set_int(void *p, struct param *params, const char *name, int val)
{
  union param_u uval = { .u_int = val };
  return mrc_params_set_type(p, params, name, PT_INT, &uval);
}

int
mrc_params_set_string(void *p, struct param *params, const char *name, const char *val)
{
  union param_u uval = { .u_string = val };
  return mrc_params_set_type(p, params, name, PT_STRING, &uval);
}

void
mrc_params_parse(void *p, struct param *params, const char *title,
		 MPI_Comm comm)
{
  for (int i = 0; params[i].name; i++) {
    union param_u *pv = p + (unsigned long) params[i].var;
    switch (params[i].type) {
    case PT_INT:
      pv->u_int = params[i].u.ini_int;
      mrc_params_get_option_int(params[i].name, &pv->u_int);
      break;
    case PT_BOOL:
      pv->u_bool = params[i].u.ini_bool;
      get_option_bool(params[i].name, &pv->u_bool);
      break;
    case PT_FLOAT:
      pv->u_float = params[i].u.ini_float;
      get_option_float(params[i].name, &pv->u_float);
      break;
    case PT_DOUBLE:
      pv->u_double = params[i].u.ini_double;
      get_option_double(params[i].name, &pv->u_double);
      break;
    case PT_STRING:
      pv->u_string = params[i].u.ini_string;
      mrc_params_get_option_string(params[i].name, &pv->u_string);
      break;
    case PT_SELECT:
      pv->u_select = params[i].u.ini_select;
      mrc_params_get_option_select(params[i].name, params[i].descr, &pv->u_select);
      break;
    }
  }
}

void
mrc_params_parse_nodefault(void *p, struct param *params, const char *title,
			   MPI_Comm comm)
{
  for (int i = 0; params[i].name; i++) {
    union param_u *pv = p + (unsigned long) params[i].var;
    switch (params[i].type) {
    case PT_INT:
      mrc_params_get_option_int(params[i].name, &pv->u_int);
      break;
    case PT_BOOL:
      get_option_bool(params[i].name, &pv->u_bool);
      break;
    case PT_FLOAT:
      get_option_float(params[i].name, &pv->u_float);
      break;
    case PT_DOUBLE:
      get_option_double(params[i].name, &pv->u_double);
      break;
    case PT_STRING:
      mrc_params_get_option_string(params[i].name, &pv->u_string);
      break;
    case PT_SELECT:
      mrc_params_get_option_select(params[i].name, params[i].descr, &pv->u_select);
      break;
    }
  }
}

void
mrc_params_print(void *p, struct param *params, const char *title, MPI_Comm comm)
{
  mpi_printf(comm, "\n");
  mpi_printf(comm, "%-20s| %s\n", "parameter", "value");
  mpi_printf(comm, "--------------------+---------------------------------------- %s\n", title);
  for (int i = 0; params[i].name; i++) {
    union param_u *pv = p + (unsigned long) params[i].var;
    switch (params[i].type) {
    case PT_INT:
      mpi_printf(comm, "%-20s| %d\n", params[i].name, pv->u_int);
      break;
    case PT_BOOL:
      mpi_printf(comm, "%-20s| %s\n", params[i].name, pv->u_bool ? "yes" : "no");
      break;
    case PT_FLOAT:
      mpi_printf(comm, "%-20s| %g\n", params[i].name, pv->u_float);
      break;
    case PT_DOUBLE:
      mpi_printf(comm, "%-20s| %g\n", params[i].name, pv->u_double);
      break;
    case PT_STRING:
      mpi_printf(comm, "%-20s| %s\n", params[i].name, pv->u_string);
      break;
    case PT_SELECT:
      mpi_printf(comm, "%-20s| %s\n", params[i].name,
		 params[i].descr[pv->u_select].str);
      break;
    }
  }
}

void
mrc_params_write(void *p, struct param *params, const char *title, struct mrc_io *io)
{
  for (int i = 0; params[i].name; i++) {
    union param_u *pv = p + (unsigned long) params[i].var;
    mrc_io_write_attr(io, title, params[i].type, params[i].name, pv);
  }
}

void
mrc_params_read(void *p, struct param *params, const char *title, struct mrc_io *io)
{
  for (int i = 0; params[i].name; i++) {
    union param_u *pv = p + (unsigned long) params[i].var;
    mrc_io_read_attr(io, title, params[i].type, params[i].name, pv);
  }
}

