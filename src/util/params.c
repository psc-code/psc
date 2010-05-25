
#include "params.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

struct option {
  const char *name;
  const char *value;
  struct option *next;
};

static struct option *option_list;

void
params_init(int argc, char **argv)
{
  struct option **p = &option_list;

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
    opt->next = NULL;
    *p = opt;
    p = &opt->next;

    // another arg left, which doesn't start with --
    if (i < argc - 1 && strncmp(argv[i+1], "--", 2) != 0) {
      opt->value = strdup(argv[i+1]);
      i++;
    }
  }
}

void
params_print_all()
{
  MPI_Comm comm = MPI_COMM_WORLD;

  mpi_printf(comm, "%-20s| %s\n", "parameter", "value");
  mpi_printf(comm, "--------------------+----------------------------------------\n");
  for (struct option *p = option_list; p; p = p->next) {
    mpi_printf(comm, "%-20s| %s\n", p->name, p->value);
  }
  mpi_printf(comm, "\n");
}

static struct option *
find_option(const char *name) {
  for (struct option *p = option_list; p; p = p->next) {
    if (strcmp(p->name, name) == 0)
      return p;
  }
  return NULL;
}

static void
get_option_int(const char *name, int *pval)
{
  struct option *p = find_option(name);
  
  if (!p)
    return;

  int rv = sscanf(p->value, "%d", pval);
  if (rv != 1) {
    fprintf(stderr, "error: cannot parse integer from '%s'\n", p->value);
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

static void
get_option_string(const char *name, const char **pval)
{
  struct option *p = find_option(name);
  
  if (!p)
    return;

  *pval = p->value;
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
params_parse_cmdline(void *p, struct param *params, const char *title,
		     MPI_Comm comm)
{
  for (int i = 0; params[i].name; i++) {
    union param_u *pv = p + (unsigned long) params[i].var;
    switch (params[i].type) {
    case PT_INT:
      pv->u_int = params[i].u.ini_int;
      get_option_int(params[i].name, &pv->u_int);
      break;
    case PT_BOOL:
      pv->u_bool = params[i].u.ini_bool;
      get_option_bool(params[i].name, &pv->u_bool);
      break;
    case PT_DOUBLE:
      pv->u_double = params[i].u.ini_double;
      get_option_double(params[i].name, &pv->u_double);
      break;
    case PT_STRING:
      pv->u_string = params[i].u.ini_string;
      get_option_string(params[i].name, &pv->u_string);
      break;
    }
  }
}

void
params_parse_cmdline_nodefault(void *p, struct param *params, const char *title,
			       MPI_Comm comm)
{
  for (int i = 0; params[i].name; i++) {
    union param_u *pv = p + (unsigned long) params[i].var;
    switch (params[i].type) {
    case PT_INT:
      get_option_int(params[i].name, &pv->u_int);
      break;
    case PT_BOOL:
      get_option_bool(params[i].name, &pv->u_bool);
      break;
    case PT_DOUBLE:
      get_option_double(params[i].name, &pv->u_double);
      break;
    case PT_STRING:
      get_option_string(params[i].name, &pv->u_string);
      break;
    }
  }
}

void
params_print(void *p, struct param *params, const char *title,
	     MPI_Comm comm)
{
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
    case PT_DOUBLE:
      mpi_printf(comm, "%-20s| %g\n", params[i].name, pv->u_double);
      break;
    case PT_STRING:
      mpi_printf(comm, "%-20s| %s\n", params[i].name, pv->u_string);
      break;
    }
  }
  mpi_printf(comm, "\n");
}

