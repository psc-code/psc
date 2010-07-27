
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

struct par_options_file {
  const char *options_file;
};

#define VAR(x) (void *)offsetof(struct par_options_file, x)

static struct param par_options_file_descr[] = {
  { "options_file"    , VAR(options_file)       , PARAM_STRING(NULL)  },
  {},
};

#undef VAR

static void
params_check_options_file(struct option **p_opt)
{
  struct par_options_file par;
  params_parse_cmdline(&par, par_options_file_descr, "options file", MPI_COMM_WORLD);
  params_print(&par, par_options_file_descr, "options file", MPI_COMM_WORLD);

  if (!par.options_file) {
    return;
  }

  FILE *file = fopen(par.options_file, "r");
  if (!file) {
    fprintf(stderr, "ERROR: cannot open options file '%s'!\n", par.options_file);
    abort();
  }
  
  while (!feof(file)) {
    char line[256];
    char *p = fgets(line, 256, file);
    if (!p)
      break;

    int l = strlen(line);
    if (l > 0 && line[l-1] == '\n')
      line[l-1] = 0;

    p = strchr(line, '#');
    if (p)
      *p = 0;

    l = strlen(line);
    if (!l)
      continue;

    char **ap, *argv[2], *s = line;
    int argc = 0;

    for (ap = argv; (*ap = strsep(&s, " \t")) != NULL;) {
      argc++;
      if (**ap != 0) {
	if (++ap >= &argv[2])
	  break;
      }
    }

    struct option *opt = malloc(sizeof(*opt));
    opt->name = strdup(argv[0] + 2);
    opt->value = NULL;
    opt->next = NULL;
    *p_opt = opt;
    p_opt = &opt->next;

    if (argc > 1) {
      opt->value = strdup(argv[1]);
    }

  }
  fclose(file);
}

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

  params_check_options_file(p);
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

static void
get_option_select(const char *name, struct param_select *descr, int *pval)
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
    case PT_SELECT:
      pv->u_select = params[i].u.ini_select;
      get_option_select(params[i].name, params[i].descr, &pv->u_select);
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
    case PT_SELECT:
      get_option_select(params[i].name, params[i].descr, &pv->u_select);
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
    case PT_SELECT:
      mpi_printf(comm, "%-20s| %s\n", params[i].name, params[i].descr[pv->u_select].str);
      break;
    }
  }
  mpi_printf(comm, "\n");
}

