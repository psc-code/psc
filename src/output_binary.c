
#include "psc.h"
#include "output_fields.h"
#include "util/profile.h"
#include "util/params.h"

#include <mpi.h>
#include <string.h>

struct psc_output_binary {
  char *data_dir;
  bool dowrite_pfield, dowrite_tfield;
  int pfield_next, tfield_next;
  int pfield_step, tfield_step;
  bool dowrite_fd[NR_EXTRA_FIELDS];
};

#define VAR(x) (void *)offsetof(struct psc_output_binary, x)

static struct param psc_binary_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(NULL)   },
  { "write_pfield"       , VAR(dowrite_pfield)       , PARAM_BOOL(1)        },
  { "pfield_first"       , VAR(pfield_next)          , PARAM_INT(0)         },
  { "pfield_step"        , VAR(pfield_step)          , PARAM_INT(10)        },
  { "write_tfield"       , VAR(dowrite_tfield)       , PARAM_BOOL(1)        },
  { "tfield_first"       , VAR(tfield_next)          , PARAM_INT(0)         },
  { "tfield_step"        , VAR(tfield_step)          , PARAM_INT(10)        },
  { "output_write_ex"    , VAR(dowrite_fd[X_EX])     , PARAM_BOOL(1)        },
  { "output_write_ey"    , VAR(dowrite_fd[X_EY])     , PARAM_BOOL(1)        },
  { "output_write_ez"    , VAR(dowrite_fd[X_EZ])     , PARAM_BOOL(1)        },
  { "output_write_bx"    , VAR(dowrite_fd[X_BX])     , PARAM_BOOL(1)        },
  { "output_write_by"    , VAR(dowrite_fd[X_BY])     , PARAM_BOOL(1)        },
  { "output_write_bz"    , VAR(dowrite_fd[X_BZ])     , PARAM_BOOL(1)        },
  { "output_write_jx"    , VAR(dowrite_fd[X_JXI])    , PARAM_BOOL(1)        },
  { "output_write_jy"    , VAR(dowrite_fd[X_JYI])    , PARAM_BOOL(1)        },
  { "output_write_jz"    , VAR(dowrite_fd[X_JZI])    , PARAM_BOOL(1)        },
  { "output_write_jxex"  , VAR(dowrite_fd[X_JXEX])   , PARAM_BOOL(1)        },
  { "output_write_jyey"  , VAR(dowrite_fd[X_JYEY])   , PARAM_BOOL(1)        },
  { "output_write_jzez"  , VAR(dowrite_fd[X_JZEZ])   , PARAM_BOOL(1)        },
  { "output_write_poyx"  , VAR(dowrite_fd[X_POYX])   , PARAM_BOOL(1)        },
  { "output_write_poyy"  , VAR(dowrite_fd[X_POYY])   , PARAM_BOOL(1)        },
  { "output_write_poyz"  , VAR(dowrite_fd[X_POYZ])   , PARAM_BOOL(1)        },
  { "output_write_e2x"   , VAR(dowrite_fd[X_E2X])    , PARAM_BOOL(1)        },
  { "output_write_e2y"   , VAR(dowrite_fd[X_E2Y])    , PARAM_BOOL(1)        },
  { "output_write_e2z"   , VAR(dowrite_fd[X_E2Z])    , PARAM_BOOL(1)        },
  { "output_write_b2x"   , VAR(dowrite_fd[X_B2X])    , PARAM_BOOL(1)        },
  { "output_write_b2y"   , VAR(dowrite_fd[X_B2Y])    , PARAM_BOOL(1)        },
  { "output_write_b2z"   , VAR(dowrite_fd[X_B2Z])    , PARAM_BOOL(1)        },
  {},
};

#undef VAR

static struct psc_output_binary psc_output_binary;

static void output_binary_create(void)
{ 
  params_parse_cmdline(&psc_output_binary, psc_binary_descr, "PSC BINARY", MPI_COMM_WORLD);
  params_print(&psc_output_binary, psc_binary_descr, "PSC BINARY", MPI_COMM_WORLD);
};


static void
binary_field_output_aux(struct psc_extra_fields *f, char *fnamehead)
{ 
  char *headstr = "PSC ";
  char *datastr = "DATA";

  // appears as "?BL?" if NO byte swapping required, ?LB? if required
  unsigned int magic_big_little = 1061962303;    
  unsigned int output_version = 1;
  
  float t_float;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char filename[200];
  sprintf(filename, "data/%s_%03d_%07d.psc", fnamehead, rank, psc.timestep);

  FILE *file = fopen(filename, "wb");

  // Header  
  fwrite(headstr, sizeof(char), 4, file);
  fwrite(&magic_big_little, sizeof(unsigned int), 1, file);
  fwrite(&output_version, sizeof(unsigned int), 1, file);

  t_float = (float) psc.dx[0];  fwrite(&t_float, sizeof(float), 1, file);
  t_float = (float) psc.dx[1];  fwrite(&t_float, sizeof(float), 1, file);
  t_float = (float) psc.dx[2];  fwrite(&t_float, sizeof(float), 1, file);
  t_float = (float) psc.dt;     fwrite(&t_float, sizeof(float), 1, file);

  // Indices on local proc
  fwrite(&psc.ilo[0], sizeof(psc.ilo[0]), 1, file);
  fwrite(&psc.ihi[0], sizeof(psc.ihi[0]), 1, file);
  fwrite(&psc.ilo[1], sizeof(psc.ilo[1]), 1, file);
  fwrite(&psc.ihi[1], sizeof(psc.ihi[1]), 1, file);
  fwrite(&psc.ilo[2], sizeof(psc.ilo[2]), 1, file);
  fwrite(&psc.ihi[2], sizeof(psc.ihi[2]), 1, file);

  // Globally saved indices (everything for now...)
  fwrite(&psc.domain.ilo[0], sizeof(psc.domain.ilo[0]), 1, file);
  fwrite(&psc.domain.ihi[0], sizeof(psc.domain.ihi[0]), 1, file);
  fwrite(&psc.domain.ilo[1], sizeof(psc.domain.ilo[1]), 1, file);
  fwrite(&psc.domain.ihi[1], sizeof(psc.domain.ihi[1]), 1, file);
  fwrite(&psc.domain.ilo[2], sizeof(psc.domain.ilo[2]), 1, file);
  fwrite(&psc.domain.ihi[2], sizeof(psc.domain.ihi[2]), 1, file);
 
  fwrite(psc_output_binary.dowrite_fd, sizeof(bool), 
           NR_EXTRA_FIELDS, file);

  fwrite(datastr, sizeof(char), 4, file);
  
  for (int m = 0; m < NR_EXTRA_FIELDS; m++) {
    if ( psc_output_binary.dowrite_fd[m] ) {
      fwrite(f->all[m], sizeof(float), f->size, file);
    }
  }

  fclose(file);
}


static void
output_binary_field()
{
  if (!psc.pfd.all[0]) {
    init_output_fields(&psc.pfd);
    init_output_fields(&psc.tfd);
  }

  static int pr;
  if (!pr) {
    pr = prof_register("binary_out_field", 1., 0, 0);
  }
  prof_start(pr);

  calculate_pfields(&psc.pfd);

  if (psc_output_binary.dowrite_pfield) {
    if (psc.timestep >= psc_output_binary.pfield_next) {
       psc_output_binary.pfield_next += psc_output_binary.pfield_step;
       binary_field_output_aux(&psc.pfd, "pfd");
    }
  }

  if (psc_output_binary.dowrite_tfield) {
    accumulate_tfields(&psc.pfd, &psc.tfd);
    if (psc.timestep >= psc_output_binary.tfield_next) {
       psc_output_binary.tfield_next += psc_output_binary.tfield_step;
       mean_tfields(&psc.tfd);
       binary_field_output_aux(&psc.tfd, "tfd");
       reset_fields(&psc.tfd);
    }
  }
  
  prof_stop(pr);
}

struct psc_output_ops psc_output_ops_binary = {
  .name           = "binary",
  .create         = output_binary_create,
  .out_field      = output_binary_field,
};

