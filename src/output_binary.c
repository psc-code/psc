
#include "psc.h"

#include "util/profile.h"
#include "util/params.h"

#include <mpi.h>
#include <string.h>

struct psc_output_binary {
  char *data_dir;
  int field_next_out;
  int field_step;
  int dowrite_ni, dowrite_ne, dowrite_nn;
  int dowrite_ex, dowrite_ey, dowrite_ez;
  int dowrite_bx, dowrite_by, dowrite_bz;
  int dowrite_jx, dowrite_jy, dowrite_jz;
};

#define VAR(x) (void *)offsetof(struct psc_output_binary, x)

static struct param psc_binary_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(NULL)   },
  { "field_first_out"    , VAR(field_next_out)       , PARAM_INT(0)        },
  { "field_step_out"     , VAR(field_step)           , PARAM_INT(10)        },
  { "output_write_ne"    , VAR(dowrite_ni)           , PARAM_BOOL(1)        },
  { "output_write_ni"    , VAR(dowrite_ne)           , PARAM_BOOL(1)        },
  { "output_write_nn"    , VAR(dowrite_nn)           , PARAM_BOOL(1)        },
  { "output_write_ex"    , VAR(dowrite_ex)           , PARAM_BOOL(1)        },
  { "output_write_ey"    , VAR(dowrite_ey)           , PARAM_BOOL(1)        },
  { "output_write_ez"    , VAR(dowrite_ez)           , PARAM_BOOL(1)        },
  { "output_write_bx"    , VAR(dowrite_bx)           , PARAM_BOOL(1)        },
  { "output_write_by"    , VAR(dowrite_by)           , PARAM_BOOL(1)        },
  { "output_write_bz"    , VAR(dowrite_bz)           , PARAM_BOOL(1)        },
  { "output_write_jx"    , VAR(dowrite_jx)           , PARAM_BOOL(1)        },
  { "output_write_jy"    , VAR(dowrite_jy)           , PARAM_BOOL(1)        },
  { "output_write_jz"    , VAR(dowrite_jz)           , PARAM_BOOL(1)        },
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
output_binary_field()
{
  char *head = "PSC ";

  // appears as "?BL?" if NO byte swapping required, ?LB? if required
  unsigned int magic_big_little = 1061962303;    
  unsigned int output_version = 1;

  int dowrite_fd[12] = { psc_output_binary.dowrite_ni, 
                         psc_output_binary.dowrite_ne,
                         psc_output_binary.dowrite_nn,
                         psc_output_binary.dowrite_ex,
                         psc_output_binary.dowrite_ey,
                         psc_output_binary.dowrite_ez,
                         psc_output_binary.dowrite_bx,
                         psc_output_binary.dowrite_by,
                         psc_output_binary.dowrite_bz,
                         psc_output_binary.dowrite_jx,
                         psc_output_binary.dowrite_jy,
                         psc_output_binary.dowrite_jz,};

  if (psc.timestep < psc_output_binary.field_next_out) {
    return;
  }
  psc_output_binary.field_next_out += psc_output_binary.field_step;

  static int pr;
  if (!pr) {
    pr = prof_register("binary_out_field", 1., 0, 0);
  }
  prof_start(pr);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char filename[200];
  sprintf(filename, "data/pfd_%03d_%07d.psc", rank, psc.timestep);

  FILE *file = fopen(filename, "wb");

  // Header  
  fwrite(head, sizeof(char), 4, file);
  fwrite(&magic_big_little, sizeof(unsigned int), 1, file);
  fwrite(&output_version, sizeof(unsigned int), 1, file);

  fwrite(&psc.dx[0], sizeof(psc.dx[0]), 1, file);
  fwrite(&psc.dx[1], sizeof(psc.dx[1]), 1, file);
  fwrite(&psc.dx[2], sizeof(psc.dx[2]), 1, file);
  fwrite(&psc.dt, sizeof(psc.dt), 1, file);

  fwrite(&psc.domain.ilo[0], sizeof(psc.domain.ilo[0]), 1, file);
  fwrite(&psc.domain.ihi[0], sizeof(psc.domain.ihi[0]), 1, file);
  fwrite(&psc.domain.ilo[1], sizeof(psc.domain.ilo[1]), 1, file);
  fwrite(&psc.domain.ihi[1], sizeof(psc.domain.ihi[1]), 1, file);
  fwrite(&psc.domain.ilo[2], sizeof(psc.domain.ilo[2]), 1, file);
  fwrite(&psc.domain.ihi[2], sizeof(psc.domain.ihi[2]), 1, file);

  fwrite(dowrite_fd, sizeof(int), 12, file);
  
  // pack buffer for output
  size_t fsize = (psc.domain.ihi[0] - psc.domain.ilo[0])
                 * (psc.domain.ihi[1] - psc.domain.ilo[1])
                 * (psc.domain.ihi[2] - psc.domain.ilo[2]);
  float fd_buf[fsize];

  for (int m = 0; m < NR_FIELDS; m++) {

    if ( dowrite_fd[m] ) {

      int j = 0;
      for (int iz = psc.ilo[2]; iz < psc.ihi[2]; iz++)
        for (int iy = psc.ilo[1]; iy < psc.ihi[1]; iy++)
           for (int ix = psc.ilo[0]; ix < psc.ihi[0]; ix++)
	     fd_buf[j++] = (float) FF3(m, ix,iy,iz);

      fwrite(fd_buf, sizeof(float), fsize, file);
    }
  }

  fclose(file);

  prof_stop(pr);
}

struct psc_output_ops psc_output_ops_binary = {
  .name           = "binary",
  .create         = output_binary_create,
  .out_field      = output_binary_field,
};

