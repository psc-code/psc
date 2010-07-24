#include "psc.h"
#include "output_fields.h"
#include "util/params.h"
#include "util/profile.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

struct psc_extra_fields {
  unsigned int size;
  unsigned int naccum;
  float *all[NR_EXTRA_FIELDS];
};

static const char *x_fldname[NR_EXTRA_FIELDS] = {
  [X_NE]   = "ne"  , [X_NI]   = "ni"  , [X_NN]   = "nn",
  [X_JXI]  = "jx"  , [X_JYI]  = "jy"  , [X_JZI]  = "jz",
  [X_EX]   = "ex"  , [X_EY]   = "ey"  , [X_EZ]   = "ez",
  [X_BX]   = "bx"  , [X_BY]   = "by"  , [X_BZ]   = "bz",
  [X_JXEX] = "jxex", [X_JYEY] = "jyey", [X_JZEZ] = "jzez",
  [X_POYX] = "poyx", [X_POYY] = "poyy", [X_POYZ] = "poyz",
  [X_E2X ] = "e2x" , [X_E2Y]  = "e2y" , [X_E2Z]  = "e2z",
  [X_B2X ] = "b2x" , [X_B2Y]  = "b2y" , [X_B2Z]  = "b2z",
};

static void
reset_fields(struct psc_extra_fields *f)
{
  f->naccum = 0;

  for (int m = 0; m < NR_EXTRA_FIELDS; m++) {
    for (int j = 0; j < f->size; j++)  {
      f->all[m][j] = 0.0;
    }
  }
}

static void
init_output_fields(struct psc_extra_fields *f)
{
  f->size = ((psc.ihi[0]-psc.ilo[0]) *
	     (psc.ihi[1]-psc.ilo[1]) *
	     (psc.ihi[2]-psc.ilo[2]));

  for (int m = 0; m < NR_EXTRA_FIELDS; m++) {
    f->all[m] = malloc(f->size * sizeof(float));
  }

  reset_fields(f);
}

#if 0 // unused
static void
free_output_fields(struct psc_extra_fields *f)
{
  for (int m = 0; m < NR_EXTRA_FIELDS ; m++)
    free(f->all[m]);
}
#endif

static void
calculate_pfields(struct psc_extra_fields *p)
{
   int j = 0;

   for (int iz = psc.ilo[2]; iz < psc.ihi[2]; iz++) {
      for (int iy = psc.ilo[1]; iy < psc.ihi[1]; iy++) {
        for (int ix = psc.ilo[0]; ix < psc.ihi[0]; ix++) {
	  p->all[X_NE][j] = FF3(NE,ix,iy,iz);
	  p->all[X_NI][j] = FF3(NI,ix,iy,iz);
	  p->all[X_NN][j] = FF3(NN,ix,iy,iz);

           p->all[X_EX][j] = .5f * ( FF3(EX,ix,iy,iz)
                                    +FF3(EX,ix-1,iy,iz));
           p->all[X_EY][j] = .5f * ( FF3(EX,ix,iy,iz)
                                    +FF3(EY,ix,iy-1,iz));
           p->all[X_EZ][j] = .5f * ( FF3(EZ,ix,iy,iz)
                                    +FF3(EZ,ix,iy,iz-1));

           p->all[X_BX][j] =  .25f * ( FF3(BX,ix,iy,iz)
                                      +FF3(BX,ix,iy-1,iz)
                                      +FF3(BX,ix,iy,iz-1) 
  			              +FF3(BX,ix,iy-1,iz-1));
           p->all[X_BY][j] =  .25f * ( FF3(BY,ix,iy,iz)
                                      +FF3(BY,ix-1,iy,iz)
                                      +FF3(BY,ix,iy,iz-1) 
				      +FF3(BY,ix-1,iy,iz-1));
           p->all[X_BZ][j] =  .25f * ( FF3(BZ,ix,iy,iz)
                                      +FF3(BZ,ix-1,iy,iz)
                                      +FF3(BZ,ix,iy-1,iz) 
                                      +FF3(BZ,ix-1,iy-1,iz));

           p->all[X_JXI][j] = .5f * ( FF3(JXI,ix,iy,iz) + FF3(JXI,ix-1,iy,iz) );
           p->all[X_JYI][j] = .5f * ( FF3(JYI,ix,iy,iz) + FF3(JYI,ix,iy-1,iz) );
           p->all[X_JZI][j] = .5f * ( FF3(JZI,ix,iy,iz) + FF3(JZI,ix,iy,iz-1) );

           p->all[X_JXEX][j] = p->all[X_JXI][j] * p->all[X_EX][j];
           p->all[X_JYEY][j] = p->all[X_JYI][j] * p->all[X_EY][j];
           p->all[X_JZEZ][j] = p->all[X_JZI][j] * p->all[X_EZ][j];

           p->all[X_POYX][j] = p->all[X_EY][j] * p->all[X_BZ][j] - p->all[X_EZ][j] * p->all[X_BY][j];
           p->all[X_POYY][j] = p->all[X_EZ][j] * p->all[X_BX][j] - p->all[X_EX][j] * p->all[X_BZ][j];
           p->all[X_POYZ][j] = p->all[X_EX][j] * p->all[X_BY][j] - p->all[X_EY][j] * p->all[X_BX][j];

           p->all[X_E2X][j] = p->all[X_EX][j]*p->all[X_EX][j];
           p->all[X_E2Y][j] = p->all[X_EY][j]*p->all[X_EY][j];
           p->all[X_E2Z][j] = p->all[X_EZ][j]*p->all[X_EZ][j];

           p->all[X_B2X][j] = p->all[X_BX][j]*p->all[X_BX][j];
           p->all[X_B2Y][j] = p->all[X_BY][j]*p->all[X_BY][j];
           p->all[X_B2Z][j] = p->all[X_BZ][j]*p->all[X_BZ][j];

           j++;
         }
      }
    }
}

static void
accumulate_tfields(struct psc_extra_fields *p, struct psc_extra_fields *t)
{
  t->naccum = t->naccum + 1;

  assert(p->size == t->size);
  for (int m = 0; m < NR_EXTRA_FIELDS; m++) {
    for (int j = 0; j < p->size; j++)  {
      t->all[m][j] += p->all[m][j];
    }
  }
}


// convert accumulated values to correct temporal mean
// (divide by naccum)
static void
mean_tfields(struct psc_extra_fields *f)
{
  assert(f->naccum > 0);
  for (int m = 0; m < NR_EXTRA_FIELDS; m++) {
    for (int j = 0; j < f->size; j++) {
      f->all[m][j] /= f->naccum;
    }
  }
}

static struct psc_output_format_ops *psc_output_format_ops_list[] = {
  &psc_output_format_ops_binary,
#ifdef HAVE_LIBHDF5
  &psc_output_format_ops_hdf5,
#endif
  NULL,
};

static struct psc_output_format_ops *
find_output_format_ops(const char *ops_name)
{
  for (int i = 0; psc_output_format_ops_list[i]; i++) {
    if (strcasecmp(psc_output_format_ops_list[i]->name, ops_name) == 0)
      return psc_output_format_ops_list[i];
  }
  fprintf(stderr, "ERROR: psc_output_format_ops '%s' not available.\n", ops_name);
  abort();
}

struct psc_output_c {
  char *data_dir;
  char *output_format;
  bool output_1proc;
  bool dowrite_pfield, dowrite_tfield;
  int pfield_next, tfield_next;
  int pfield_step, tfield_step;
  bool dowrite_fd[NR_EXTRA_FIELDS];

  // storage for output
  struct psc_extra_fields pfd, tfd;

  struct psc_output_format_ops *format_ops;
};

#define VAR(x) (void *)offsetof(struct psc_output_c, x)

static struct param psc_output_c_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(NULL)   },
  { "output_format"      , VAR(output_format)        , PARAM_STRING("binary") },
  { "output_1proc"       , VAR(output_1proc)         , PARAM_BOOL(0)        },
  { "write_pfield"       , VAR(dowrite_pfield)       , PARAM_BOOL(1)        },
  { "pfield_first"       , VAR(pfield_next)          , PARAM_INT(0)         },
  { "pfield_step"        , VAR(pfield_step)          , PARAM_INT(10)        },
  { "write_tfield"       , VAR(dowrite_tfield)       , PARAM_BOOL(1)        },
  { "tfield_first"       , VAR(tfield_next)          , PARAM_INT(0)         },
  { "tfield_step"        , VAR(tfield_step)          , PARAM_INT(10)        },
  { "output_write_ne"    , VAR(dowrite_fd[X_NE])     , PARAM_BOOL(1)        },
  { "output_write_ni"    , VAR(dowrite_fd[X_NI])     , PARAM_BOOL(1)        },
  { "output_write_nn"    , VAR(dowrite_fd[X_NN])     , PARAM_BOOL(1)        },
  { "output_write_jx"    , VAR(dowrite_fd[X_JXI])    , PARAM_BOOL(1)        },
  { "output_write_jy"    , VAR(dowrite_fd[X_JYI])    , PARAM_BOOL(1)        },
  { "output_write_jz"    , VAR(dowrite_fd[X_JZI])    , PARAM_BOOL(1)        },
  { "output_write_ex"    , VAR(dowrite_fd[X_EX])     , PARAM_BOOL(1)        },
  { "output_write_ey"    , VAR(dowrite_fd[X_EY])     , PARAM_BOOL(1)        },
  { "output_write_ez"    , VAR(dowrite_fd[X_EZ])     , PARAM_BOOL(1)        },
  { "output_write_bx"    , VAR(dowrite_fd[X_BX])     , PARAM_BOOL(1)        },
  { "output_write_by"    , VAR(dowrite_fd[X_BY])     , PARAM_BOOL(1)        },
  { "output_write_bz"    , VAR(dowrite_fd[X_BZ])     , PARAM_BOOL(1)        },
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

static struct psc_output_c psc_output_c;

// ----------------------------------------------------------------------
// output_c_create

static void output_c_create(void)
{ 
  struct psc_output_c *out = &psc_output_c;
  params_parse_cmdline(out, psc_output_c_descr, "PSC output C", MPI_COMM_WORLD);
  params_print(out, psc_output_c_descr, "PSC output C", MPI_COMM_WORLD);

  out->format_ops = find_output_format_ops(out->output_format);
  if (out->format_ops->create) {
    out->format_ops->create();
  }
};

// ----------------------------------------------------------------------
// make_fields_list

static void
make_fields_list(struct psc_fields_list *list, struct psc_extra_fields *f,
		 bool *dowrite_fd)
{
  list->nr_flds = 0;
  for (int m = 0; m < NR_EXTRA_FIELDS; m++) {
    if (!dowrite_fd[m])
      continue;

    struct psc_field *fld = &list->flds[list->nr_flds++];
    fld->data = f->all[m];
    fld->size = f->size;
    fld->name = x_fldname[m];
    for (int d = 0; d < 3; d++) {
      fld->ilo[d] = psc.ilo[d];
      fld->ihi[d] = psc.ihi[d];
    }
  }
  list->dowrite_fd = dowrite_fd;
}

// ----------------------------------------------------------------------
// copy_to_global helper

static void
copy_to_global(float *fld, float *buf, int *ilo, int *ihi, int *ilg, int *img)
{
  int my = psc.ghi[1] - psc.glo[1];
  int mx = psc.ghi[0] - psc.glo[0];

  for (int iz = ilo[2]; iz < ihi[2]; iz++) {
    for (int iy = ilo[1]; iy < ihi[1]; iy++) {
      for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	fld[((iz - psc.glo[2]) * my + iy - psc.glo[1]) * mx + ix - psc.glo[0]] =
	  buf[((iz - ilg[2]) * img[1] + iy - ilg[1]) * img[0] + ix - ilg[0]];
      }
    }
  }
}

// ----------------------------------------------------------------------
// write_fields_1proc

void
write_fields_1proc(struct psc_output_format_ops *format_ops,
		   struct psc_fields_list *list, const char *prefix)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  void *ctx;
  if (rank == 0) {
    char datadir[] = "data";
    char filename[strlen(datadir) + 30];
    sprintf(filename, "%s/%s_%07d%s", datadir, prefix, psc.timestep, format_ops->ext);
    printf("[%d] write_fields_1proc: %s\n", rank, filename);
    format_ops->open(list, filename, &ctx);
  }

  /* printf("glo %d %d %d ghi %d %d %d\n", glo[0], glo[1], glo[2], */
  /* 	     ghi[0], ghi[1], ghi[2]); */

  for (int m = 0; m < list->nr_flds; m++) {
    int s_ilo[3], s_ihi[3], s_ilg[3], s_img[3];
    float *s_data = list->flds[m].data;

    for (int d = 0; d < 3; d++) {
      s_ilo[d] = list->flds[m].ilo[d];
      s_ihi[d] = list->flds[m].ihi[d];
      s_ilg[d] = s_ilo[d];
      s_img[d] = s_ihi[d] - s_ilo[d];
    }
    
    if (rank != 0) {
      MPI_Send(s_ilo, 3, MPI_INT, 0, 100, MPI_COMM_WORLD);
      MPI_Send(s_ihi, 3, MPI_INT, 0, 101, MPI_COMM_WORLD);
      MPI_Send(s_ilg, 3, MPI_INT, 0, 102, MPI_COMM_WORLD);
      MPI_Send(s_img, 3, MPI_INT, 0, 103, MPI_COMM_WORLD);
      MPI_Send(s_data, list->flds[m].size, MPI_FLOAT, 0, 104, MPI_COMM_WORLD);
    } else { // rank == 0
      struct psc_field fld;
      fld.name = list->flds[m].name;
      fld.size = 1;
      for (int d = 0; d < 3; d++) {
	fld.ilo[d] = psc.glo[d];
	fld.ihi[d] = psc.ghi[d];
	fld.size *= (psc.ghi[d] - psc.glo[d]);
      }
      fld.data = calloc(fld.size, sizeof(*fld.data));

      for (int n = 0; n < size; n++) {
	int ilo[3], ihi[3], ilg[3], img[3];
	float *buf;
	
	if (n == 0) {
	  for (int d = 0; d < 3; d++) {
	    ilo[d] = s_ilo[d];
	    ihi[d] = s_ihi[d];
	    ilg[d] = s_ilg[d];
	    img[d] = s_img[d];
	  }
	  buf = s_data;
	} else {
	  MPI_Recv(ilo, 3, MPI_INT, n, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(ihi, 3, MPI_INT, n, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(ilg, 3, MPI_INT, n, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(img, 3, MPI_INT, n, 103, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  int ntot = img[0] * img[1] * img[2];
	  buf = calloc(ntot, sizeof(*buf));
	  MPI_Recv(buf, ntot, MPI_FLOAT, n, 104, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	/* printf("[%d] ilo %d %d %d ihi %d %d %d\n", rank, ilo[0], ilo[1], ilo[2], */
	/*        ihi[0], ihi[1], ihi[2]); */
	copy_to_global(fld.data, buf, ilo, ihi, ilg, img);
	if (n != 0) {
	  free(buf);
	}
      }
      format_ops->write_field(ctx, &fld);
      free(fld.data);
    }
  }

  if (rank == 0) {
    format_ops->close(ctx);
  }
}

// ----------------------------------------------------------------------
// write_fields

static void
write_fields(struct psc_output_c *out, struct psc_fields_list *list,
	     const char *prefix)
{
  if (out->output_1proc) {
    return write_fields_1proc(out->format_ops, list, prefix);
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char datadir[] = "data";
  char filename[strlen(datadir) + 30];
  sprintf(filename, "data/%s_%06d_%07d%s", prefix, rank, psc.timestep,
	  out->format_ops->ext);

  void *ctx;
  out->format_ops->open(list, filename, &ctx);

  for (int m = 0; m < list->nr_flds; m++) {
    out->format_ops->write_field(ctx, &list->flds[m]);
  }
  
  out->format_ops->close(ctx);
}

// ----------------------------------------------------------------------
// output_c_field

static void
output_c_field()
{
  struct psc_output_c *out = &psc_output_c;

  if (!out->pfd.all[0]) {
    init_output_fields(&out->pfd);
    init_output_fields(&out->tfd);
  }

  static int pr;
  if (!pr) {
    pr = prof_register("output_c_field", 1., 0, 0);
  }
  prof_start(pr);

  calculate_pfields(&out->pfd);

  if (out->dowrite_pfield) {
    if (psc.timestep >= out->pfield_next) {
       out->pfield_next += out->pfield_step;
       struct psc_fields_list flds_list;
       make_fields_list(&flds_list, &out->pfd, out->dowrite_fd);
       write_fields(out, &flds_list, "pfd");
    }
  }

  if (out->dowrite_tfield) {
    accumulate_tfields(&out->pfd, &out->tfd);
    if (psc.timestep >= out->tfield_next) {
      out->tfield_next += out->tfield_step;
      mean_tfields(&out->tfd);
      struct psc_fields_list flds_list;
      make_fields_list(&flds_list, &out->tfd, out->dowrite_fd);
      write_fields(out, &flds_list, "tfd");
      reset_fields(&out->tfd);
    }
  }
  
  prof_stop(pr);
}

// ======================================================================
// psc_output_ops_c

struct psc_output_ops psc_output_ops_c = {
  .name           = "c",
  .create         = output_c_create,
  .out_field      = output_c_field,
};

