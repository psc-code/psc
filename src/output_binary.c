
#include "psc.h"
#include "output_fields.h"
#include "util/profile.h"
#include "util/params.h"

#include <mpi.h>
#include <string.h>

struct binary_ctx {
  FILE *file;
};

static void
binary_open(struct psc_fields_list *list, const char *filename, void **pctx)
{
  const char headstr[] = "PSC ";
  const char datastr[] = "DATA";

  // appears as "?BL?" if NO byte swapping required, ?LB? if required
  unsigned int magic_big_little = 1061962303;    
  unsigned int output_version = 1;
  
  float t_float;

  struct binary_ctx *binary = malloc(sizeof(*binary));
  
  binary->file = fopen(filename, "wb");
  FILE *file = binary->file;

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
  
  assert(list->dowrite_fd);
  fwrite(list->dowrite_fd, sizeof(bool), NR_EXTRA_FIELDS, file);

  fwrite(datastr, sizeof(char), 4, file);
  
  *pctx = binary;
}

static void
binary_close(void *ctx)
{
  struct binary_ctx *binary = ctx;

  fclose(binary->file);
}

static void
binary_write_field(void *ctx, struct psc_field *fld)
{
  struct binary_ctx *binary = ctx;

  fwrite(fld->data, sizeof(float), fld->size, binary->file);
}

// ======================================================================
// psc_output_format_ops_binary

struct psc_output_format_ops psc_output_format_ops_binary = {
  .name         = "binary",
  .ext          = ".psc",
  .open         = binary_open,
  .close        = binary_close,
  .write_field  = binary_write_field,
};


