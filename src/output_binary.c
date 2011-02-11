
#include "psc.h"
#include "output_fields.h"
#include "util/profile.h"
#include <mrc_params.h>

#include <mpi.h>
#include <string.h>

static void
binary_write_field(FILE *file, fields_base_t *fld)
{
  int *ilo = psc.ilo, *ihi = psc.ihi;

  unsigned int sz = (ihi[0] - ilo[0]) * (ihi[1] - ilo[1]) * (ihi[2] - ilo[2]);

  // convert to float, drop ghost points
  float *data = calloc(sz, sizeof(float));
  int i = 0;
  for (int iz = ilo[2]; iz < ihi[2]; iz++) {
    for (int iy = ilo[1]; iy < ihi[1]; iy++) {
      for (int ix = ilo[0]; ix < ihi[0]; ix++) {
	data[i++] = XF3_BASE(fld,0, ix,iy,iz);
      }
    }
  }
    
  fwrite(data, sizeof(*data), sz, file);

  free(data);
}

static void
binary_write_fields(struct psc_output_c *out, struct psc_fields_list *list,
		    const char *pfx)
{
  const char headstr[] = "PSC ";
  const char datastr[] = "DATA";

  // appears as "?BL?" if NO byte swapping required, ?LB? if required
  unsigned int magic_big_little = 1061962303;    
  unsigned int output_version = 1;
  
  float t_float;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  char filename[strlen(out->data_dir) + 30];
  sprintf(filename, "%s/%s_%06d_%07d.psc", out->data_dir, pfx, rank, psc.timestep);

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

  fwrite(&list->nr_flds, sizeof(list->nr_flds), 1, file);
  for (int i = 0; i < list->nr_flds; i++) {
    char fldname[8] = {};
    snprintf(fldname, 8, "%s", list->flds[i].name[0]);
    fwrite(fldname, 8, 1, file); 
  }

  fwrite(datastr, sizeof(char), 4, file);
  
  for (int m = 0; m < list->nr_flds; m++) {
    binary_write_field(file, &list->flds[m]);
  }

  fclose(file);
}

// ======================================================================
// psc_output_format_ops_binary

struct psc_output_format_ops psc_output_format_ops_binary = {
  .name         = "binary",
  .write_fields = binary_write_fields,
};


