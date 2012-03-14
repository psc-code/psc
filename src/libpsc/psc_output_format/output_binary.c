
#include "psc.h"
#include "psc_output_fields_c.h"
#include "psc_output_format_private.h"
#include <mrc_profile.h>
#include <mrc_params.h>

#include <mpi.h>
#include <string.h>

static void
binary_write_field(FILE *file, fields_c_t *fld)
{
  struct psc_patch *patch = &ppsc->patch[0];
  int *ilo = patch->off, ihi[3] = { patch->off[0] + patch->ldims[0],
				    patch->off[1] + patch->ldims[1],
				    patch->off[2] + patch->ldims[2] };

  unsigned int sz = (ihi[0] - ilo[0]) * (ihi[1] - ilo[1]) * (ihi[2] - ilo[2]);

  // convert to float, drop ghost points
  float *data = calloc(sz, sizeof(float));
  int i = 0;
  psc_foreach_patch(ppsc, patch) {
    psc_foreach_3d(ppsc, patch, ix, iy, iz, 0, 0) {
      data[i++] = F3_C(fld,0, ix,iy,iz);
    } foreach_3d_end;
  }
    
  fwrite(data, sizeof(*data), sz, file);

  free(data);
}

static void
psc_output_format_binary_write_fields(struct psc_output_format *format,
				      struct psc_output_fields_c *out,
				      struct psc_fields_list *list,
				      const char *pfx)
{
  struct psc_patch *patch = &ppsc->patch[0];
  int ihi[3] = { patch->off[0] + patch->ldims[0],
		 patch->off[1] + patch->ldims[1],
		 patch->off[2] + patch->ldims[2] };
  const char headstr[] = "PSC ";
  const char datastr[] = "DATA";

  // appears as "?BL?" if NO byte swapping required, ?LB? if required
  unsigned int magic_big_little = 1061962303;    
  unsigned int output_version = 1;
  
  float t_float;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  char filename[strlen(out->data_dir) + 30];
  sprintf(filename, "%s/%s_%06d_%07d.psc", out->data_dir, pfx, rank, ppsc->timestep);

  FILE *file = fopen(filename, "wb");

  // Header  
  fwrite(headstr, sizeof(char), 4, file);
  fwrite(&magic_big_little, sizeof(unsigned int), 1, file);
  fwrite(&output_version, sizeof(unsigned int), 1, file);

  t_float = (float) ppsc->dx[0];  fwrite(&t_float, sizeof(float), 1, file);
  t_float = (float) ppsc->dx[1];  fwrite(&t_float, sizeof(float), 1, file);
  t_float = (float) ppsc->dx[2];  fwrite(&t_float, sizeof(float), 1, file);
  t_float = (float) ppsc->dt;     fwrite(&t_float, sizeof(float), 1, file);

  // Indices on local proc
  fwrite(&patch->off[0], sizeof(patch->off[0]), 1, file);
  fwrite(&ihi[0], sizeof(ihi[0]), 1, file);
  fwrite(&patch->off[1], sizeof(patch->off[1]), 1, file);
  fwrite(&ihi[1], sizeof(ihi[1]), 1, file);
  fwrite(&patch->off[2], sizeof(patch->off[2]), 1, file);
  fwrite(&ihi[2], sizeof(ihi[2]), 1, file);

  // Globally saved indices (everything for now...)
  int glo[3] = {};
  fwrite(&glo[0], sizeof(glo[0]), 1, file);
  fwrite(&ppsc->domain.gdims[0], sizeof(ppsc->domain.gdims[0]), 1, file);
  fwrite(&glo[1], sizeof(glo[1]), 1, file);
  fwrite(&ppsc->domain.gdims[1], sizeof(ppsc->domain.gdims[1]), 1, file);
  fwrite(&glo[2], sizeof(glo[2]), 1, file);
  fwrite(&ppsc->domain.gdims[2], sizeof(ppsc->domain.gdims[2]), 1, file);

  fwrite(&list->nr_flds, sizeof(list->nr_flds), 1, file);
  for (int i = 0; i < list->nr_flds; i++) {
    char fldname[8] = {};
    snprintf(fldname, 8, "%s", psc_mfields_comp_name(list->flds[i], 0));
    fwrite(fldname, 8, 1, file); 
  }

  fwrite(datastr, sizeof(char), 4, file);
  
  for (int m = 0; m < list->nr_flds; m++) {
    binary_write_field(file, psc_mfields_get_patch_c(list->flds[m], 0));
  }

  fclose(file);
}

// ======================================================================
// psc_output_format: subclass "binary"

struct psc_output_format_ops psc_output_format_binary_ops = {
  .name                  = "binary",
  .write_fields          = psc_output_format_binary_write_fields,
};



