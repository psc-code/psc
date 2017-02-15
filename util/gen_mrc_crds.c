
#include <stdio.h>
#include <stdbool.h>

#include "mrc_params.h"
#include "mrc_domain.h"
#include "mrc_crds.h"
#include "mrc_fld.h"

#define MAX_STRLEN 91

int
main(int argc, char **argv)
{
  const char *run = "test";
  char grid2_fname[MAX_STRLEN];
  char hgrid2_fname[MAX_STRLEN];
  struct mrc_domain *mrc_domain;
  struct mrc_crds *mrc_crds;
  struct mrc_ndarray *nd;
  FILE *fout1;
  FILE *fout2;

  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  mrc_params_get_option_string("run", &run);
  snprintf(grid2_fname, MAX_STRLEN - 1, "%s.grid2", run);
  snprintf(hgrid2_fname, MAX_STRLEN - 1, "%s.hgrid2", run);

  mrc_domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_crds = mrc_domain_get_crds(mrc_domain);
  mrc_domain_set_from_options(mrc_domain);
  mrc_domain_setup(mrc_domain);

  // mrc_domain_view(mrc_domain);

  ///////////////////////////////////////////
  // write out a 'high precision' grid2 file
  fout1 = fopen(grid2_fname, "w");
  fout2 = fopen(hgrid2_fname, "w");
  for (int d=0; d < 3; d++) {
    nd = mrc_crds->global_crd[d];
    
    fprintf(fout1, "%d\n", nd->dims.vals[0]);
    fprintf(fout2, "%d\n", nd->dims.vals[0]);

    for (int i = nd->nd_acc.beg[0]; i < nd->nd_acc.end[0]; i++) {
      fprintf(fout1, "%16.10g\n", MRC_D2(nd, i, 0));
      fprintf(fout2, "%16.10g\n", MRC_D2(nd, i, 1));
    }
  }
  fclose(fout1);
  fclose(fout2);

  mrc_domain_destroy(mrc_domain);
  MPI_Finalize();
  return 0;
}
