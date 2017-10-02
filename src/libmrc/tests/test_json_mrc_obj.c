
#include <stdio.h>
#include <assert.h>

#include <mrc_domain.h>
#include <mrc_json.h>

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  
  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);

  mrc_domain_view(domain);
  mrc_json_print(MRC_OBJ_TO_JSON(domain), 0);

  MPI_Finalize();
  
  return 0;
}
