
#include <stdio.h>
#include <assert.h>

#include <mrc_json.h>

#include <mrc_domain.h>

struct mrc_json_value *
mrc_obj_to_json(struct mrc_obj *obj)
{
  return &obj->json;
}

#define MRC_OBJ_TO_JSON(obj) mrc_obj_to_json((struct mrc_obj *) obj)

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  
  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);

  mrc_json_print(MRC_OBJ_TO_JSON(domain), 0);

  MPI_Finalize();
  
  return 0;
}
