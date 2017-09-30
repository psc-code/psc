
#include <stdio.h>
#include <assert.h>

#include <mrc_json.h>

int
main(int argc, char **argv)
{
  // integer
  struct mrc_json_value json_int = {
    .type = MRC_JSON_INTEGER,
    .u.integer = 7,
  };
  mrc_json_print(&json_int, 0);

  // double
  struct mrc_json_value json_dbl = {
    .type = MRC_JSON_DOUBLE,
    .u.dbl = 3.33,
  };
  mrc_json_print(&json_dbl, 0);

  // object
  struct mrc_json_object_entry json_obj_entries[] = {
    { .name = "my_integer", .value = &json_int },
    { .name = "my_double", .value = &json_dbl },
  };
  struct mrc_json_value json_obj = {
    .type = MRC_JSON_OBJECT,
    .u.obj.entries = json_obj_entries,
    .u.obj.length  = 2,
  };
  mrc_json_print(&json_obj, 0);

  return 0;
}
