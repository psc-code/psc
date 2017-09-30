
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

  // string
  struct mrc_json_value json_str = {
    .type = MRC_JSON_STRING,
    .u.str = "my test string",
  };
  mrc_json_print(&json_str, 0);

  // object
  struct mrc_json_object_entry json_obj_entries[] = {
    { .name = "my_integer", .value = &json_int },
    { .name = "my_double", .value = &json_dbl },
    { .name = "my_string", .value = &json_str },
  };
  struct mrc_json_value json_obj = {
    .type = MRC_JSON_OBJECT,
    .u.obj.entries = json_obj_entries,
    .u.obj.length  = 3,
  };
  mrc_json_print(&json_obj, 0);

  // array
  struct mrc_json_value *json_arr_entries[] = {
    &json_int,
    &json_dbl,
    &json_str,
  };
  struct mrc_json_value json_arr = {
    .type = MRC_JSON_ARRAY,
    .u.arr.entries = json_arr_entries,
    .u.arr.length  = 3,
  };
  mrc_json_print(&json_arr, 0);

  return 0;
}
