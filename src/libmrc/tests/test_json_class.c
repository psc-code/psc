
#include <stdio.h>
#include <assert.h>

#include <mrc_json.h>

struct test {
  int i;
  int j;
  double d;
  struct mrc_json_value *sub;
  
  struct mrc_json_value json;
};

struct mrc_json_class_entry test_descr[] = {
  { .name = "i"  , .off = offsetof(struct test, i  ), .value.type = MRC_JSON_INTEGER },
  { .name = "j"  , .off = offsetof(struct test, j  ), .value.type = MRC_JSON_INTEGER },
  { .name = "d"  , .off = offsetof(struct test, d  ), .value.type = MRC_JSON_DOUBLE  },
  { .name = "sub", .off = offsetof(struct test, sub), .value.type = MRC_JSON_REF_CLASS },
  {},
};

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

  // class
  struct test test = {
    .i = 99,
    .j = 22,
    .d = 1.1,
    .sub = &json_obj,

    .json.type = MRC_JSON_OBJECT_CLASS,
    .json.u.cls.off = offsetof(struct test, json),
    .json.u.cls.descr = test_descr,
  };

  mrc_json_print(&test.json, 0);
  
  return 0;
}
