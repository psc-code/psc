
#include "mrc_json.h"

#include <mrc_obj.h>

#include <stdio.h>
#include <assert.h>

enum {
  MRC_JSON_INT3    = 100,
  MRC_JSON_DOUBLE3,
};

static struct mrc_json_ops mrc_json_mrc_obj_ops;

static unsigned int
mrc_descr_length(struct param *params)
{
  if (!params) {
    return 0;
  }

  int cnt;
  for (cnt = 0; params[cnt].name; cnt++)
    ;

  return cnt;
}

static mrc_json_t
mrc_descr_to_mrc_json_t(struct param *param, char *p)
{
  char *var = p + (unsigned long) param->var;
  switch (param->type) {
  case PT_INT:
  case MRC_VAR_INT: {
    mrc_json_t json = {
      .u.mrc.v.integer = * (int *) var,
      .u.mrc.type = MRC_JSON_INTEGER,
      .ops = &mrc_json_mrc_obj_ops,
    };
    return json;
  }
  case PT_BOOL: {
    mrc_json_t  json = {
      .u.mrc.v.boolean = * (bool *) var,
      .u.mrc.type = MRC_JSON_BOOLEAN,
      .ops = &mrc_json_mrc_obj_ops,
    };
    return json;
  }
  case PT_DOUBLE:
  case MRC_VAR_DOUBLE: {
    mrc_json_t  json = {
      .u.mrc.v.dbl = * (double *) var,
      .u.mrc.type = MRC_JSON_DOUBLE,
      .ops = &mrc_json_mrc_obj_ops,
    };
    return json;
  }
  case PT_INT3: {
    mrc_json_t  json = {
      .u.mrc.v.int3_ptr = (int *) var,
      .u.mrc.type = MRC_JSON_INT3,
      .ops = &mrc_json_mrc_obj_ops,
    };
    return json;
  }
  case PT_DOUBLE3:
  case MRC_VAR_DOUBLE3: {
    mrc_json_t  json = {
      .u.mrc.v.dbl3_ptr = (double *) var,
      .u.mrc.type = MRC_JSON_DOUBLE3,
      .ops = &mrc_json_mrc_obj_ops,
    };
    return json;
  }
  case MRC_VAR_OBJ: {
    mrc_json_t  json = {
      .u.mrc.v.obj = *(struct mrc_obj **) var,
      .u.mrc.type = MRC_JSON_OBJECT,
      .ops = &mrc_json_mrc_obj_ops,
    };
    return json;
  }
  case PT_OBJ: // FIXME? problems with infinite recursion
  case PT_SELECT: // FIXME
  case PT_INT_ARRAY: { // FIXME
    mrc_json_t  json = {
      .u.mrc.type = MRC_JSON_NONE,
      .ops = &mrc_json_mrc_obj_ops,
    };
    return json;
  }
  default:
    fprintf(stderr, "mrc_descr_to_mrc_json_t: unhandled type %d\n", param->type);
    assert(0);
  }
}

static int
mrc_json_mrc_obj_get_type(mrc_json_t json)
{
  switch (json.u.mrc.type) {
  case MRC_JSON_NONE:
  case MRC_JSON_OBJECT:
  case MRC_JSON_ARRAY:
  case MRC_JSON_INTEGER:
  case MRC_JSON_DOUBLE:
  case MRC_JSON_STRING:
  case MRC_JSON_BOOLEAN:
    return json.u.mrc.type;

  case MRC_JSON_INT3:
  case MRC_JSON_DOUBLE3:
    return MRC_JSON_ARRAY;

  default:
    fprintf(stderr, "mrc_json_mrc_obj_get_type: unhandled type %d\n", json.u.mrc.type);
    assert(0);
  }
}

static int
mrc_json_mrc_obj_get_integer(mrc_json_t json)
{
  assert(json.u.mrc.type == MRC_JSON_INTEGER);
  return json.u.mrc.v.integer;
}

static double
mrc_json_mrc_obj_get_double(mrc_json_t json)
{
  assert(json.u.mrc.type == MRC_JSON_DOUBLE);
  return json.u.mrc.v.dbl;
}

static const char *
mrc_json_mrc_obj_get_string(mrc_json_t json)
{
  assert(json.u.mrc.type == MRC_JSON_STRING);
  return json.u.mrc.v.string;
}

static bool
mrc_json_mrc_obj_get_boolean(mrc_json_t json)
{
  assert(json.u.mrc.type == MRC_JSON_BOOLEAN);
  return json.u.mrc.v.boolean;
}

static unsigned int
mrc_json_mrc_obj_get_object_length(mrc_json_t json)
{
  struct mrc_obj *obj = json.u.mrc.v.obj;

  int cnt = 2; // name and type

  if (obj->ops) { // subtype
    cnt++;
  }

  cnt += mrc_descr_length(obj->cls->param_descr);

  if (obj->ops) {
    cnt += mrc_descr_length(obj->ops->param_descr);
  }
  
  return cnt;
}

static const char *
mrc_json_mrc_obj_get_object_entry_name(mrc_json_t json, unsigned int i)
{
  struct mrc_obj *obj = json.u.mrc.v.obj;

  if (i == 0) { // name
    return "mrc_obj_name";
  }
  i--;

  if (i == 0) { // type
    return "mrc_obj_type";
  }
  i--;

  if (obj->ops) {
    if (i == 0) { // subtype
      return "mrc_obj_subtype";
    }
    i--;
  }

  int len = mrc_descr_length(obj->cls->param_descr);
  if (i < len) {
    return obj->cls->param_descr[i].name;
  }
  i -= len;

  if (obj->ops) {
    int len = mrc_descr_length(obj->ops->param_descr);
    if (i < len) {
      return obj->ops->param_descr[i].name;
    }
    i -= len;
  }

  assert(0);
}

static mrc_json_t
mrc_json_mrc_obj_get_object_entry_value(mrc_json_t json, unsigned int i)
{
  struct mrc_obj *obj = json.u.mrc.v.obj;

  if (i == 0) { // name
    return (mrc_json_t) {
      .u.mrc.v.string = obj->name,
      .u.mrc.type     = MRC_JSON_STRING,
      .ops            = &mrc_json_mrc_obj_ops,
    };
  }
  i--;

  if (i == 0) { // type
    return (mrc_json_t) {
      .u.mrc.v.string = obj->cls->name,
      .u.mrc.type     = MRC_JSON_STRING,
      .ops            = &mrc_json_mrc_obj_ops,
    };
  }
  i--;

  if (obj->ops) {
    if (i == 0) { // subtype
      return (mrc_json_t) {
	.u.mrc.v.string = obj->ops->name,
	.u.mrc.type     = MRC_JSON_STRING,
	.ops            = &mrc_json_mrc_obj_ops,
      };
    }
    i--;
  }

  int len = mrc_descr_length(obj->cls->param_descr);
  if (i < len) {
    char *p = (char *) obj + obj->cls->param_offset;
    return mrc_descr_to_mrc_json_t(&obj->cls->param_descr[i], p);
  }
  i -= len;

  if (obj->ops) {
    int len = mrc_descr_length(obj->ops->param_descr);
    if (i < len) {
      char *p = (char *) obj->subctx + obj->ops->param_offset;
      return mrc_descr_to_mrc_json_t(&obj->ops->param_descr[i], p);
    }
    i -= len;
  }

  assert(0);
}

static unsigned int
mrc_json_mrc_obj_get_array_length(mrc_json_t json)
{
  switch (json.u.mrc.type) {
  case MRC_JSON_INT3:
  case MRC_JSON_DOUBLE3:
    return 3;
  default:
    fprintf(stderr, "mrc_json_mrc_obj_get_array_length: unhandled type %d\n", json.u.mrc.type);
    assert(0);
  }
}

static mrc_json_t
mrc_json_mrc_obj_get_array_entry(mrc_json_t json, unsigned int i)
{
  switch (json.u.mrc.type) {
  case MRC_JSON_INT3:
    return (mrc_json_t) {
      .u.mrc.type      = MRC_JSON_INTEGER,
      .u.mrc.v.integer = json.u.mrc.v.int3_ptr[i],
      .ops             = &mrc_json_mrc_obj_ops,
    };
  case MRC_JSON_DOUBLE3:
    return (mrc_json_t) {
      .u.mrc.type      = MRC_JSON_DOUBLE,
      .u.mrc.v.dbl     = json.u.mrc.v.dbl3_ptr[i],
      .ops             = &mrc_json_mrc_obj_ops,
    };
  default:
    fprintf(stderr, "mrc_json_mrc_obj_get_array_entry: unhandled type %d\n", json.u.mrc.type);
    assert(0);
  }
}

static struct mrc_json_ops mrc_json_mrc_obj_ops = {
  .get_type               = mrc_json_mrc_obj_get_type,
  .get_integer            = mrc_json_mrc_obj_get_integer,
  .get_double             = mrc_json_mrc_obj_get_double,
  .get_string             = mrc_json_mrc_obj_get_string,
  .get_boolean            = mrc_json_mrc_obj_get_boolean,
  .get_object_length      = mrc_json_mrc_obj_get_object_length,
  .get_object_entry_name  = mrc_json_mrc_obj_get_object_entry_name,
  .get_object_entry_value = mrc_json_mrc_obj_get_object_entry_value,
  .get_array_length       = mrc_json_mrc_obj_get_array_length,
  .get_array_entry        = mrc_json_mrc_obj_get_array_entry,
};

mrc_json_t
mrc_obj_to_json(struct mrc_obj *obj)
{
  return (mrc_json_t) {
    .u.mrc.v.obj = obj,
    .u.mrc.type  = MRC_JSON_OBJECT,
    .ops         = &mrc_json_mrc_obj_ops,
  };
}

