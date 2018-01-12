
#include "mrc_json.h"

#include <json-builder.h> // for conversion to JSON buffer

#include <assert.h>
#include <stdio.h>
#include <string.h>

static struct mrc_json_ops mrc_json_parser_ops;

// ----------------------------------------------------------------------
// helper to get jason_value * out of mrc_json_t

static inline json_value *
get_json_value(mrc_json_t json)
{
  return json.u.parser.value;
}

// ======================================================================
// mrc_json_t based on json-parser library

static int
mrc_json_parser_get_type(mrc_json_t json)
{
  json_value *value = get_json_value(json);
  return value->type;
}

static int
mrc_json_parser_get_integer(mrc_json_t json)
{
  json_value *value = get_json_value(json);
  assert(value->type == json_integer);
  return value->u.integer;
}

static double
mrc_json_parser_get_double(mrc_json_t json)
{
  json_value *value = get_json_value(json);
  assert(value->type == json_double);
  return value->u.dbl;
}

static const char *
mrc_json_parser_get_string(mrc_json_t json)
{
  json_value *value = get_json_value(json);
  assert(value->type == json_string);
  return value->u.string.ptr;
}

static bool
mrc_json_parser_get_boolean(mrc_json_t json)
{
  json_value *value = get_json_value(json);
  assert(value->type == json_boolean);
  return value->u.boolean;
}

static unsigned int
mrc_json_parser_get_object_length(mrc_json_t json)
{
  json_value *value = get_json_value(json);
  assert(value->type == json_object);
  return value->u.object.length;
}

static const char *
mrc_json_parser_get_object_entry_name(mrc_json_t json, unsigned int i)
{
  json_value *value = get_json_value(json);
  assert(value->type == json_object);
  assert(i < value->u.object.length);
  return value->u.object.values[i].name;
}

static mrc_json_t
mrc_json_parser_get_object_entry_value(mrc_json_t json, unsigned int i)
{
  json_value *value = get_json_value(json);
  assert(value->type == json_object);
  assert(i < value->u.object.length);
  return (mrc_json_t) {
    .u.parser.value = value->u.object.values[i].value,
    .ops = &mrc_json_parser_ops,
  };
}

static unsigned int
mrc_json_parser_get_array_length(mrc_json_t json)
{
  json_value *value = get_json_value(json);
  assert(value->type == json_array);
  return value->u.array.length;
}

static mrc_json_t
mrc_json_parser_get_array_entry(mrc_json_t json, unsigned int i)
{
  json_value *value = get_json_value(json);
  assert(value->type == json_array);
  assert(i < value->u.array.length);
  return (mrc_json_t) {
    .u.parser.value = value->u.array.values[i],
    .ops = &mrc_json_parser_ops,
  };
}

char *
mrc_json_to_string(mrc_json_t json)
{
  char *buf = malloc(json_measure(json.u.parser.value));
  json_serialize(buf, json.u.parser.value);
  return buf;
}

static struct mrc_json_ops mrc_json_parser_ops = {
  .get_type               = mrc_json_parser_get_type,
  .get_integer            = mrc_json_parser_get_integer,
  .get_double             = mrc_json_parser_get_double,
  .get_string             = mrc_json_parser_get_string,
  .get_boolean            = mrc_json_parser_get_boolean,
  .get_object_length      = mrc_json_parser_get_object_length,
  .get_object_entry_name  = mrc_json_parser_get_object_entry_name,
  .get_object_entry_value = mrc_json_parser_get_object_entry_value,
  .get_array_length       = mrc_json_parser_get_array_length,
  .get_array_entry        = mrc_json_parser_get_array_entry,
};


mrc_json_t
mrc_json_from_json_parser(void *value)
{
  mrc_json_t json = {
    .u.parser.value = value,
    .ops            = &mrc_json_parser_ops,
  };
  return json;
}

mrc_json_t
mrc_json_parse(const char *buf)
{
  json_settings settings = {};
  settings.value_extra = json_builder_extra;  /* space for json-builder state */

  char error[json_error_max];
  mrc_json_t json = {
    .u.parser.value = json_parse_ex(&settings, buf, strlen(buf), error),
    .ops            = &mrc_json_parser_ops,
  };
  if (!json.u.parser.value) {
    fprintf(stderr, "mrc_json_parse: %s\n", error);
    assert(0);
  }
  return json;
}

mrc_json_t
mrc_json_object_new(unsigned int length)
{
  return mrc_json_from_json_parser(json_object_new(length));
}

mrc_json_t
mrc_json_array_new(unsigned int length)
{
  return mrc_json_from_json_parser(json_array_new(length));
}

mrc_json_t
mrc_json_integer_new(int integer)
{
  return mrc_json_from_json_parser(json_integer_new(integer));
}

mrc_json_t
mrc_json_double_new(double dbl)
{
  return mrc_json_from_json_parser(json_double_new(dbl));
}

mrc_json_t
mrc_json_string_new(const char *str)
{
  return mrc_json_from_json_parser(json_string_new(str));
}

mrc_json_t
mrc_json_boolean_new(bool boolean)
{
  return mrc_json_from_json_parser(json_boolean_new(boolean));
}

mrc_json_t
mrc_json_integer_array_new(unsigned int length, int *arr)
{
  mrc_json_t array = mrc_json_from_json_parser(json_array_new(length));
  for (int i = 0; i < length; i++) {
    mrc_json_array_push(array, mrc_json_integer_new(arr[i]));
  }
  return array;
}

mrc_json_t
mrc_json_double_array_new(unsigned int length, double *arr)
{
  mrc_json_t array = mrc_json_from_json_parser(json_array_new(length));
  for (int i = 0; i < length; i++) {
    mrc_json_array_push_double(array, arr[i]);
  }
  return array;
}

void
mrc_json_object_push(mrc_json_t obj, const char *name, mrc_json_t entry)
{
  assert(obj.ops == &mrc_json_parser_ops);
  assert(entry.ops == &mrc_json_parser_ops);

  json_object_push(get_json_value(obj), name, get_json_value(entry));
}

void
mrc_json_object_push_integer(mrc_json_t obj, const char *name, int integer)
{
  mrc_json_object_push(obj, name, mrc_json_integer_new(integer));
}

void
mrc_json_object_push_double(mrc_json_t obj, const char *name, double dbl)
{
  mrc_json_object_push(obj, name, mrc_json_double_new(dbl));
}

void
mrc_json_object_push_boolean(mrc_json_t obj, const char *name, bool boolean)
{
  mrc_json_object_push(obj, name, mrc_json_boolean_new(boolean));
}

void
mrc_json_object_push_integer_array(mrc_json_t obj, const char *name, unsigned int length, int *arr)
{
  mrc_json_object_push(obj, name, mrc_json_integer_array_new(length, arr));
}

void
mrc_json_object_push_double_array(mrc_json_t obj, const char *name, unsigned int length, double *arr)
{
  mrc_json_object_push(obj, name, mrc_json_double_array_new(length, arr));
}


void
mrc_json_array_push(mrc_json_t arr, mrc_json_t entry)
{
  assert(arr.ops == &mrc_json_parser_ops);
  assert(entry.ops == &mrc_json_parser_ops);

  json_array_push(get_json_value(arr), get_json_value(entry));
}

void
mrc_json_array_push_integer(mrc_json_t arr, int integer)
{
  mrc_json_array_push(arr, mrc_json_integer_new(integer));
}

void
mrc_json_array_push_double(mrc_json_t arr, double dbl)
{
  mrc_json_array_push(arr, mrc_json_double_new(dbl));
}

void
mrc_json_array_push_integer_array(mrc_json_t arr, unsigned int length, int *int_arr)
{
  mrc_json_array_push(arr, mrc_json_integer_array_new(length, int_arr));
}

void
mrc_json_array_push_double_array(mrc_json_t arr, unsigned int length, double *dbl_arr)
{
  mrc_json_array_push(arr, mrc_json_double_array_new(length, dbl_arr));
}
