
#include "mrc_json.h"

#include <json-builder.h> // for conversion to JSON buffer

#include <assert.h>

static inline json_value *
get_jason_value(mrc_json_t json)
{
  return json.value;
}

// ======================================================================
// mrc_json_t based on json-parser library

int
mrc_json_get_type(mrc_json_t json)
{
  json_value *value = get_jason_value(json);
  return value->type;
}

int
mrc_json_get_integer(mrc_json_t json)
{
  json_value *value = get_jason_value(json);
  assert(value->type == json_integer);
  return value->u.integer;
}

int
mrc_json_get_double(mrc_json_t json)
{
  json_value *value = get_jason_value(json);
  assert(value->type == json_double);
  return value->u.dbl;
}

const char *
mrc_json_get_string(mrc_json_t json)
{
  json_value *value = get_jason_value(json);
  assert(value->type == json_string);
  return value->u.string.ptr;
}

unsigned int
mrc_json_get_object_length(mrc_json_t json)
{
  json_value *value = get_jason_value(json);
  assert(value->type == json_object);
  return value->u.object.length;
}

const char *
mrc_json_get_object_entry_name(mrc_json_t json, unsigned int i)
{
  json_value *value = get_jason_value(json);
  assert(value->type == json_object);
  assert(i < value->u.object.length);
  return value->u.object.values[i].name;
}

mrc_json_t
mrc_json_get_object_entry_value(mrc_json_t json, unsigned int i)
{
  json_value *value = get_jason_value(json);
  assert(value->type == json_object);
  assert(i < value->u.object.length);
  return (mrc_json_t) { .value = value->u.object.values[i].value };
}

unsigned int
mrc_json_get_array_length(mrc_json_t json)
{
  json_value *value = get_jason_value(json);
  assert(value->type == json_array);
  return value->u.array.length;
}

mrc_json_t
mrc_json_get_array_entry(mrc_json_t json, unsigned int i)
{
  json_value *value = get_jason_value(json);
  assert(value->type == json_array);
  assert(i < value->u.array.length);
  return (mrc_json_t) { .value = value->u.array.values[i] };
}

char *
mrc_json_to_string(mrc_json_t json)
{
  char *buf = malloc(json_measure(json.value));
  json_serialize(buf, json.value);
  return buf;
}



