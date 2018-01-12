
#include "mrc_json.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>

// ======================================================================
// mrc_json dispatch

int
mrc_json_get_type(mrc_json_t json)
{
  assert(json.ops && json.ops->get_type);
  return json.ops->get_type(json);
}

int
mrc_json_get_integer(mrc_json_t json)
{
  assert(json.ops && json.ops->get_integer);
  return json.ops->get_integer(json);
}

double
mrc_json_get_double(mrc_json_t json)
{
  assert(json.ops && json.ops->get_double);
  return json.ops->get_double(json);
}

const char *
mrc_json_get_string(mrc_json_t json)
{
  assert(json.ops && json.ops->get_string);
  return json.ops->get_string(json);
}

bool
mrc_json_get_boolean(mrc_json_t json)
{
  assert(json.ops && json.ops->get_boolean);
  return json.ops->get_boolean(json);
}

unsigned int
mrc_json_get_object_length(mrc_json_t json)
{
  assert(json.ops && json.ops->get_object_length);
  return json.ops->get_object_length(json);
}

const char *
mrc_json_get_object_entry_name(mrc_json_t json, unsigned int i)
{
  assert(json.ops && json.ops->get_object_entry_name);
  return json.ops->get_object_entry_name(json, i);
}

mrc_json_t
mrc_json_get_object_entry_value(mrc_json_t json, unsigned int i)
{
  assert(json.ops && json.ops->get_object_entry_value);
  return json.ops->get_object_entry_value(json, i);
}

unsigned int
mrc_json_get_array_length(mrc_json_t json)
{
  assert(json.ops && json.ops->get_array_length);
  return json.ops->get_array_length(json);
}

mrc_json_t
mrc_json_get_array_entry(mrc_json_t json, unsigned int i)
{
  assert(json.ops && json.ops->get_array_entry);
  return json.ops->get_array_entry(json, i);
}

// ======================================================================
// mrc_json_print

static void
print_indent(int depth)
{
  for (int j = 0; j < depth; j++) {
    printf("  ");
  }
}

static void
mrc_json_print_object(mrc_json_t value, int depth)
{
  print_indent(depth);
  printf("{\n");

  int length = mrc_json_get_object_length(value);
  for (int i = 0; i < length; i++) {
    print_indent(depth+1);
    printf("(name) %s :\n", mrc_json_get_object_entry_name(value, i));
    mrc_json_print(mrc_json_get_object_entry_value(value, i), depth+2);
  }

  print_indent(depth);
  printf("}\n");
}

static void
mrc_json_print_array(mrc_json_t value, int depth)
{
  print_indent(depth);
  printf("[\n");

  int length = mrc_json_get_array_length(value);
  for (int i = 0; i < length; i++) {
    mrc_json_print(mrc_json_get_array_entry(value, i), depth+1);
  }

  print_indent(depth);
  printf("]\n");
}

void
mrc_json_print(mrc_json_t value, unsigned int depth)
{
  int type = mrc_json_get_type(value);

  switch (type) {
  case MRC_JSON_NONE:
    print_indent(depth);
    printf("(none)\n");
    break;
  case MRC_JSON_OBJECT:
    mrc_json_print_object(value, depth + 1);
    break;
  case MRC_JSON_ARRAY:
    mrc_json_print_array(value, depth + 1);
    break;
  case MRC_JSON_INTEGER:
    print_indent(depth);
    printf("(int) %d\n", mrc_json_get_integer(value));
    break;
  case MRC_JSON_DOUBLE:
    print_indent(depth);
    printf("(double) %g\n", mrc_json_get_double(value));
    break;
  case MRC_JSON_STRING:
    print_indent(depth);
    printf("(string) \"%s\"\n", mrc_json_get_string(value));
    break;
  case MRC_JSON_BOOLEAN:
    print_indent(depth);
    printf("(boolean) %s\n", mrc_json_get_boolean(value) ? "true" : "false");
    break;
  default:
    fprintf(stderr, "json_print: unknown type %d\n", type);
    assert(0);
  };
}

// ======================================================================
// helpers

void
mrc_json_get_int3(mrc_json_t json, int arr[3])
{
  for (int d = 0; d < 3; d++) {
    arr[d] = mrc_json_get_array_entry_integer(json, d);
  }
}

void
mrc_json_get_double3(mrc_json_t json, double arr[3])
{
  for (int d = 0; d < 3; d++) {
    arr[d] = mrc_json_get_array_entry_double(json, d);
  }
}

void
mrc_json_get_float3(mrc_json_t json, float arr[3])
{
  for (int d = 0; d < 3; d++) {
    arr[d] = mrc_json_get_array_entry_double(json, d);
  }
}

mrc_json_t
mrc_json_get_object_entry(mrc_json_t json, const char *name)
{
  int len = mrc_json_get_object_length(json);
  for (int i = 0; i < len; i++) {
    const char *entry_name = mrc_json_get_object_entry_name(json, i);
    if (strcmp(name, entry_name) == 0) {
      return mrc_json_get_object_entry_value(json, i);
    }
  }
  return (mrc_json_t) {};
}

int
mrc_json_get_object_entry_integer(mrc_json_t json, const char *name)
{
  mrc_json_t entry = mrc_json_get_object_entry(json, name);
  if (!entry.ops) {
    fprintf(stderr, "%s: entry \"%s\" not found!\n", __func__, name);
    assert(0);
  }
  return mrc_json_get_integer(entry);
}


double
mrc_json_get_object_entry_double(mrc_json_t json, const char *name)
{
  mrc_json_t entry = mrc_json_get_object_entry(json, name);
  if (!entry.ops) {
    fprintf(stderr, "%s: entry \"%s\" not found!\n", __func__, name);
    assert(0);
  }
  return mrc_json_get_double(entry);
}

void
mrc_json_get_object_entry_int3(mrc_json_t json, const char *name, int arr[3])
{
  mrc_json_t entry = mrc_json_get_object_entry(json, name);
  if (!entry.ops) {
    fprintf(stderr, "%s: entry \"%s\" not found!\n", __func__, name);
    assert(0);
  }
  mrc_json_get_int3(entry, arr);
}

void
mrc_json_get_object_entry_double3(mrc_json_t json, const char *name, double arr[3])
{
  mrc_json_t entry = mrc_json_get_object_entry(json, name);
  if (!entry.ops) {
    fprintf(stderr, "%s: entry \"%s\" not found!\n", __func__, name);
    assert(0);
  }
  mrc_json_get_double3(entry, arr);
}

int
mrc_json_get_array_entry_integer(mrc_json_t json, unsigned int i)
{
  return mrc_json_get_integer(mrc_json_get_array_entry(json, i));
}

double
mrc_json_get_array_entry_double(mrc_json_t json, unsigned int i)
{
  return mrc_json_get_double(mrc_json_get_array_entry(json, i));
}


