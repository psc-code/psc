
#include "mrc_json.h"

#include <assert.h>
#include <stdio.h>

int
mrc_json_get_integer(struct mrc_json_value *value)
{
  if (value->type == MRC_JSON_INTEGER) {
    return value->u.integer;
  } else {
    assert(0);
  }
}

double
mrc_json_get_double(struct mrc_json_value *value)
{
  if (value->type == MRC_JSON_DOUBLE) {
    return value->u.dbl;
  } else {
    assert(0);
  }
}

int
mrc_json_type(struct mrc_json_value *value)
{
  return value->type & MRC_JSON_TYPE_MASK;
}

unsigned int
mrc_json_object_length(struct mrc_json_value *value)
{
  if (value->type == MRC_JSON_OBJECT) {
    return value->u.obj.length;
  } else if (value->type == MRC_JSON_OBJECT_CLASS) {
    int i;
    for (i = 0; value->u.cls.descr[i].name; i++)
      ;
    return i;
  } else {
    assert(0);
  }
}

struct mrc_json_value *
mrc_json_object_entry(struct mrc_json_value *value, unsigned int i)
{
  if (value->type == MRC_JSON_OBJECT) {
    assert(i < value->u.obj.length);
    return value->u.obj.entries[i].value;
  } else if (value->type == MRC_JSON_OBJECT_CLASS) {
    struct mrc_json_class_entry *d = &value->u.cls.descr[i];
    void *base = (void *) value - value->u.cls.off;
    switch (d->value.type) {
    case MRC_JSON_INTEGER:
      d->value.u.integer = * (int *) (base + d->off);
      break;
    case MRC_JSON_DOUBLE:
      d->value.u.dbl = * (double *) (base + d->off);
      break;
    case MRC_JSON_REF_CLASS:
      return * (struct mrc_json_value **) (base + d->off);
    default:
      assert(0);
    }
    return &d->value;
  } else {
    assert(0);
  }
}

const char *
mrc_json_object_entry_name(struct mrc_json_value *value, unsigned int i)
{
  if (value->type == MRC_JSON_OBJECT) {
    assert(i < value->u.obj.length);
    return value->u.obj.entries[i].name;
  } else if (value->type == MRC_JSON_OBJECT_CLASS) {
    struct mrc_json_class_entry *d = &value->u.cls.descr[i];
    return d->name;
  } else {
    assert(0);
  }
}

// ======================================================================
// mrc_json_print

static void
print_indent(int depth)
{
  for (int j = 0; j < depth; j++) {
    printf(" ");
  }
}

static void
mrc_json_print_object(struct mrc_json_value* value, int depth)
{
  assert(value);

  print_indent(depth);
  printf("{\n");

  int length = mrc_json_object_length(value);
  for (int i = 0; i < length; i++) {
    print_indent(depth+2);
    printf("(name) %s : \n", mrc_json_object_entry_name(value, i));
    mrc_json_print(mrc_json_object_entry(value, i), depth+4);
  }

  print_indent(depth);
  printf("}\n");
}

void
mrc_json_print(struct mrc_json_value *value, int depth)
{
  assert(value);

  int type = mrc_json_type(value);
  switch (type) {
  case MRC_JSON_OBJECT:
    mrc_json_print_object(value, depth+1);
    break;
  case MRC_JSON_INTEGER:
    print_indent(depth);
    printf("(int) %d\n", mrc_json_get_integer(value));
    break;
  case MRC_JSON_DOUBLE:
    print_indent(depth);
    printf("(double) %g\n", mrc_json_get_double(value));
    break;
  default:
    fprintf(stderr, "MRC_JSON: unhandled type = %d\n", type);
    assert(0);
  }
}

