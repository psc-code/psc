
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
  return mrc_json_get_integer(mrc_json_get_object_entry(json, name));
}


double
mrc_json_get_object_entry_double(mrc_json_t json, const char *name)
{
  return mrc_json_get_double(mrc_json_get_object_entry(json, name));
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


#if 0

int
mrc_json_get_integer(struct mrc_json_value *value)
{
  assert(value->type == MRC_JSON_INTEGER);
  return value->u.integer;
}

double
mrc_json_get_double(struct mrc_json_value *value)
{
  assert(value->type == MRC_JSON_DOUBLE);
  return value->u.dbl;
}

const char *
mrc_json_get_string(struct mrc_json_value *value)
{
  assert(value->type == MRC_JSON_STRING);
  return value->u.str;
}

int
mrc_json_type(struct mrc_json_value *value)
{
  return value->type & MRC_JSON_TYPE_MASK;
}

// ======================================================================
// mrc_json_object (basic)

static unsigned int
mrc_json_object_basic_length(struct mrc_json_value *value)
{
  return value->u.obj.length;
}

static const char *
mrc_json_object_basic_entry_name(struct mrc_json_value *value, unsigned int i)
{
  assert(i < value->u.obj.length);
  return value->u.obj.entries[i].name;
}

static struct mrc_json_value *
mrc_json_object_basic_entry(struct mrc_json_value *value, unsigned int i)
{
  assert(i < value->u.obj.length);
  return value->u.obj.entries[i].value;
}

// ======================================================================
// mrc_json_object (class)

static unsigned int
mrc_json_object_class_length(struct mrc_json_value *value)
{
  int cnt;
  for (cnt = 0; value->u.cls.descr[cnt].name; cnt++)
    ;
  return cnt;
}

static const char *
mrc_json_object_class_entry_name(struct mrc_json_value *value, unsigned int i)
{
  struct mrc_json_class_entry *d = &value->u.cls.descr[i];
  return d->name;
}

static struct mrc_json_value *
mrc_json_object_class_entry(struct mrc_json_value *value, unsigned int i)
{
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
}

// ======================================================================

static unsigned int mrc_json_object_mrc_obj_length(struct mrc_json_value *value);
static const char *mrc_json_object_mrc_obj_entry_name(struct mrc_json_value *value, unsigned int i);
static struct mrc_json_value *mrc_json_object_mrc_obj_entry(struct mrc_json_value *value, unsigned int i);

unsigned int
mrc_json_object_length(struct mrc_json_value *value)
{
  if (value->type == MRC_JSON_OBJECT) {
    return mrc_json_object_basic_length(value);
  } else if (value->type == MRC_JSON_OBJECT_CLASS) {
    return mrc_json_object_class_length(value);
  } else if (value->type == MRC_JSON_OBJECT_MRC_OBJ) {
    return mrc_json_object_mrc_obj_length(value);
  } else {
    assert(0);
  }
}

const char *
mrc_json_object_entry_name(struct mrc_json_value *value, unsigned int i)
{
  if (value->type == MRC_JSON_OBJECT) {
    return mrc_json_object_basic_entry_name(value, i);
  } else if (value->type == MRC_JSON_OBJECT_CLASS) {
    return mrc_json_object_class_entry_name(value, i);
  } else if (value->type == MRC_JSON_OBJECT_MRC_OBJ) {
    return mrc_json_object_mrc_obj_entry_name(value, i);
  } else {
    assert(0);
  }
}

struct mrc_json_value *
mrc_json_object_entry(struct mrc_json_value *value, unsigned int i)
{
  if (value->type == MRC_JSON_OBJECT) {
    return mrc_json_object_basic_entry(value, i);
  } else if (value->type == MRC_JSON_OBJECT_CLASS) {
    return mrc_json_object_class_entry(value, i);
  } else if (value->type == MRC_JSON_OBJECT_MRC_OBJ) {
    return mrc_json_object_mrc_obj_entry(value, i);
  } else {
    assert(0);
  }
}

unsigned int
mrc_json_array_length(struct mrc_json_value *value)
{
  assert(value->type == MRC_JSON_ARRAY);
  return value->u.arr.length;
}

struct mrc_json_value *
mrc_json_array_entry(struct mrc_json_value *value, unsigned int i)
{
  assert(value->type == MRC_JSON_ARRAY);
  assert(i < value->u.arr.length);
  return value->u.arr.entries[i];
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

static void
mrc_json_print_array(struct mrc_json_value* value, int depth)
{
  assert(value);

  print_indent(depth);
  printf("[\n");

  int length = mrc_json_array_length(value);
  for (int i = 0; i < length; i++) {
    mrc_json_print(mrc_json_array_entry(value, i), depth+2);
  }

  print_indent(depth);
  printf("]\n");
}

void
mrc_json_print(struct mrc_json_value *value, int depth)
{
  assert(value);

  int type = mrc_json_type(value);
  switch (type) {
  case MRC_JSON_NONE:
    print_indent(depth);
    printf("(none)\n");
    break;
  case MRC_JSON_OBJECT:
    mrc_json_print_object(value, depth+1);
    break;
  case MRC_JSON_ARRAY:
    mrc_json_print_array(value, depth+1);
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
  default:
    fprintf(stderr, "MRC_JSON: unhandled type = %d\n", type);
    assert(0);
  }
}

// ======================================================================
// mrc_json_object (mrc_obj)

#include <mrc_obj.h>

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

static struct mrc_json_value *
mrc_descr_entry(struct param *param, char *p)
{
  // FIXME, this mrc_json_value is per-class, not per-object,
  // so if we look at multiple objects of the same class at the same
  // time, bugs will happen :(
  struct mrc_json_value *v = &param->json;
  
  switch (param->type) {
  case MRC_VAR_OBJ:
    v->type = MRC_JSON_NONE; // FIXME
    break;
  case PT_SELECT:
    v->type = MRC_JSON_NONE; // FIXME
    break;
  case PT_INT3: {
    int *int3 = (int *) p;
    static struct mrc_json_value entry[3];
    static struct mrc_json_value *entries[3] = {
      &entry[0], &entry[1], &entry[2],
    };
    for (int d = 0; d < 3; d++) {
      entry[d].type = MRC_JSON_INTEGER;
      entry[d].u.integer = int3[d];
    }
    v->type = MRC_JSON_ARRAY; // FIXME
    v->u.arr.length = 3;
    v->u.arr.entries = entries;
    break;
  }
  default:
    fprintf(stderr, "unhandled type: %d\n", param->type);
    assert(0);
  }

  return v;
}

static unsigned int
mrc_json_object_mrc_obj_length(struct mrc_json_value *value)
{
  struct mrc_obj *obj = container_of(value, struct mrc_obj, json);
  struct mrc_class *cls = obj->cls;

  int cnt = 2; // mrc_obj type and name

  if (obj->ops) { // mrc_obj subtype
    cnt++;
  }
  
  cnt += mrc_descr_length(cls->param_descr);

  if (obj->ops) {
    cnt += mrc_descr_length(obj->ops->param_descr);
  }
  
  printf("length = %d\n", cnt);
  return cnt;
}

static const char *
mrc_json_object_mrc_obj_entry_name(struct mrc_json_value *value, unsigned int i)
{
  struct mrc_obj *obj = container_of(value, struct mrc_obj, json);
  struct mrc_class *cls = obj->cls;

  if (i == 0) {
    return "mrc_obj_type";
  }
  i--;

  if (i == 0) {
    return "mrc_obj_name";
  }
  i--;

  if (obj->ops) {
    if (i == 0) {
      return "mrc_obj_subtype";
    }
    i--;
  }
  
  int len = mrc_descr_length(cls->param_descr);
  if (i < len) {
    return cls->param_descr[i].name;
  }
  i -= len;

  if (obj->ops) {
    int len = mrc_descr_length(obj->ops->param_descr);
    if (i < len) {
      return obj->ops->param_descr[i].name;
    }
    i-= len;
  }

  assert(0);
}

static struct mrc_json_value *
mrc_json_object_mrc_obj_entry(struct mrc_json_value *value, unsigned int i)
{
  struct mrc_obj *obj = container_of(value, struct mrc_obj, json);
  struct mrc_class *cls = obj->cls;

  if (i == 0) {
    static struct mrc_json_value v_type;
    v_type.type = MRC_JSON_STRING;
    v_type.u.str = obj->cls->name;
    return &v_type;
  }
  i--;

  if (i == 0) {
    static struct mrc_json_value v_name;
    v_name.type = MRC_JSON_STRING;
    v_name.u.str = obj->name;
    return &v_name;
  }
  i--;
  
  if (obj->ops) {
    if (i == 0) {
      static struct mrc_json_value v_name;
      v_name.type = MRC_JSON_STRING;
      v_name.u.str = obj->ops->name;
      return &v_name;
    }
    i--;
  }
  
  int len = mrc_descr_length(cls->param_descr);
  if (i < len) {
    return mrc_descr_entry(&cls->param_descr[i], (char *) obj + cls->param_offset);
  }
  i -= len;

  if (obj->ops) {
    int len = mrc_descr_length(obj->ops->param_descr);
    if (i < len) {
      return mrc_descr_entry(&obj->ops->param_descr[i], (char *) obj->subctx + obj->ops->param_offset);
    }
    i-= len;
  }

  assert(0);

}

#endif
