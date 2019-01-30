
#ifndef MRC_JSON_H
#define MRC_JSON_H

#include <unistd.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif
#if 0 // hack to fix indentation
}
#endif

enum mrc_json_type {
  // FIXME? This have to match the json.h types, but I don't to want to include that here
  MRC_JSON_NONE,
  MRC_JSON_OBJECT,
  MRC_JSON_ARRAY,
  MRC_JSON_INTEGER,
  MRC_JSON_DOUBLE,
  MRC_JSON_STRING,
  MRC_JSON_BOOLEAN,
  MRC_JSON_NULL,
};

typedef struct mrc_json mrc_json_t;

struct mrc_json_ops {
  int (*get_type)(mrc_json_t json);
  int (*get_integer)(mrc_json_t json);
  double (*get_double)(mrc_json_t json);
  const char * (*get_string)(mrc_json_t json);
  bool (*get_boolean)(mrc_json_t json);
  unsigned int (*get_object_length)(mrc_json_t json);
  const char * (*get_object_entry_name)(mrc_json_t json, unsigned int i);
  mrc_json_t (*get_object_entry_value)(mrc_json_t json, unsigned int i);
  unsigned int (*get_array_length)(mrc_json_t json);
  mrc_json_t (*get_array_entry)(mrc_json_t json, unsigned int i);
};

struct mrc_json_parser {
  void *value;
};

struct mrc_json_mrc_obj {
  union {
    struct mrc_obj *obj;
    int integer;
    double dbl;
    bool boolean;
    const char *string;
    int *int3_ptr;
    double *dbl3_ptr;
  } v;
  int type;
};

struct mrc_json {
  union {
    struct mrc_json_parser parser;
    struct mrc_json_mrc_obj mrc;
  } u;
  struct mrc_json_ops *ops;
};

int mrc_json_get_type(mrc_json_t value);

int mrc_json_get_integer(mrc_json_t value);
double mrc_json_get_double(mrc_json_t value);
const char *mrc_json_get_string(mrc_json_t value);
bool mrc_json_get_boolean(mrc_json_t value);

void mrc_json_get_int3(mrc_json_t json, int arr[3]);
void mrc_json_get_double3(mrc_json_t json, double arr[3]);
void mrc_json_get_float3(mrc_json_t json, float arr[3]); // will convert double -> float

unsigned int mrc_json_get_object_length(mrc_json_t value);
const char *mrc_json_get_object_entry_name(mrc_json_t value, unsigned int i);
mrc_json_t mrc_json_get_object_entry_value(mrc_json_t value, unsigned int i);
mrc_json_t mrc_json_get_object_entry(mrc_json_t value, const char *name);
int mrc_json_get_object_entry_integer(mrc_json_t value, const char *name);
double mrc_json_get_object_entry_double(mrc_json_t value, const char *name);
void mrc_json_get_object_entry_int3(mrc_json_t value, const char *name, int arr[3]);
void mrc_json_get_object_entry_double3(mrc_json_t value, const char *name, double arr[3]);

unsigned int mrc_json_get_array_length(mrc_json_t value);
mrc_json_t mrc_json_get_array_entry(mrc_json_t value, unsigned int i);
int mrc_json_get_array_entry_integer(mrc_json_t value, unsigned int i);
double mrc_json_get_array_entry_double(mrc_json_t value, unsigned int i);

mrc_json_t mrc_json_object_new(unsigned int length);
mrc_json_t mrc_json_array_new(unsigned int length);
mrc_json_t mrc_json_integer_new(int integer);
mrc_json_t mrc_json_double_new(double dbl);
mrc_json_t mrc_json_string_new(const char *str);
mrc_json_t mrc_json_boolean_new(bool boolean);
mrc_json_t mrc_json_integer_array_new(unsigned int length, int *arr);
mrc_json_t mrc_json_double_array_new(unsigned int length, double *arr);
void mrc_json_object_push(mrc_json_t obj, const char *name, mrc_json_t entry);
void mrc_json_object_push_integer(mrc_json_t obj, const char *name, int integer);
void mrc_json_object_push_double(mrc_json_t obj, const char *name, double dbl);
void mrc_json_object_push_boolean(mrc_json_t obj, const char *name, bool boolean);
void mrc_json_object_push_integer_array(mrc_json_t obj, const char *name, unsigned int length, int *arr);
void mrc_json_object_push_double_array(mrc_json_t obj, const char *name, unsigned int length, double *arr);
void mrc_json_array_push(mrc_json_t arr, mrc_json_t entry);
void mrc_json_array_push_integer(mrc_json_t arr, int integer);
void mrc_json_array_push_double(mrc_json_t arr, double dbl);
void mrc_json_array_push_integer_array(mrc_json_t arr, unsigned int length, int *int_arr);
void mrc_json_array_push_double_array(mrc_json_t arr, unsigned int length, double *dbl_arr);

void mrc_json_print(mrc_json_t value, unsigned int depth);
char *mrc_json_to_string(mrc_json_t json);

// parse from actual JSON string
mrc_json_t mrc_json_parse(const char *buf);

// create mrc_json_t wrapper from json_value *
mrc_json_t mrc_json_from_json_parser(void *value);

// create mrc_json_t wrapper from mrc_obj
mrc_json_t mrc_obj_to_json(struct mrc_obj *obj);
// FIXME, maybe instead of this ugly macro we could add something to mrc_obj.h?
#define MRC_OBJ_TO_JSON(obj) mrc_obj_to_json((struct mrc_obj *) obj)

#if 0
{
#endif
#ifdef __cplusplus
}
#endif

#endif
