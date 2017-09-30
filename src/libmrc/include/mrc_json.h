
#ifndef MRC_JSON_H
#define MRC_JSON_H

#include <unistd.h>

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

#include "json.h"

typedef struct {
  json_value *value;
} mrc_json_t;

int mrc_json_get_type(mrc_json_t value);

int mrc_json_get_integer(mrc_json_t value);
int mrc_json_get_double(mrc_json_t value);
const char *mrc_json_get_string(mrc_json_t value);

unsigned int mrc_json_get_object_length(mrc_json_t value);
const char *mrc_json_get_object_entry_name(mrc_json_t value, unsigned int i);
mrc_json_t mrc_json_get_object_entry_value(mrc_json_t value, unsigned int i);

unsigned int mrc_json_get_array_length(mrc_json_t value);
mrc_json_t mrc_json_get_array_entry(mrc_json_t value, unsigned int i);

void mrc_json_print(mrc_json_t value, unsigned int depth);
char *mrc_json_to_string(mrc_json_t json);

#endif
