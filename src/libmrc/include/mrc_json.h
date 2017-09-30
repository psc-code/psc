
#ifndef MRC_JSON_H
#define MRC_JSON_H

#include <unistd.h>

enum mrc_json_type {
  MRC_JSON_NONE,
  MRC_JSON_INTEGER,
  MRC_JSON_DOUBLE,
  MRC_JSON_OBJECT,

  MRC_JSON_REF,

  MRC_JSON_TYPE_MASK = 0xff,
  MRC_JSON_FLAG_CLASS = 0x100,
  MRC_JSON_OBJECT_CLASS = MRC_JSON_OBJECT | MRC_JSON_FLAG_CLASS,
  MRC_JSON_REF_CLASS = MRC_JSON_REF | MRC_JSON_FLAG_CLASS,
};

struct mrc_json_value {
  int type;
  union {
    int integer;

    double dbl;

    struct {
      struct mrc_json_object_entry *entries;
      unsigned int length;
    } obj;

    struct {
      ssize_t off;
      struct mrc_json_class_entry *descr;
    } cls;
  } u;
};

struct mrc_json_object_entry {
  const char *name;
  struct mrc_json_value *value;
};

struct mrc_json_class_entry {
  const char *name;
  ssize_t off;
  struct mrc_json_value value;
};

#ifndef offsetof
#define offsetof(TYPE, MEMBER) ((size_t) &((TYPE *)0)->MEMBER)
#endif

void mrc_json_print(struct mrc_json_value* value, int depth);



#endif
