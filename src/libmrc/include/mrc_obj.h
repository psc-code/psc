
#ifndef MRC_OBJ_H
#define MRC_OBJ_H

#include <mrc_common.h>
#include <mrc_list.h>
#include <mrc_params.h>

#include <mpi.h>
#include <stdbool.h>

BEGIN_C_DECLS

struct mrc_io;

struct mrc_obj {
  MPI_Comm comm;
  char *name;
  struct mrc_obj_ops *ops;
  struct mrc_class *cls;
  void *subctx;
  int refcount;
  list_t instance_entry; //< keep track of all instances of mrc_obj's
  list_t child_entry; //< an mrc_obj can be child of exactly one parent mrc_obj
  list_t children_list; //< this is the list where a parent keeps track of children
  list_t dict_list; //< keep track of additional, dynamic, properties
  bool view_flag; //< if true, call ::view() at the end of ::setup()
  bool is_setup; //< keep track of whether ::setup() was already called
  void *creation_trace;
};

struct mrc_dict_entry {
  struct param prm;
  union param_u val;

  list_t entry; //< entry in obj::dict_list
};

typedef void (*mrc_void_func_t)(void);

struct mrc_obj_method {
  const char *name;
  mrc_void_func_t func;
};

#define MRC_OBJ_METHOD(name, func)              \
  { name, (mrc_void_func_t) func }

#define MRC_SUBCLASS_OPS(obj_type)                      \
  list_t list;                                          \
  const char *name;                                     \
  size_t size;                                          \
  struct param *param_descr;                            \
  size_t param_offset;                                  \
  struct mrc_obj_method *methods;                       \
  void (*create)(obj_type *);                           \
  void (*destroy)(obj_type *);                          \
  void (*set_from_options)(obj_type *);                 \
  void (*view)(obj_type *);                             \
  void (*setup)(obj_type *);                            \
  void (*write)(obj_type *, struct mrc_io *);           \
  void (*read)(obj_type *, struct mrc_io *)

struct mrc_obj_ops {
  MRC_SUBCLASS_OPS(struct mrc_obj);
};

#define DECLARE_STRUCT_MRC_CLASS(sfx, obj_type)         \
  struct mrc_class ## sfx {                             \
    const char *name;                                   \
    list_t subclasses;                                  \
    list_t instances;                                   \
    int nr_instances;         \
    int nr_instances_last_printed; \
    int nr_creates_last_printed;   \
    int nr_creates;           \
    int nr_frees;          \
    list_t active_classes;    \
    bool watch;          \
    size_t size;                                        \
    struct param *param_descr;                          \
    struct mrc_obj_method *methods;                     \
    size_t param_offset;                                \
    void (*init)(void);                                 \
    bool initialized;                                   \
    void (*create)(obj_type *);                         \
    void (*destroy)(obj_type *);                        \
    void (*set_from_options)(obj_type *);               \
    void (*view)(obj_type *);                           \
    void (*setup)(obj_type *);                          \
    void (*read)(obj_type *, struct mrc_io *);          \
    void (*write)(obj_type *, struct mrc_io *);         \
  }

DECLARE_STRUCT_MRC_CLASS(, struct mrc_obj);

struct mrc_io;

struct mrc_obj *mrc_obj_create(MPI_Comm comm);
struct mrc_obj *__mrc_obj_create(MPI_Comm comm, struct mrc_class *cls);
struct mrc_obj *mrc_obj_get(struct mrc_obj *obj);
void mrc_obj_put(struct mrc_obj *obj);
void mrc_obj_destroy(struct mrc_obj *obj);
MPI_Comm mrc_obj_comm(struct mrc_obj *obj);
const char *mrc_obj_type(struct mrc_obj *obj);
const char *mrc_obj_name(struct mrc_obj *obj);
void mrc_obj_set_name(struct mrc_obj *obj, const char *name);
void mrc_obj_set_type(struct mrc_obj *obj, const char *type);
void mrc_obj_set_from_options(struct mrc_obj *obj);
void mrc_obj_set_param_int(struct mrc_obj *obj, const char *name, int val);
void mrc_obj_set_param_float(struct mrc_obj *obj, const char *name, float val);
void mrc_obj_set_param_double(struct mrc_obj *obj, const char *name, double val);
void mrc_obj_set_param_string(struct mrc_obj *obj, const char *name, const char *val);
void mrc_obj_set_param_select(struct mrc_obj *obj, const char *name, int val);
void mrc_obj_set_param_bool(struct mrc_obj *obj, const char *name, bool val);
void mrc_obj_set_param_int3(struct mrc_obj *obj, const char *name, const int val[3]);
void mrc_obj_set_param_float3(struct mrc_obj *obj, const char *name, const float val[3]);
void mrc_obj_set_param_double3(struct mrc_obj *obj, const char *name, const double val[3]);
void mrc_obj_set_param_int_array(struct mrc_obj *obj, const char *name, int nr_vals, const int val[]);
void mrc_obj_set_param_float_array(struct mrc_obj *obj, const char *name, int nr_vals, const float val[]);
void mrc_obj_set_param_ptr(struct mrc_obj *obj, const char *name, void* val);
void mrc_obj_set_param_obj(struct mrc_obj *obj, const char *name, void* val);
int mrc_obj_get_param_bool(struct mrc_obj *obj, const char *name, bool *pval);
int mrc_obj_get_param_int(struct mrc_obj *obj, const char *name, int *pval);
int mrc_obj_get_param_float(struct mrc_obj *obj, const char *name, float *pval);
int mrc_obj_get_param_double(struct mrc_obj *obj, const char *name, double *pval);
int mrc_obj_get_param_string(struct mrc_obj *obj, const char *name, const char **val);
int mrc_obj_get_param_int3(struct mrc_obj *obj, const char *name, int *pval);
int mrc_obj_get_param_float3(struct mrc_obj *obj, const char *name, float *pval);
int mrc_obj_get_param_double3(struct mrc_obj *obj, const char *name, double *pval);
int mrc_obj_get_param_obj(struct mrc_obj *obj, const char *name, struct mrc_obj **pval);
int mrc_obj_get_param_ptr(struct mrc_obj *obj, const char *name, void **pval);
int mrc_obj_get_param_float_array_nr_vals(struct mrc_obj *obj, const char *name, int *nr_vals);
int mrc_obj_get_param_float_array(struct mrc_obj *obj, const char *name, float *pval);

int mrc_obj_get_var(struct mrc_obj *obj, const char *name, union param_u **pv);
int mrc_obj_get_var_double(struct mrc_obj *obj, const char *name, double *pval);
struct mrc_obj *mrc_obj_get_var_obj(struct mrc_obj *obj, const char *name);

void mrc_obj_dict_add(struct mrc_obj *obj, int type, const char *name,
                      union param_u *pv);
void mrc_obj_dict_add_int(struct mrc_obj *obj, const char *name, int val);
void mrc_obj_dict_add_bool(struct mrc_obj *obj, const char *name, bool val);
void mrc_obj_dict_add_float(struct mrc_obj *obj, const char *name, float val);
void mrc_obj_dict_add_double(struct mrc_obj *obj, const char *name, double val);
void mrc_obj_dict_add_string(struct mrc_obj *obj, const char *name, const char *val);
void mrc_obj_dict_add_obj(struct mrc_obj *obj, const char *name, struct mrc_obj *val);
void mrc_obj_dict_add_float_array(struct mrc_obj *obj, const char *name, float *vals, int nr_vals);

void mrc_obj_view(struct mrc_obj *obj);
void mrc_obj_setup(struct mrc_obj *obj);
void mrc_obj_setup_super(struct mrc_obj *obj);
void mrc_obj_setup_member_objs(struct mrc_obj *obj);
void mrc_obj_setup_member_objs_sub(struct mrc_obj *obj);
void mrc_obj_setup_children(struct mrc_obj *obj);
void mrc_obj_add_child(struct mrc_obj *obj, struct mrc_obj *child);
struct mrc_obj *mrc_obj_find_child(struct mrc_obj *obj, const char *name);
void mrc_obj_write(struct mrc_obj *obj, struct mrc_io *io);
struct mrc_obj *mrc_obj_read(struct mrc_io *io, const char *path, struct mrc_class *cls);
struct mrc_obj *mrc_obj_read_comm(struct mrc_io *io, const char *path, struct mrc_class *cls,
				  MPI_Comm comm);
void mrc_obj_read_super(struct mrc_obj *obj, struct mrc_io *io);
void mrc_obj_read_member_objs(struct mrc_obj *obj, struct mrc_io *io);
void mrc_obj_read_member_objs_sub(struct mrc_obj *obj, struct mrc_io *io);
void mrc_obj_read_children(struct mrc_obj *obj, struct mrc_io *io);
bool mrc_obj_is_setup(struct mrc_obj *obj);
mrc_void_func_t mrc_obj_get_method(struct mrc_obj *obj, const char *name);

struct mrc_obj_ops *mrc_obj_find_subclass_ops(struct mrc_class *cls, const char *subclass);

enum {
  CLASS_INFO_VERB_NONE = 0,
  CLASS_INFO_VERB_DIFF,
  CLASS_INFO_VERB_ACTIVE,
  CLASS_INFO_VERB_FULL
};
int mrc_obj_print_class_info(int verbosity);

#ifdef __cplusplus
#define MRC_CLASS_DECLARE_1(pfx, obj_type)				\
  obj_type;                                                             \
  DECLARE_STRUCT_MRC_CLASS(_ ## pfx, obj_type);				\
  extern struct mrc_class_ ##pfx##_ mrc_class_ ## pfx;			\
  
#else
#define MRC_CLASS_DECLARE_1(pfx, obj_type)				\
  obj_type;                                                             \
  DECLARE_STRUCT_MRC_CLASS(_ ## pfx, obj_type);				\
  extern struct mrc_class_ ##pfx mrc_class_ ## pfx;			\
  
#endif

#define MRC_CLASS_DECLARE(pfx, obj_type)                                \
  MRC_CLASS_DECLARE_1(pfx, obj_type)					\
  static inline obj_type *                                              \
  pfx ## _create(MPI_Comm comm)                                         \
  {                                                                     \
    return (obj_type *)                                                 \
      __mrc_obj_create(comm, (struct mrc_class *) &mrc_class_ ## pfx);  \
  }                                                                     \
                                                                        \
  static inline struct mrc_obj *                                        \
  pfx ## _to_mrc_obj(obj_type *obj)                                     \
  {                                                                     \
    return (struct mrc_obj *) obj;                                      \
  }                                                                     \
                                                                        \
  static inline obj_type *                                              \
  pfx ## _get(obj_type *obj)                                            \
  {                                                                     \
    return (obj_type *)                                                 \
      mrc_obj_get((struct mrc_obj *)obj);                               \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _destroy(obj_type *obj)                                        \
  {                                                                     \
    mrc_obj_destroy((struct mrc_obj *)obj);                             \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _put(obj_type *obj)                                            \
  {                                                                     \
    mrc_obj_put((struct mrc_obj *)obj);                                 \
  }                                                                     \
                                                                        \
  static inline MPI_Comm                                                \
  pfx ## _comm(obj_type *obj)                                           \
  {                                                                     \
    return mrc_obj_comm((struct mrc_obj *)obj);                         \
  }                                                                     \
                                                                        \
  static inline const char *                                            \
  pfx ## _name(obj_type *obj)                                           \
  {                                                                     \
    return mrc_obj_name((struct mrc_obj *)obj);                         \
  }                                                                     \
                                                                        \
  static inline const char *                                            \
  pfx ## _type(obj_type *obj)                                           \
  {                                                                     \
    return mrc_obj_type((struct mrc_obj *)obj);                         \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_name(obj_type *obj, const char *name)                     \
  {                                                                     \
    mrc_obj_set_name((struct mrc_obj *)obj, name);                      \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_type(obj_type *obj, const char *name)                     \
  {                                                                     \
    mrc_obj_set_type((struct mrc_obj *)obj, name);                      \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_int(obj_type *obj, const char *name, int val)       \
  {                                                                     \
    mrc_obj_set_param_int((struct mrc_obj *)obj, name, val);            \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_float(obj_type *obj, const char *name, float val)   \
  {                                                                     \
    mrc_obj_set_param_float((struct mrc_obj *)obj, name, val);          \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_double(obj_type *obj, const char *name, double val) \
  {                                                                     \
    mrc_obj_set_param_double((struct mrc_obj *)obj, name, val);         \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_string(obj_type *obj, const char *name, const char *val) \
  {                                                                     \
    mrc_obj_set_param_string((struct mrc_obj *)obj, name, val);         \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_select(obj_type *obj, const char *name, int val)    \
  {                                                                     \
    mrc_obj_set_param_select((struct mrc_obj *)obj, name, val);         \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_bool(obj_type *obj, const char *name, bool val)     \
  {                                                                     \
    mrc_obj_set_param_bool((struct mrc_obj *)obj, name, val);           \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_int3(obj_type *obj, const char *name, const int val[3])     \
  {                                                                     \
    mrc_obj_set_param_int3((struct mrc_obj *)obj, name, val);           \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_float3(obj_type *obj, const char *name, const float val[3]) \
  {                                                                     \
    mrc_obj_set_param_float3((struct mrc_obj *)obj, name, val);         \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_double3(obj_type *obj, const char *name, const double val[3])       \
  {                                                                     \
    mrc_obj_set_param_double3((struct mrc_obj *)obj, name, val);        \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_int_array(obj_type *obj, const char *name,          \
                              int nr_vals, const int val[])             \
  {                                                                     \
    mrc_obj_set_param_int_array((struct mrc_obj *)obj, name, nr_vals, val); \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_float_array(obj_type *obj, const char *name,	\
				int nr_vals, const float val[])		\
  {                                                                     \
    mrc_obj_set_param_float_array((struct mrc_obj *)obj, name, nr_vals, val); \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_ptr(obj_type *obj, const char *name, void* val)     \
  {                                                                     \
    mrc_obj_set_param_ptr((struct mrc_obj *)obj, name, val);            \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_param_obj(obj_type *obj, const char *name, void* val)     \
  {                                                                     \
    mrc_obj_set_param_obj((struct mrc_obj *)obj, name, val);            \
  }                                                                     \
                                                                        \
  static inline int                                                     \
  pfx ## _get_param_bool(obj_type *obj, const char *name, bool *pval)   \
  {                                                                     \
    return mrc_obj_get_param_bool((struct mrc_obj *)obj, name, pval);   \
  }                                                                     \
                                                                        \
  static inline int                                                     \
  pfx ## _get_param_int(obj_type *obj, const char *name, int *pval)     \
  {                                                                     \
    return mrc_obj_get_param_int((struct mrc_obj *)obj, name, pval);    \
  }                                                                     \
                                                                        \
  static inline int                                                     \
  pfx ## _get_param_float(obj_type *obj, const char *name, float *pval) \
  {                                                                     \
    return mrc_obj_get_param_float((struct mrc_obj *)obj, name, pval);  \
  }                                                                     \
                                                                        \
  static inline int                                                     \
  pfx ## _get_param_double(obj_type *obj, const char *name, double *pval) \
  {                                                                     \
    return mrc_obj_get_param_double((struct mrc_obj *)obj, name, pval); \
  }                                                                     \
                                                                        \
  static inline int                                                     \
  pfx ## _get_param_string(obj_type *obj, const char *name, const char **val) \
  {                                                                     \
    return mrc_obj_get_param_string((struct mrc_obj *)obj, name, val);  \
  }                                                                     \
                                                                        \
  static inline int                                                     \
  pfx ## _get_param_int3(obj_type *obj, const char *name, int *pval)    \
  {                                                                     \
    return mrc_obj_get_param_int3((struct mrc_obj *)obj, name, pval);   \
  }                                                                     \
                                                                        \
  static inline int                                                     \
  pfx ## _get_param_float3(obj_type *obj, const char *name, float *pval) \
  {                                                                     \
    return mrc_obj_get_param_float3((struct mrc_obj *)obj, name, pval); \
  }                                                                     \
                                                                        \
  static inline int                                                     \
  pfx ## _get_param_double3(obj_type *obj, const char *name, double *pval) \
  {                                                                     \
    return mrc_obj_get_param_double3((struct mrc_obj *)obj, name, pval); \
  }                                                                     \
                                                                        \
  static inline int                                                     \
  pfx ## _get_param_obj(obj_type *obj, const char *name, void *pval)    \
  {                                                                     \
    return mrc_obj_get_param_obj((struct mrc_obj *)obj, name,           \
                                 (struct mrc_obj **) pval);             \
  }                                                                     \
                                                                        \
  static inline int                                                     \
  pfx ## _get_param_float_array_nr_vals(obj_type *obj, const char *name,\
      int *nr_vals)                                                     \
  {                                                                     \
    return mrc_obj_get_param_float_array_nr_vals((struct mrc_obj *)obj, \
        name, nr_vals);                                                 \
  }                                                                     \
                                                                        \
  static inline int                                                     \
  pfx ## _get_param_float_array(obj_type *obj, const char *name,        \
      float *vals)                                                      \
  {                                                                     \
    return mrc_obj_get_param_float_array((struct mrc_obj *)obj, name, vals); \
  }                                                                     \
  static inline int                                                    \
  pfx ## _get_param_ptr(obj_type *obj, const char *name, void **pval)  \
  {                                                                    \
    return mrc_obj_get_param_ptr((struct mrc_obj *)obj, name, pval);    \
  }                                                                     \
  static inline int                                                     \
  pfx ## _get_var(obj_type *obj, const char *name, union param_u **pv)  \
  {                                                                     \
    return mrc_obj_get_var((struct mrc_obj *)obj, name, pv);            \
  }                                                                     \
                                                                        \
  static inline int							\
  pfx ## _get_var_double(obj_type *obj, const char *name, double *pval)	\
  {                                                                     \
    return mrc_obj_get_var_double((struct mrc_obj *)obj, name, pval);	\
  }                                                                     \
                                                                        \
  /* we should return a subtype of mrc_obj, but that's not */           \
  /* easily possible without having the user do the casting */          \
  static inline void *                                                  \
  pfx ## _get_var_obj(obj_type *obj, const char *name)                  \
  {                                                                     \
    return mrc_obj_get_var_obj((struct mrc_obj *)obj, name);            \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _dict_add(obj_type *obj, int type, const char *name,           \
                   union param_u *pv)                                   \
  {                                                                     \
    mrc_obj_dict_add((struct mrc_obj *)obj, type, name, pv);            \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _dict_add_int(obj_type *obj, const char *name, int val)        \
  {                                                                     \
    mrc_obj_dict_add_int((struct mrc_obj *)obj, name, val);             \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _dict_add_bool(obj_type *obj, const char *name, bool val)      \
  {                                                                     \
    mrc_obj_dict_add_bool((struct mrc_obj *)obj, name, val);            \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _dict_add_float(obj_type *obj, const char *name, float val)    \
  {                                                                     \
    mrc_obj_dict_add_float((struct mrc_obj *)obj, name, val);           \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _dict_add_double(obj_type *obj, const char *name, double val)  \
  {                                                                     \
    mrc_obj_dict_add_double((struct mrc_obj *)obj, name, val);          \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _dict_add_string(obj_type *obj, const char *name, const char *val) \
  {                                                                     \
    mrc_obj_dict_add_string((struct mrc_obj *)obj, name, val);          \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _dict_add_obj(obj_type *obj, const char *name, void *val)      \
  {                                                                     \
    mrc_obj_dict_add_obj((struct mrc_obj *)obj, name,                   \
                         (struct mrc_obj *)val);                        \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _dict_add_float_array(obj_type *obj, const char *name,         \
      float *vals, int nr_vals)                                         \
  {                                                                     \
    mrc_obj_dict_add_float_array((struct mrc_obj *)obj, name,           \
                                  vals, nr_vals);                       \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _set_from_options(obj_type *obj)                               \
  {                                                                     \
    mrc_obj_set_from_options((struct mrc_obj *)obj);                    \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _view(obj_type *obj)                                           \
  {                                                                     \
    mrc_obj_view((struct mrc_obj *)obj);                                \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _setup(obj_type *obj)                                          \
  {                                                                     \
    mrc_obj_setup((struct mrc_obj *)obj);                               \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _setup_super(obj_type *obj)                                    \
  {                                                                     \
    mrc_obj_setup_super((struct mrc_obj *)obj);                         \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _setup_member_objs(obj_type *obj)                              \
  {                                                                     \
    mrc_obj_setup_member_objs((struct mrc_obj *)obj);                   \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _setup_member_objs_sub(obj_type *obj)                          \
  {                                                                     \
    mrc_obj_setup_member_objs_sub((struct mrc_obj *)obj);               \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _setup_children(obj_type *obj)                                 \
  {                                                                     \
    mrc_obj_setup_children((struct mrc_obj *)obj);                      \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _add_child(obj_type *obj, struct mrc_obj *child)               \
  {                                                                     \
    mrc_obj_add_child((struct mrc_obj *)obj, child);                    \
  }                                                                     \
                                                                        \
  static inline struct mrc_obj *                                        \
  pfx ## _find_child(obj_type *obj, const char *name)                   \
  {                                                                     \
    return mrc_obj_find_child((struct mrc_obj *)obj, name);             \
  }                                                                     \
                                                                        \
  static inline obj_type *                                              \
  pfx ## _read(struct mrc_io *io, const char *path)                     \
  {                                                                     \
    return (obj_type *) mrc_obj_read(io, path,                          \
                                     (struct mrc_class *) &mrc_class_ ## pfx);\
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _read_super(obj_type *obj, struct mrc_io *io)                  \
  {                                                                     \
    mrc_obj_read_super((struct mrc_obj *)obj, io);                      \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _read_children(obj_type *obj, struct mrc_io *io)               \
  {                                                                     \
    mrc_obj_read_children((struct mrc_obj *) obj, io);                  \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _read_member_objs(obj_type *obj, struct mrc_io *io)            \
  {                                                                     \
    mrc_obj_read_member_objs((struct mrc_obj *) obj, io);               \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _read_member_objs_sub(obj_type *obj, struct mrc_io *io)        \
  {                                                                     \
    mrc_obj_read_member_objs_sub((struct mrc_obj *) obj, io);           \
  }                                                                     \
                                                                        \
  static inline void                                                    \
  pfx ## _write(obj_type *obj, struct mrc_io *io)                       \
  {                                                                     \
    mrc_obj_write((struct mrc_obj *)obj, io);                           \
  }                                                                     \
                                                                        \
  static inline bool                                                    \
  pfx ## _is_setup(obj_type *obj)                                       \
  {                                                                     \
    return mrc_obj_is_setup((struct mrc_obj *)obj);                     \
  }                                                                     \
                                                                        \
  static inline mrc_void_func_t                                         \
  pfx ## _get_method(obj_type *obj, const char *name)                   \
  {                                                                     \
    return mrc_obj_get_method((struct mrc_obj *)obj, name);             \
  }                                                                     \
                                                                        \
  /* force a semicolon following the macro use */                       \
  struct __dummy


// for direct instantiation of (non-derived) mrc_obj
extern struct mrc_class mrc_class_mrc_obj;

// use a macro here to do the casting to mrc_obj_ops

#define mrc_class_register_subclass(cls, ops) \
  __mrc_class_register_subclass((struct mrc_class *)(cls), (struct mrc_obj_ops *)(ops))

void __mrc_class_register_subclass(struct mrc_class *cls,
                                   struct mrc_obj_ops *ops);

#define mrc_to_subobj(o, subobj_type) ((subobj_type *)((o)->obj.subctx))

#define mrc_obj_for_each_child(item, parent, type) \
  __list_for_each_entry(item, &parent->obj.children_list, obj.child_entry, type)

END_C_DECLS

#endif
