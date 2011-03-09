
#ifndef MRC_OBJ_H
#define MRC_OBJ_H

#include <mrc_common.h>
#include <mrc_list.h>

#include <mpi.h>
#include <stdbool.h>

struct mrc_io;

struct mrc_obj {
  MPI_Comm comm;
  char *name;
  struct mrc_obj_ops *ops;
  struct mrc_class *class;
  void *subctx;
  int refcount;
  list_t instance_entry;
};

#define MRC_OBJ_OPS				 \
  list_t list;					 \
  const char *name;				 \
  size_t size;					 \
  struct param *param_descr;			 \
  size_t param_offset;				 \
  void (*create)(struct mrc_obj *);		 \
  void (*destroy)(struct mrc_obj *);		 \
  void (*set_from_options)(struct mrc_obj *);	 \
  void (*view)(struct mrc_obj *);		 \
  void (*setup)(struct mrc_obj *);		 \
  void (*write)(struct mrc_obj *, struct mrc_io *)

struct mrc_obj_ops {
  MRC_OBJ_OPS;
};

struct mrc_class {
  const char *name;
  list_t *subclasses;
  list_t instances;
  size_t size;
  struct param *param_descr;
  size_t param_offset;
  void (*init)(void);
  bool initialized;
  void (*create)(struct mrc_obj *);
  void (*destroy)(struct mrc_obj *);
  void (*set_from_options)(struct mrc_obj *);
  void (*view)(struct mrc_obj *);
  void (*setup)(struct mrc_obj *);
  void (*read)(struct mrc_obj *, struct mrc_io *);
  void (*write)(struct mrc_obj *, struct mrc_io *);
};

struct mrc_io;

struct mrc_obj *mrc_obj_create(MPI_Comm comm, struct mrc_class *class);
struct mrc_obj *mrc_obj_get(struct mrc_obj *obj);
void mrc_obj_put(struct mrc_obj *obj);
void mrc_obj_destroy(struct mrc_obj *obj);
MPI_Comm mrc_obj_comm(struct mrc_obj *obj);
const char *mrc_obj_name(struct mrc_obj *obj);
void mrc_obj_set_name(struct mrc_obj *obj, const char *name);
void mrc_obj_set_type(struct mrc_obj *obj, const char *type);
void mrc_obj_set_from_options(struct mrc_obj *obj);
void mrc_obj_set_param_int(struct mrc_obj *obj, const char *name, int val);
void mrc_obj_set_param_float(struct mrc_obj *obj, const char *name, float val);
void mrc_obj_set_param_string(struct mrc_obj *obj, const char *name, const char *val);
void mrc_obj_set_param_select(struct mrc_obj *obj, const char *name, int val);
void mrc_obj_set_param_int3(struct mrc_obj *obj, const char *name, int val[3]);
void mrc_obj_set_param_float3(struct mrc_obj *obj, const char *name, float val[3]);
void mrc_obj_get_param_int(struct mrc_obj *obj, const char *name, int *pval);
void mrc_obj_get_param_string(struct mrc_obj *obj, const char *name, const char **val);
void mrc_obj_get_param_int3(struct mrc_obj *obj, const char *name, int *pval);
void mrc_obj_view(struct mrc_obj *obj);
void mrc_obj_setup(struct mrc_obj *obj);
void mrc_obj_setup_sub(struct mrc_obj *obj);
void mrc_obj_write(struct mrc_obj *obj, struct mrc_io *io);
struct mrc_obj *mrc_obj_read(struct mrc_io *io, const char *name, struct mrc_class *class);

#define MRC_OBJ_DEFINE_STANDARD_METHODS(pfx, class_type)		\
  static inline class_type *						\
  pfx ## _create(MPI_Comm comm)						\
  {									\
    return (class_type *)						\
      mrc_obj_create(comm, &mrc_class_ ## pfx);				\
  }									\
									\
  static inline class_type *						\
  pfx ## _get(class_type *obj)						\
  {									\
    return (class_type *)						\
      mrc_obj_get((struct mrc_obj *)obj);				\
  }									\
									\
  static inline void 							\
  pfx ## _destroy(class_type *obj)					\
  {									\
    mrc_obj_destroy((struct mrc_obj *)obj);				\
  }									\
									\
  static inline void 							\
  pfx ## _put(class_type *obj)						\
  {									\
    mrc_obj_put((struct mrc_obj *)obj);					\
  }									\
									\
  static inline MPI_Comm						\
  pfx ## _comm(class_type *obj)						\
  {									\
    return mrc_obj_comm((struct mrc_obj *)obj);				\
  }									\
									\
  static inline const char *						\
  pfx ## _name(class_type *obj)						\
  {									\
    return mrc_obj_name((struct mrc_obj *)obj);				\
  }									\
									\
  static inline void 							\
  pfx ## _set_name(class_type *obj, const char *name)			\
  {									\
    mrc_obj_set_name((struct mrc_obj *)obj, name);			\
  }									\
									\
  static inline void 							\
  pfx ## _set_type(class_type *obj, const char *name)			\
  {									\
    mrc_obj_set_type((struct mrc_obj *)obj, name);			\
  }									\
									\
  static inline void 							\
  pfx ## _set_param_int(class_type *obj, const char *name, int val)	\
  {									\
    mrc_obj_set_param_int((struct mrc_obj *)obj, name, val);		\
  }									\
									\
  static inline void 							\
  pfx ## _set_param_float(class_type *obj, const char *name, float val)	\
  {									\
    mrc_obj_set_param_float((struct mrc_obj *)obj, name, val);		\
  }									\
									\
  static inline void 							\
  pfx ## _set_param_string(class_type *obj, const char *name, const char *val) \
  {									\
    mrc_obj_set_param_string((struct mrc_obj *)obj, name, val);		\
  }									\
									\
  static inline void 							\
  pfx ## _set_param_select(class_type *obj, const char *name, int val)	\
  {									\
    mrc_obj_set_param_select((struct mrc_obj *)obj, name, val);		\
  }									\
									\
  static inline void 							\
  pfx ## _set_param_int3(class_type *obj, const char *name, int val[3])	\
  {									\
    mrc_obj_set_param_int3((struct mrc_obj *)obj, name, val);		\
  }									\
									\
  static inline void 							\
  pfx ## _set_param_float3(class_type *obj, const char *name, float val[3])	\
  {									\
    mrc_obj_set_param_float3((struct mrc_obj *)obj, name, val);		\
  }									\
									\
  static inline void 							\
  pfx ## _get_param_int(class_type *obj, const char *name, int *pval)	\
  {									\
    mrc_obj_get_param_int((struct mrc_obj *)obj, name, pval);		\
  }									\
									\
  static inline void 							\
  pfx ## _get_param_int3(class_type *obj, const char *name, int *pval)	\
  {									\
    mrc_obj_get_param_int3((struct mrc_obj *)obj, name, pval);		\
  }									\
									\
  static inline void 							\
  pfx ## _set_from_options(class_type *obj)				\
  {									\
    mrc_obj_set_from_options((struct mrc_obj *)obj);			\
  }									\
									\
  static inline void 							\
  pfx ## _view(class_type *obj)						\
  {									\
    mrc_obj_view((struct mrc_obj *)obj);				\
  }									\
									\
  static inline void 							\
  pfx ## _setup(class_type *obj)					\
  {									\
    mrc_obj_setup((struct mrc_obj *)obj);				\
  }									\
									\
  static inline class_type *						\
  pfx ## _read(struct mrc_io *io, const char *name)			\
  {									\
    return (class_type *) mrc_obj_read(io, name, &mrc_class_ ## pfx);	\
  }									\
									\
  static inline void 							\
  pfx ## _write(class_type *obj, struct mrc_io *io)			\
  {									\
    mrc_obj_write((struct mrc_obj *)obj, io);				\
  }									\
									\

#endif
