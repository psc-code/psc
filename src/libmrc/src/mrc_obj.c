
#include <mrc_obj.h>
#include <mrc_params.h>
#include <mrc_io.h>
#include <mrc_list.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <stdbool.h>
#include <execinfo.h>

static list_t active_classes_head;

int mrc_view_level = 0;

void
mrc_view_printf(MPI_Comm comm, const char *fmt, ...)
{
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank != 0)
    return;

  printf("%*.*s", 3 * (mrc_view_level - 1), 3 * (mrc_view_level - 1), "");
  va_list argp;
  va_start(argp, fmt);
  vprintf(fmt, argp);
  va_end(argp);
}

static void
create_member_objs(MPI_Comm comm, void *p, struct param *descr)
{
  for (int i = 0; descr && descr[i].name; i++) {
    if (descr[i].type == MRC_VAR_OBJ) {
      union param_u *pv = p + (unsigned long) descr[i].var;
      pv->u_obj = __mrc_obj_create(comm, descr[i].u.mrc_obj.cls);
    }
  }
}

static void
destroy_member_objs(void *p, struct param *descr)
{
  for (int i = 0; descr && descr[i].name; i++) {
    if (descr[i].type == MRC_VAR_OBJ) {
      union param_u *pv = p + (unsigned long) descr[i].var;
      mrc_obj_destroy(pv->u_obj);
    }
    if (descr[i].type == PT_INT_ARRAY) { 
      union param_u *pv = p + (unsigned long) descr[i].var;
      free(pv->u_int_array.vals);
      pv->u_int_array.nr_vals = 0;
    }
    if (descr[i].type == PT_FLOAT_ARRAY) { 
      union param_u *pv = p + (unsigned long) descr[i].var;
      free(pv->u_float_array.vals);
      pv->u_float_array.nr_vals = 0;
    }
  }
}

static void
set_from_options_member_objs(void *p, struct param *descr)
{
  for (int i = 0; descr && descr[i].name; i++) {
    if (descr[i].type == MRC_VAR_OBJ) {
      union param_u *pv = p + (unsigned long) descr[i].var;
      mrc_obj_set_from_options(pv->u_obj);
    }
  }
}

static void
view_member_objs(void *p, struct param *descr)
{
  for (int i = 0; descr && descr[i].name; i++) {
    if (descr[i].type == MRC_VAR_OBJ) {
      union param_u *pv = p + (unsigned long) descr[i].var;
      mrc_obj_view(pv->u_obj);
    }
  }
}

static void
setup_member_objs(void *p, struct param *descr)
{
  for (int i = 0; descr && descr[i].name; i++) {
    if (descr[i].type == MRC_VAR_OBJ) {
      union param_u *pv = p + (unsigned long) descr[i].var;
      mrc_obj_setup(pv->u_obj);
    }
  }
}

static void
read_member_objs(struct mrc_io *io, struct mrc_obj *obj, void *p, struct param *descr)
{
  for (int i = 0; descr && descr[i].name; i++) {
    if (descr[i].type == MRC_VAR_OBJ) {
      union param_u *pv = p + (unsigned long) descr[i].var;
      pv->u_obj = __mrc_io_read_ref(io, obj, descr[i].name, descr[i].u.mrc_obj.cls);
    }
  }
}

static struct mrc_obj *
obj_create(MPI_Comm comm, struct mrc_class *cls, bool basic_only)
{
  // this is for printing a backtrace if creating a watched class
  void *backtrace_addr[4];
  size_t size, frame;
  char **backtrace_str;

  assert_collective(comm);

  assert(cls->size >= sizeof(struct mrc_obj));
  struct mrc_obj *obj = calloc(1, cls->size);
  if (comm == MPI_COMM_NULL) {
    obj->comm = MPI_COMM_NULL;
  } else {
    MPI_Comm_dup(comm, &obj->comm);
  }

  obj->cls = cls;
  obj->refcount = 1;
  INIT_LIST_HEAD(&obj->children_list);
  INIT_LIST_HEAD(&obj->dict_list);
  mrc_obj_set_name(obj, cls->name);

  if (cls->param_descr) {
    char *p = (char *) obj + cls->param_offset;
    mrc_params_set_default(p, cls->param_descr);
  }

  // no ref count for this reference, will be deleted on final destroy()
#pragma omp critical
  list_add_tail(&obj->instance_entry, &cls->instances);

  // keep track what classes have been instantiated
  if (active_classes_head.next == NULL) {
    INIT_LIST_HEAD(&active_classes_head);
  }
  if (cls->active_classes.next == NULL) {
    list_add(&cls->active_classes, &active_classes_head);
  }
  cls->nr_instances++;
  cls->nr_creates++;

  // print extra info if we're creating a watched class
  if (cls->watch) {
    size = backtrace(backtrace_addr, 4);
    backtrace_str = backtrace_symbols(backtrace_addr, size);
    frame = size >= 4 ? 3 : size - 1;  // pick frame 4, or less...
    mprintf("Class Watcher:: (%s) created near (%s)\n",
            cls->name, backtrace_str[frame]);
    obj->creation_trace = backtrace_addr[frame];
  }

  // if we're read()'ing this, we don't call the actual create() functions
  if (basic_only) {
    return obj;
  }

  create_member_objs(comm, obj, cls->param_descr);

  if (cls->create) {
    cls->create(obj);
  }

  // set subclass to first in list as default
  if (!obj->ops && !list_empty(&cls->subclasses)) {
    struct mrc_obj_ops *ops =
      list_entry(cls->subclasses.next, struct mrc_obj_ops, list);
    mrc_obj_set_type(obj, ops->name);
  }

  return obj;
}

struct mrc_obj_ops *
mrc_obj_find_subclass_ops(struct mrc_class *cls, const char *subclass)
{
  if (!subclass)
    return NULL;

  if (list_empty(&cls->subclasses)) {
    mpi_printf(MPI_COMM_WORLD,
               "ERROR: requested subclass '%s', but class '%s' has no subclasses!\n",
               subclass, cls->name);
  }

  struct mrc_obj_ops *ops;
  __list_for_each_entry(ops, &cls->subclasses, list, struct mrc_obj_ops) {
    assert(ops->name);
    if (strcmp(subclass, ops->name) == 0) {
      return ops;
    }
  }

  mpi_printf(MPI_COMM_WORLD, "ERROR: unknown subclass '%s' of class '%s'\n", subclass,
          cls->name);
  mpi_printf(MPI_COMM_WORLD, "valid choices are:\n");
  __list_for_each_entry(ops, &cls->subclasses, list, struct mrc_obj_ops) {
    mpi_printf(MPI_COMM_WORLD, "- %s\n", ops->name);
  }
  abort();
}

static void
init_class(struct mrc_class *cls)
{
  if (!cls->initialized) {
    cls->initialized = true;
    if (!cls->subclasses.next) {
      INIT_LIST_HEAD(&cls->subclasses);
    }
    INIT_LIST_HEAD(&cls->instances);
    if (cls->init) {
      cls->init();
    }
  }
}

struct mrc_obj *
__mrc_obj_create(MPI_Comm comm, struct mrc_class *cls)
{
  init_class(cls);
  return obj_create(comm, cls, false);
}

struct mrc_class mrc_class_mrc_obj = {
  .name          = "mrc_obj",
  .size          = sizeof(struct mrc_obj),
};

struct mrc_obj *
mrc_obj_create(MPI_Comm comm)
{
  return __mrc_obj_create(comm, &mrc_class_mrc_obj);
}

struct mrc_obj *
mrc_obj_get(struct mrc_obj *obj)
{
  assert(obj->refcount > 0);
  obj->refcount++;
  return obj;
}

void
mrc_obj_put(struct mrc_obj *obj)
{
  if (--obj->refcount > 0)
    return;

  struct mrc_class *cls = obj->cls;
#pragma omp critical
  list_del(&obj->instance_entry);

  if (obj->ops) {
    if (obj->ops->destroy) {
      obj->ops->destroy(obj);
    }

    char *p = (char *) obj->subctx + obj->ops->param_offset;
    destroy_member_objs(p, obj->ops->param_descr);
  }

  free(obj->subctx);

  if (cls->destroy) {
    cls->destroy(obj);
  }

  char *p = (char *) obj + obj->cls->param_offset;
  destroy_member_objs(p, obj->cls->param_descr);

  while (!list_empty(&obj->children_list)) {
    struct mrc_obj *child = list_entry(obj->children_list.next, struct mrc_obj,
                                       child_entry);
    list_del(&child->child_entry);
    mrc_obj_destroy(child);
  }

  while (!list_empty(&obj->dict_list)) {
    struct mrc_dict_entry *p =
      list_entry(obj->dict_list.next, struct mrc_dict_entry, entry);
    if (p->prm.type == PT_STRING) {
      free((char *)p->val.u_string);
    } else if (p->prm.type == PT_INT_ARRAY) {
      if (p->val.u_int_array.vals) {
        free((int *)p->val.u_int_array.vals);
      }
    } else if (p->prm.type == PT_FLOAT_ARRAY) {
      if (p->val.u_float_array.vals) {
        free((float *)p->val.u_float_array.vals);
      }
    }
    free((char *)p->prm.name);
    list_del(&p->entry);
    free(p);
  }

  if (obj->comm != MPI_COMM_NULL) {
    MPI_Comm_free(&obj->comm);
  }

  if (obj->name != obj->cls->name) {
    free(obj->name);
  }

  // in the event we want to note when watched classes are freed
  // if (obj->cls->watch) {
  // }

  obj->cls->nr_instances -= 1;
  obj->cls->nr_frees += 1;
  // uncommenting this suppresses class info for classes with no more instances,
  // but keeping the list adds no extra overhead
  // if (obj->cls->nr_instances == 0) {
  //   list_del(&obj->cls->active_classes);
  //   obj->cls->active_classes.next = NULL;
  //   obj->cls->active_classes.prev = NULL;
  // }

  free(obj);
}

void
mrc_obj_destroy(struct mrc_obj *obj)
{
  if (!obj)
    return;

  mrc_obj_put(obj);
}

MPI_Comm
mrc_obj_comm(struct mrc_obj *obj)
{
  return obj->comm;
}

const char *
mrc_obj_name(struct mrc_obj *obj)
{
  return obj->name;
}

const char *
mrc_obj_type(struct mrc_obj *obj)
{
  if (!obj->ops) {
    return NULL;
  }
  return obj->ops->name;
}

void
mrc_obj_set_name(struct mrc_obj *obj, const char *name)
{
  if (obj->name) {
    if (strcmp(name, obj->name) == 0)
      return;

    if (obj->name != obj->cls->name) {
      free(obj->name);
      obj->name = NULL;
    }
  }
  char *new_name = strdup(name);
#if 0
  char *new_name;
  for (int i = 0; ; i++) {
    if (i == 0) {
      new_name = strdup(name);
    } else {
      new_name = malloc(strlen(name) + 10);
      sprintf(new_name, "%s_%d", name, i);
    }

    bool unique = true;
    struct mrc_obj *p;
    list_for_each_entry(p, &obj->cls->instances, instance_entry) {
      if (p->name && strcmp(p->name, new_name) == 0) {
        unique = false;
        break;
      }
    }

    if (unique) {
      break;
    }
    free(new_name);
  }
  if (strcmp(new_name, name) != 0) {
    mprintf("WARNING: renaming '%s' -> '%s'!\n", name, new_name);
  }
#endif
  if (strcmp(new_name, obj->cls->name) == 0) {
    obj->name = (char *) obj->cls->name;
    free(new_name);
  } else {
    obj->name = new_name;
  }
}

static void
obj_set_type(struct mrc_obj *obj, const char *subclass, bool basic_only)
{
  /* if (obj->ops) { */
  /*   printf("set_type: %s -> %s\n", obj->ops->name, subclass); */
  /* } */
  if (obj->ops && strcmp(obj->ops->name, subclass) == 0)
    return;

  if (obj->ops && obj->ops->destroy) {
    obj->ops->destroy(obj);

    char *p = (char *) obj->subctx + obj->ops->param_offset;
    destroy_member_objs(p, obj->ops->param_descr);
  }

  free(obj->subctx);
  obj->subctx = NULL;
  
  struct mrc_obj_ops *ops = mrc_obj_find_subclass_ops(obj->cls, subclass);
  assert(ops);
  obj->ops = ops;

  if (ops->size) {
    obj->subctx = calloc(1, ops->size);
    if (ops->param_descr) {
      char *p = (char *) obj->subctx + ops->param_offset;
      mrc_params_set_default(p, ops->param_descr);
    }
  }

  if (basic_only) {
    return;
  }

  char *p = (char *) obj->subctx + ops->param_offset;
  create_member_objs(obj->comm, p, ops->param_descr);

  if (ops->create) {
    ops->create(obj);
  }
}

void
mrc_obj_set_type(struct mrc_obj *obj, const char *subclass)
{
  obj_set_type(obj, subclass, false);
}

static void
mrc_obj_set_from_options_this(struct mrc_obj *obj)
{
  struct mrc_class *cls = obj->cls;

  const char help_type[] = "Choose the type (subclass) for the given class.";
  const char *type;
  char option[strlen(mrc_obj_name(obj)) + 6];
  sprintf(option, "%s_type", mrc_obj_name(obj));
  if (mrc_params_get_option_string_help(option, &type, help_type) == 0) {
    mrc_obj_set_type(obj, type);
  }

  const char help_view[] = "Print out info about the object.";
  sprintf(option, "%s_view", mrc_obj_name(obj));
  mrc_params_get_option_bool_help(option, &obj->view_flag, help_view);

  if (cls->param_descr) {
    char *p = (char *) obj + cls->param_offset;
    mrc_params_parse_nodefault(p, cls->param_descr, mrc_obj_name(obj), obj->comm);
    mrc_params_parse_pfx(p, cls->param_descr, mrc_obj_name(obj), obj->comm);
  }

  if (obj->ops) {
    if (obj->ops->param_descr) {
      char *p = (char *) obj->subctx + obj->ops->param_offset;
      mrc_params_parse_nodefault(p, obj->ops->param_descr, obj->ops->name, obj->comm);
      mrc_params_parse_pfx(p, obj->ops->param_descr, mrc_obj_name(obj), obj->comm);
    }
    
    if (obj->ops->set_from_options) {
      obj->ops->set_from_options(obj);
    }
  }

  if (cls->set_from_options) {
    cls->set_from_options(obj);
  }
}

void
mrc_obj_set_from_options(struct mrc_obj *obj)
{
  mrc_obj_set_from_options_this(obj);

  char *p = (char *) obj + obj->cls->param_offset;
  set_from_options_member_objs(p, obj->cls->param_descr);

  if (obj->ops) {
    char *p = (char *) obj->subctx + obj->ops->param_offset;
    set_from_options_member_objs(p, obj->ops->param_descr);
  }

  struct mrc_obj *child;
  __list_for_each_entry(child, &obj->children_list, child_entry, struct mrc_obj) {
    mrc_obj_set_from_options(child);
  }
}

void
mrc_obj_set_param_type(struct mrc_obj *obj, const char *name,
                       int type, union param_u *uval)
{
  struct mrc_class *cls = obj->cls;
  if (cls->param_descr) {
    char *p = (char *) obj + cls->param_offset;
    if (mrc_params_set_type(p, cls->param_descr, name, type, uval) == 0)
      return;
  }
  struct mrc_obj_ops *ops = obj->ops;
  if (ops && ops->param_descr) {
    char *p = (char *) obj->subctx + ops->param_offset;
    if (mrc_params_set_type(p, ops->param_descr, name, type, uval) == 0)
      return;
  }

  struct mrc_dict_entry *e;
  __list_for_each_entry(e, &obj->dict_list, entry, struct mrc_dict_entry) {
    if (strcmp(e->prm.name, name) == 0) {
      assert(e->prm.type == type);
      e->val = *uval;
      return;
    }
  }

  fprintf(stderr, "ERROR: option '%s' not found (type %d)!\n", name, type);
  abort();
}

int
mrc_obj_get_param_type(struct mrc_obj *obj, const char *name,
                       int type, union param_u *uval)
{
  struct mrc_class *cls = obj->cls;
  if (cls->param_descr) {
    char *p = (char *) obj + cls->param_offset;
    if (mrc_params_get_type(p, cls->param_descr, name, type, uval) == 0)
      return 0;
  }
  struct mrc_obj_ops *ops = obj->ops;
  if (ops && ops->param_descr) {
    char *p = (char *) obj->subctx + ops->param_offset;
    if (mrc_params_get_type(p, ops->param_descr, name, type, uval) == 0)
      return 0;
  }

  struct mrc_dict_entry *e;
  __list_for_each_entry(e, &obj->dict_list, entry, struct mrc_dict_entry) {
    if (strcmp(e->prm.name, name) == 0) {
      assert(e->prm.type == type);
      *uval = e->val;
      return 0;
    }
  }

  fprintf(stderr, "WARNING: option '%s' not found (type %d)!\n", name, type);
  return -1;
}

void
mrc_obj_set_param_int(struct mrc_obj *obj, const char *name, int val)
{
  union param_u uval = { .u_int = val };
  mrc_obj_set_param_type(obj, name, PT_INT, &uval);
}

void
mrc_obj_set_param_float(struct mrc_obj *obj, const char *name, float val)
{
  union param_u uval = { .u_float = val };
  mrc_obj_set_param_type(obj, name, PT_FLOAT, &uval);
}

void
mrc_obj_set_param_double(struct mrc_obj *obj, const char *name, double val)
{
  union param_u uval = { .u_double = val };
  mrc_obj_set_param_type(obj, name, PT_DOUBLE, &uval);
}

void
mrc_obj_set_param_string(struct mrc_obj *obj, const char *name, const char *val)
{
  union param_u uval = { .u_string = strdup(val) };
  mrc_obj_set_param_type(obj, name, PT_STRING, &uval);
}

void
mrc_obj_set_param_bool(struct mrc_obj *obj, const char *name, bool val)
{
  union param_u uval = { .u_bool = val };
  mrc_obj_set_param_type(obj, name, PT_BOOL, &uval);
}

void
mrc_obj_set_param_select(struct mrc_obj *obj, const char *name, int val)
{
  union param_u uval = { .u_select = val };
  mrc_obj_set_param_type(obj, name, PT_SELECT, &uval);
}

void
mrc_obj_set_param_int3(struct mrc_obj *obj, const char *name, const int val[3])
{
  union param_u uval = { .u_int3 = { val[0], val[1], val[2] } };
  mrc_obj_set_param_type(obj, name, PT_INT3, &uval);
}

void
mrc_obj_set_param_float3(struct mrc_obj *obj, const char *name, const float val[3])
{
  union param_u uval = { .u_float3 = { val[0], val[1], val[2] } };
  mrc_obj_set_param_type(obj, name, PT_FLOAT3, &uval);
}

void
mrc_obj_set_param_double3(struct mrc_obj *obj, const char *name, const double val[3])
{
  union param_u uval = { .u_double3 = { val[0], val[1], val[2] } };
  mrc_obj_set_param_type(obj, name, PT_DOUBLE3, &uval);
}

void
mrc_obj_set_param_int_array(struct mrc_obj *obj, const char *name,
                            int nr_vals, const int val[])
{
  union param_u uval = { .u_int_array = { nr_vals, (int *) val } };
  mrc_obj_set_param_type(obj, name, PT_INT_ARRAY, &uval);
}

void
mrc_obj_set_param_float_array(struct mrc_obj *obj, const char *name,
                            int nr_vals, const float val[])
{
  union param_u uval = { .u_float_array = { nr_vals, (float *) val } };
  mrc_obj_set_param_type(obj, name, PT_FLOAT_ARRAY, &uval);
}

void
mrc_obj_set_param_ptr(struct mrc_obj *obj, const char *name, void* val)
{
  union param_u uval = { .u_ptr = val };
  mrc_obj_set_param_type(obj, name, PT_PTR, &uval);
}

void
mrc_obj_set_param_obj(struct mrc_obj *obj, const char *name, void* val)
{
  union param_u uval = { .u_obj = val };
  mrc_obj_set_param_type(obj, name, PT_OBJ, &uval);
}

int
mrc_obj_get_param_bool(struct mrc_obj *obj, const char *name, bool *pval)
{
  union param_u uval;
  int err = mrc_obj_get_param_type(obj, name, PT_BOOL, &uval);
  *pval = uval.u_bool;
  return err;
}

int
mrc_obj_get_param_int(struct mrc_obj *obj, const char *name, int *pval)
{
  union param_u uval;
  int err = mrc_obj_get_param_type(obj, name, PT_INT, &uval);
  *pval = uval.u_int;
  return err;
}

int
mrc_obj_get_param_float(struct mrc_obj *obj, const char *name, float *pval)
{
  union param_u uval;
  int err = mrc_obj_get_param_type(obj, name, PT_FLOAT, &uval);
  *pval = uval.u_float;
  return err;
}

int
mrc_obj_get_param_double(struct mrc_obj *obj, const char *name, double *pval)
{
  union param_u uval;
  int err = mrc_obj_get_param_type(obj, name, PT_DOUBLE, &uval);
  *pval = uval.u_double;
  return err;
}

int
mrc_obj_get_param_string(struct mrc_obj *obj, const char *name, const char **val)
{
  union param_u uval;
  int err = mrc_obj_get_param_type(obj, name, PT_STRING, &uval);
  *val = uval.u_string;
  return err;
}

int
mrc_obj_get_param_int3(struct mrc_obj *obj, const char *name, int *pval)
{
  union param_u uval;
  int err = mrc_obj_get_param_type(obj, name, PT_INT3, &uval);
  for (int d = 0; d < 3; d++) {
    pval[d] = uval.u_int3[d];
  }
  return err;
}

int
mrc_obj_get_param_float3(struct mrc_obj *obj, const char *name, float *pval)
{
  union param_u uval;
  int err = mrc_obj_get_param_type(obj, name, PT_FLOAT3, &uval);
  for (int d = 0; d < 3; d++) {
    pval[d] = uval.u_float3[d];
  }
  return err;
}

int
mrc_obj_get_param_double3(struct mrc_obj *obj, const char *name, double *pval)
{
  union param_u uval;
  int err = mrc_obj_get_param_type(obj, name, PT_DOUBLE3, &uval);
  for (int d = 0; d < 3; d++) {
    pval[d] = uval.u_double3[d];
  }
  return err;
}

int
mrc_obj_get_param_obj(struct mrc_obj *obj, const char *name, struct mrc_obj **pval)
{
  union param_u uval;
  int rc = mrc_obj_get_param_type(obj, name, PT_OBJ, &uval);
  if (rc < 0) {
    return rc;
  }

  *pval = uval.u_obj;
  return 0;
}

int
mrc_obj_get_param_ptr(struct mrc_obj *obj, const char *name, void **val)
{
  union param_u uval;
  int err = mrc_obj_get_param_type(obj, name, PT_PTR, &uval);
  *val = uval.u_ptr;
  return err;
}

int
mrc_obj_get_param_float_array_nr_vals(struct mrc_obj *obj, const char *name, int *nr_vals)
{
  union param_u uval;
  int err = mrc_obj_get_param_type(obj, name, PT_FLOAT_ARRAY, &uval);
  *nr_vals = uval.u_float_array.nr_vals;
  return err;
}

int
mrc_obj_get_param_float_array(struct mrc_obj *obj, const char *name, float *pval)
{
  union param_u uval;
  int err = mrc_obj_get_param_type(obj, name, PT_FLOAT_ARRAY, &uval);
  for (int d = 0; d < uval.u_float_array.nr_vals; d++) {
    pval[d] = uval.u_float_array.vals[d];
  }
  return err;
}

static void
mrc_obj_view_this(struct mrc_obj *obj)
{
  struct mrc_class *cls = obj->cls;
  MPI_Comm comm = obj->comm;

  mrc_view_printf(comm, "==================================================== class == %s\n",
                  mrc_obj_name(obj));

  if (cls->param_descr || !list_empty(&obj->dict_list) || 
      (obj->ops && obj->ops->param_descr)) {
    mrc_view_printf(comm, "%-20s| %s\n", "parameter", "value");
  }

  if (cls->param_descr) {
    mrc_view_printf(comm, "--------------------+----------------------------------------\n");
    char *p = (char *) obj + cls->param_offset;
    for (int i = 0; cls->param_descr[i].name; i++) {
      mrc_params_print_one(p, &cls->param_descr[i], comm);
    }
  }

  struct mrc_dict_entry *e;
  __list_for_each_entry(e, &obj->dict_list, entry, struct mrc_dict_entry) {
    mrc_params_print_one(&e->val, &e->prm, comm);
  }

  if (obj->ops) {
    mrc_view_printf(comm, "--------------------+-------------------------------- type -- %s\n",
               obj->ops->name);
    if (obj->ops->param_descr) {
      char *p = (char *) obj->subctx + obj->ops->param_offset;
      for (int i = 0; obj->ops->param_descr[i].name; i++) {
        mrc_params_print_one(p, &obj->ops->param_descr[i], comm);
      }
    } 
  }

  if (cls->view) {
    cls->view(obj);
  }

  if (obj->ops && obj->ops->view) {
    obj->ops->view(obj);
  }

  mrc_view_printf(comm, "\n");
}

void
mrc_obj_view(struct mrc_obj *obj)
{
  mrc_view_level++;
  mrc_obj_view_this(obj);

  char *p = (char *) obj + obj->cls->param_offset;
  view_member_objs(p, obj->cls->param_descr);

  if (obj->ops) {
    char *p = (char *) obj->subctx + obj->ops->param_offset;
    view_member_objs(p, obj->ops->param_descr);
  }

  struct mrc_obj *child;
  __list_for_each_entry(child, &obj->children_list, child_entry, struct mrc_obj) {
    mrc_obj_view(child);
  }
  mrc_view_level--;
}

void
mrc_obj_setup_member_objs(struct mrc_obj *obj)
{
  char *p = (char *) obj + obj->cls->param_offset;
  setup_member_objs(p, obj->cls->param_descr);
}

void
mrc_obj_setup_member_objs_sub(struct mrc_obj *obj)
{
  assert(obj->ops);
  char *p = (char *) obj->subctx + obj->ops->param_offset;
  setup_member_objs(p, obj->ops->param_descr);
}

void
mrc_obj_read_member_objs(struct mrc_obj *obj, struct mrc_io *io)
{
  char *p = (char *) obj + obj->cls->param_offset;
  read_member_objs(io, obj, p, obj->cls->param_descr);
}

void
mrc_obj_read_member_objs_sub(struct mrc_obj *obj, struct mrc_io *io)
{
  assert(obj->ops);
  char *p = (char *) obj->subctx + obj->ops->param_offset;
  read_member_objs(io, obj, p, obj->ops->param_descr);
}

void
mrc_obj_setup_children(struct mrc_obj *obj)
{
  struct mrc_obj *child;
  __list_for_each_entry(child, &obj->children_list, child_entry, struct mrc_obj) {
    mrc_obj_setup(child);
  }
}

static void
mrc_obj_setup_default(struct mrc_obj *obj)
{
  mrc_obj_setup_member_objs(obj);
  mrc_obj_setup_children(obj);
}

// to be called internally from subclass's setup() to do its superclass setup

void
mrc_obj_setup_super(struct mrc_obj *obj)
{
  struct mrc_class *cls = obj->cls;

  if (cls->setup) {
    cls->setup(obj);
  } else {
    mrc_obj_setup_default(obj);
  }
  if (obj->view_flag) {
    obj->view_flag = false;
    mrc_obj_view(obj);
  }
}

bool
mrc_obj_is_setup(struct mrc_obj *obj)
{
  return obj->is_setup;
}

void
mrc_obj_setup(struct mrc_obj *obj)
{
  if (obj->is_setup) {
    mprintf("WARNING: %s/%p is set up twice!\n\n", obj->cls->name, obj);
    assert(0);
  }
  obj->is_setup = true;

  if (obj->ops && obj->ops->setup) {
    obj->ops->setup(obj);
  } else {
    if (obj->ops) {
      mrc_obj_setup_member_objs_sub(obj);
    }
    mrc_obj_setup_super(obj);
  }
}

void
mrc_obj_add_child(struct mrc_obj *obj, struct mrc_obj *child)
{
  // make sure we're only one parent's child
  assert(!child->child_entry.next);

  // we're not taking a reference count explicitly,
  // we are taking over the reference of the child object passed
  // in instead, so the caller relinquishes control of the object and
  // should not call ::destroy()
  list_add_tail(&child->child_entry, &obj->children_list);
}

struct mrc_obj *
mrc_obj_find_child(struct mrc_obj *obj, const char *name)
{
  struct mrc_obj *child;
  __list_for_each_entry(child, &obj->children_list, child_entry, struct mrc_obj) {
    if (strcmp(child->name, name) == 0)
      return child;
  }
  return NULL;
}

void
mrc_obj_read_super(struct mrc_obj *obj, struct mrc_io *io)
{
  struct mrc_class *cls = obj->cls;

  if (cls->read) {
    cls->read(obj, io);
  } else {
    mrc_obj_read_member_objs(obj, io);
    mrc_obj_read_children(obj, io);
    // FIXME, ugly: basically the same as mrc_obj_setup(), but skipping the children
    // setup
    obj->is_setup = true;

    if (obj->ops && !obj->ops->read && obj->ops->setup) {
      obj->ops->setup(obj);
    } else {
      if (cls->setup) {
        cls->setup(obj);
      }
    }
  }
}

static void
mrc_obj_read_params(struct mrc_obj *obj, void *p, struct param *params,
                    const char *path, struct mrc_io *io)
{
  for (int i = 0; params[i].name; i++) {
    struct param *prm = &params[i];
    union param_u *pv = p + (unsigned long) prm->var;
    if (prm->type == PT_OBJ) {
      pv->u_obj = __mrc_io_read_ref(io, obj, prm->name, prm->u.mrc_obj.cls);
    } else if (prm->type == MRC_VAR_OBJ) {
      // do nothing, member objects are read explicitly by calling 
      // mrc_obj_read_member_objs()
    } else {
      mrc_io_read_attr(io, path, prm->type, prm->name, pv);
    }
  }
}

static void
mrc_obj_read_dict(struct mrc_obj *obj, const char *path, struct mrc_io *io)
{
  int cnt;
  mrc_io_read_attr_int(io, path, "mrc_obj_dict_count", &cnt);

  for (int i = 0; i < cnt; i++) {
    char s[30], *name;
    int type;
    snprintf(s, 30, "mrc_obj_dict_name_%d", i);
    mrc_io_read_attr_string(io, path, s, &name);
    snprintf(s, 30, "mrc_obj_dict_type_%d", i);
    mrc_io_read_attr_int(io, path, s, &type);

    union param_u pv;
    if (type == PT_OBJ || type == MRC_VAR_OBJ) {
      mpi_printf(mrc_io_comm(io),
                 "!!! WARNING: cannot read back dictionary object %s (type %d) in %s!!!\n",
                 name, type, mrc_obj_name(obj));
      //pv->u_obj = __mrc_io_read_ref(io, obj, name, descr[i].u.mrc_obj.cls);
    } else {
      mrc_io_read_attr(io, path, type, name, &pv);
      mrc_obj_dict_add(obj, type, name, &pv);
    }

    free(name);
  }
}

static void
mrc_obj_read2(struct mrc_obj *obj, struct mrc_io *io, const char *path)
{
  struct mrc_class *cls = obj->cls;

  mrc_io_add_obj(io, obj, path);

  char *s;
  mrc_io_read_attr_string(io, path, "mrc_obj_name", &s);
  mrc_obj_set_name(obj, s);
  free(s);
  mrc_io_read_attr_string(io, path, "mrc_obj_class", &s);
  assert(strcmp(cls->name, s) == 0);
  free(s);

  if (cls->param_descr) {
    char *p = (char *) obj + cls->param_offset;
    mrc_obj_read_params(obj, p, cls->param_descr, path, io);
  }
  if (!list_empty(&cls->subclasses)) {
    char *type;
    mrc_io_read_attr_string(io, path, "mrc_obj_type", &type);
    obj_set_type(obj, type, true);
    free(type);
  }
  if (obj->ops && obj->ops->param_descr) {
    char *p = (char *) obj->subctx + obj->ops->param_offset;
    mrc_obj_read_params(obj, p, obj->ops->param_descr, path, io);
  }

  mrc_obj_read_dict(obj, path, io);

  if (obj->ops && obj->ops->read) {
    obj->ops->read(obj, io);
  } else {
    // This change to call order makes it more consistent with C++,
    // which implicitly calls base class constructors before 
    // derived class constructors
    mrc_obj_read_super(obj, io);
    if (obj->ops && obj->ops->create) {
      obj->ops->create(obj);
    }
  }
}

void
mrc_obj_read_children(struct mrc_obj *obj, struct mrc_io *io)
{
  struct mrc_obj *child;
  int cnt = 0;
  __list_for_each_entry(child, &obj->children_list, child_entry, struct mrc_obj) {
    char name[20]; sprintf(name, "child%d", cnt++);
    char *s;
    mrc_io_read_string(io, obj, name, &s);
    mrc_obj_read2(child, io, s);
    free(s);
  }
}

struct mrc_obj *
mrc_obj_read(struct mrc_io *io, const char *path, struct mrc_class *cls)
{
  init_class(cls);

  struct mrc_obj *obj = mrc_io_find_obj(io, path);
  if (obj) {
    assert(obj->cls == cls);
    return obj;
  }

  obj = obj_create(mrc_io_comm(io), cls, true);
  mrc_obj_read2(obj, io, path);
  return obj;
}

struct mrc_obj *
mrc_obj_read_comm(struct mrc_io *io, const char *path, struct mrc_class *cls,
		  MPI_Comm comm)
{
  init_class(cls);

  struct mrc_obj *obj = mrc_io_find_obj(io, path);
  if (obj) {
    assert(obj->cls == cls);
    return obj;
  }

  obj = obj_create(comm, cls, true);
  mrc_obj_read2(obj, io, path);
  return obj;
}

static void
mrc_obj_write_params(struct mrc_obj *obj, void *p, struct param *params,
                     const char *path, struct mrc_io *io)
{
  for (int i = 0; params[i].name; i++) {
    struct param *prm = &params[i];
    union param_u *pv = (union param_u *) (p + (unsigned long) prm->var);
    if (prm->type == PT_OBJ ||
        prm->type == MRC_VAR_OBJ) {
      mrc_io_write_ref(io, obj, prm->name, pv->u_obj);
    } else {
      mrc_io_write_attr(io, path, prm->type, prm->name, pv);
    }
  }
}

static void
mrc_obj_write_dict(struct mrc_obj *obj, const char *path, struct mrc_io *io)
{
  int cnt = 0;
  struct mrc_dict_entry *e;
  __list_for_each_entry(e, &obj->dict_list, entry, struct mrc_dict_entry) {
    char s[30];
    snprintf(s, 30, "mrc_obj_dict_name_%d", cnt);
    mrc_io_write_attr_string(io, path, s, e->prm.name);
    // FIXME, this could (possibly) be done more elegantly by getting the kind
    // from the attribute directly on read
    snprintf(s, 30, "mrc_obj_dict_type_%d", cnt);
    mrc_io_write_attr_int(io, path, s, e->prm.type);
    cnt++;
  }
  mrc_io_write_attr_int(io, path, "mrc_obj_dict_count", cnt);

  __list_for_each_entry(e, &obj->dict_list, entry, struct mrc_dict_entry) {
    if (e->prm.type == PT_OBJ ||
        e->prm.type == MRC_VAR_OBJ) {
      mrc_io_write_ref(io, obj, e->prm.name, e->val.u_obj);
    } else {
      mrc_io_write_attr(io, path, e->prm.type, e->prm.name, &e->val);
    }
  }
}

static unsigned long
mrc_obj_uid(struct mrc_obj *obj)
{
  unsigned long uid = (unsigned long) obj;
  if (obj->comm != MPI_COMM_NULL) {
    MPI_Bcast(&uid, 1, MPI_LONG, 0, obj->comm);
  }
  return uid;
}

void
mrc_obj_write(struct mrc_obj *obj, struct mrc_io *io)
{
  char path[strlen(obj->name) + 20];
  sprintf(path, "%s-uid-%#lx", obj->name, mrc_obj_uid(obj));
  if (mrc_io_add_obj(io, obj, path) == 1) // exists?
    return;

  struct mrc_class *cls = obj->cls;

  mrc_io_write_attr_string(io, path, "mrc_obj_name", obj->name);
  mrc_io_write_attr_string(io, path, "mrc_obj_class", cls->name);
  if (cls->param_descr) {
    char *p = (char *) obj + cls->param_offset;
    mrc_obj_write_params(obj, p, cls->param_descr, path, io);
  }

  mrc_obj_write_dict(obj, path, io);

  if (obj->ops) {
    mrc_io_write_attr_string(io, path, "mrc_obj_type", obj->ops->name);
    if (obj->ops->param_descr) {
      char *p = (char *) obj->subctx + obj->ops->param_offset;
      mrc_obj_write_params(obj, p, obj->ops->param_descr, path, io);
    }
  }
  if (cls->write) {
    cls->write(obj, io);
  }
  if (obj->ops && obj->ops->write) {
    obj->ops->write(obj, io);
  }

  struct mrc_obj *child;
  int cnt = 0;
  __list_for_each_entry(child, &obj->children_list, child_entry, struct mrc_obj) {
    char name[20]; sprintf(name, "child%d", cnt++);
    mrc_io_write_path(io, path, name, child);
  }
}

void
__mrc_class_register_subclass(struct mrc_class *cls, struct mrc_obj_ops *ops)
{
  if (!cls->subclasses.next) {
    INIT_LIST_HEAD(&cls->subclasses);
  }
  list_add_tail(&ops->list, &cls->subclasses);
}

void
mrc_obj_dict_add(struct mrc_obj *obj, int type, const char *name,
                 union param_u *pv)
{
  struct mrc_dict_entry *p = calloc(1, sizeof(*p));
  p->prm.type = type;
  p->prm.name = strdup(name);
  if (type == PT_STRING) {
    p->val.u_string = strdup(pv->u_string);
  } else {
    p->val = *pv;
  }
  list_add_tail(&p->entry, &obj->dict_list);
}

void
mrc_obj_dict_add_int(struct mrc_obj *obj, const char *name, int val)
{
  union param_u uval;
  uval.u_int = val;
  mrc_obj_dict_add(obj, PT_INT, name, &uval);
}

void
mrc_obj_dict_add_bool(struct mrc_obj *obj, const char *name, bool val)
{
  union param_u uval;
  uval.u_bool = val;
  mrc_obj_dict_add(obj, PT_BOOL, name, &uval);
}

void
mrc_obj_dict_add_float(struct mrc_obj *obj, const char *name, float val)
{
  union param_u uval;
  uval.u_float = val;
  mrc_obj_dict_add(obj, PT_FLOAT, name, &uval);
}

void
mrc_obj_dict_add_double(struct mrc_obj *obj, const char *name, double val)
{
  union param_u uval;
  uval.u_double = val;
  mrc_obj_dict_add(obj, PT_DOUBLE, name, &uval);
}

void
mrc_obj_dict_add_string(struct mrc_obj *obj, const char *name, const char *val)
{
  union param_u uval;
  uval.u_string = val;
  mrc_obj_dict_add(obj, PT_STRING, name, &uval);
}

void
mrc_obj_dict_add_obj(struct mrc_obj *obj, const char *name, struct mrc_obj *val)
{
  union param_u uval;
  uval.u_obj = val;
  mrc_obj_dict_add(obj, PT_OBJ, name, &uval);
}

void
mrc_obj_dict_add_float_array(struct mrc_obj *obj, const char *name, 
    float *vals, int nr_vals)
{
  union param_u uval;
  uval.u_float_array.nr_vals = nr_vals;
  if (nr_vals == 0) {
    uval.u_float_array.vals = NULL;
  } else {
    uval.u_float_array.vals = calloc(nr_vals, sizeof(float));
  }
  for (int i = 0; i < nr_vals; i++) {
    uval.u_float_array.vals[i] = vals[i];
  }
  mrc_obj_dict_add(obj, PT_FLOAT_ARRAY, name, &uval);
}

static void
get_var(void *p, struct param *descr, const char *name, int type, union param_u **pv)
{
  for (int i = 0; descr && descr[i].name; i++) {
    if (strcmp(descr[i].name, name) == 0) {
      assert(type < 0 || descr[i].type == type);
      *pv = p + (unsigned long) descr[i].var;
      return;
    }
  }
  *pv = NULL;
}

static int
mrc_obj_get_var_type(struct mrc_obj *obj, const char *name, int type,
                     union param_u **pv)
{
  // try to find variable 'name' in the class
  if (obj->cls->param_descr) {
    char *p = (char *) obj + obj->cls->param_offset;
    get_var(p, obj->cls->param_descr, name, type, pv);
    if (*pv)
      return 0;
  }

  // try to find variable 'name' in the subclass
  if (obj->ops && obj->ops->param_descr) {
    char *p = (char *) obj->subctx + obj->ops->param_offset;
    get_var(p, obj->ops->param_descr, name, type, pv);
    if (*pv)
      return 0;
  }
  
  // finally, try the dict
  struct mrc_dict_entry *e;
  __list_for_each_entry(e, &obj->dict_list, entry, struct mrc_dict_entry) {
    if (strcmp(e->prm.name, name) == 0) {
      *pv = &e->val;
      return 0;
    }
  }
  mprintf("WARNING: var '%s' not found!\n", name);
  *pv = NULL;
  return -1;
}

int
mrc_obj_get_var(struct mrc_obj *obj, const char *name, union param_u **pv)
{
  return mrc_obj_get_var_type(obj, name, -1, pv);
}

struct mrc_obj *
mrc_obj_get_var_obj(struct mrc_obj *obj, const char *name)
{
  union param_u *pv;
  mrc_obj_get_var(obj, name, &pv);
  if (!pv)
    return NULL;

  return pv->u_obj;
}

int
mrc_obj_get_var_double(struct mrc_obj *obj, const char *name, double *pval)
{
  union param_u *pv;
  int rc = mrc_obj_get_var_type(obj, name, MRC_VAR_DOUBLE, &pv);
  if (rc) {
    return rc;
  }
  *pval = pv->u_double;
  return 0;
}

mrc_void_func_t
mrc_obj_get_method(struct mrc_obj *obj, const char *name)
{
  if (obj->ops && obj->ops->methods) {
    struct mrc_obj_method *methods = obj->ops->methods;
    for (int i = 0; methods[i].name; i++) {
      if (strcmp(name, methods[i].name) == 0) {
        return methods[i].func;
      }
    }
  }

  if (obj->cls->methods) {
    struct mrc_obj_method *methods = obj->cls->methods;
    for (int i = 0; methods[i].name; i++) {
      if (strcmp(name, methods[i].name) == 0) {
        return methods[i].func;
      }
    }
  }

  return NULL;
}

// print classes that have been created, how the number of instances has
// changed since the last call, and how many times that class has been
// created / destroyed.
// Returns the total number of instances of all classes
// @verbosity: CLASS_INFO_VERB_NONE => print nothing
//             CLASS_INFO_VERB_DIFF => print only the changes since last call
//             CLASS_INFO_VERB_ACTIVE => print all classes with >= 1 instance
//             CLASS_INFO_VERB_FULL => print all classes that have been created
int
mrc_obj_print_class_info(int verbosity) {
  int total_instances = 0;
  struct mrc_class *cls;
  struct mrc_obj *instance;
  bool do_print;

  list_for_each_entry(cls, &active_classes_head, active_classes) {
    total_instances += cls->nr_instances;

    do_print = (verbosity == CLASS_INFO_VERB_DIFF && 
                ((cls->nr_instances != cls->nr_instances_last_printed) || 
                 (cls->nr_creates != cls->nr_creates_last_printed)));
    do_print |= (verbosity == CLASS_INFO_VERB_ACTIVE && cls->nr_instances > 0);
    do_print |= (verbosity >= CLASS_INFO_VERB_FULL);

    if (do_print) {
      mprintf("Class Info:: (%s):  %d (%d change)  %d+  %d-\n",
              cls->name, 
              cls->nr_instances, cls->nr_instances - cls->nr_instances_last_printed,
              cls->nr_creates, cls->nr_frees);
      cls->nr_instances_last_printed = cls->nr_instances;
      cls->nr_creates_last_printed = cls->nr_creates;    
      assert(cls->nr_instances == cls->nr_creates - cls->nr_frees);
    }

    // print watcher information
    if (verbosity >= CLASS_INFO_VERB_ACTIVE && cls->watch) {  
      list_for_each_entry(instance, &cls->instances, instance_entry) {
        mprintf("Class Watcher:: I remember a (%s, %s) was created near (%p)\n",
                cls->name, instance->name, instance->creation_trace);
      }
    }
  }

  return total_instances;
}
