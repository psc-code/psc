
#include <mrc_obj.h>
#include <mrc_params.h>
#include <mrc_io.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

static struct mrc_obj *
obj_create(MPI_Comm comm, struct mrc_class *cls)
{
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
  mrc_obj_set_name(obj, cls->name);

  if (cls->param_descr) {
    char *p = (char *) obj + cls->param_offset;
    mrc_params_set_default(p, cls->param_descr);
  }

  if (cls->create) {
    cls->create(obj);
  }

  // set subclass to first in list as default
  if (!obj->ops && !list_empty(&cls->subclasses)) {
    struct mrc_obj_ops *ops =
      list_entry(cls->subclasses.next, struct mrc_obj_ops, list);
    mrc_obj_set_type(obj, ops->name);
  }

  // no ref count for this reference, will be deleted on final destroy()
  list_add_tail(&obj->instance_entry, &cls->instances);

  return obj;
}

static struct mrc_obj_ops *
find_subclass_ops(struct mrc_class *cls, const char *subclass)
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
mrc_obj_create(MPI_Comm comm, struct mrc_class *cls)
{
  init_class(cls);
  return obj_create(comm, cls);
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
  list_del(&obj->instance_entry);

  if (obj->ops) {
    if (obj->ops->destroy) {
      obj->ops->destroy(obj);
    }
  }

  free(obj->subctx);

  if (cls->destroy) {
    cls->destroy(obj);
  }

  while (!list_empty(&obj->children_list)) {
    struct mrc_obj *child = list_entry(obj->children_list.next, struct mrc_obj,
				       child_entry);
    list_del(&child->child_entry);
    mrc_obj_destroy(child);
  }

  if (obj->comm != MPI_COMM_NULL) {
    MPI_Comm_free(&obj->comm);
  }

  if (obj->name != obj->cls->name) {
    free(obj->name);
  }
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

void
mrc_obj_set_type(struct mrc_obj *obj, const char *subclass)
{
  //  printf("set_type: %s -> %s\n", obj->ops->name, subclass);
  if (obj->ops && strcmp(obj->ops->name, subclass) == 0)
    return;

  if (obj->ops && obj->ops->destroy) {
    obj->ops->destroy(obj);
  }

  free(obj->subctx);
  obj->subctx = NULL;
  
  struct mrc_obj_ops *ops = find_subclass_ops(obj->cls, subclass);
  assert(ops);
  obj->ops = ops;

  if (ops->size) {
    obj->subctx = calloc(1, ops->size);
    if (ops->param_descr) {
      char *p = (char *) obj->subctx + ops->param_offset;
      mrc_params_set_default(p, ops->param_descr);
    }
  }

  if (ops->create) {
    ops->create(obj);
  }
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

  fprintf(stderr, "ERROR: option '%s' not found (type %d)!\n", name, type);
  abort();
}

void
mrc_obj_get_param_type(struct mrc_obj *obj, const char *name,
		       int type, union param_u *uval)
{
  struct mrc_class *cls = obj->cls;
  if (cls->param_descr) {
    char *p = (char *) obj + cls->param_offset;
    if (mrc_params_get_type(p, cls->param_descr, name, type, uval) == 0)
      return;
  }
  struct mrc_obj_ops *ops = obj->ops;
  if (ops && ops->param_descr) {
    char *p = (char *) obj->subctx + ops->param_offset;
    if (mrc_params_get_type(p, ops->param_descr, name, type, uval) == 0)
      return;
  }

  fprintf(stderr, "ERROR: option '%s' not found (type %d)!\n", name, type);
  abort();
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
  union param_u uval = { .u_string = val };
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

void mrc_obj_set_param_ptr(struct mrc_obj *obj, const char *name, void* val)
{
  union param_u uval = { .u_ptr = val };
  mrc_obj_set_param_type(obj, name, PT_PTR, &uval);
}

void
mrc_obj_get_param_int(struct mrc_obj *obj, const char *name, int *pval)
{
  union param_u uval;
  mrc_obj_get_param_type(obj, name, PT_INT, &uval);
  *pval = uval.u_int;
}

void
mrc_obj_get_param_string(struct mrc_obj *obj, const char *name, const char **val)
{
  union param_u uval;
  mrc_obj_get_param_type(obj, name, PT_STRING, &uval);
  *val = uval.u_string;
}

void
mrc_obj_get_param_int3(struct mrc_obj *obj, const char *name, int *pval)
{
  union param_u uval;
  mrc_obj_get_param_type(obj, name, PT_INT3, &uval);
  for (int d = 0; d < 3; d++) {
    pval[d] = uval.u_int3[d];
  }
}

static void
mrc_obj_view_this(struct mrc_obj *obj)
{
  struct mrc_class *cls = obj->cls;
  if (cls->param_descr) {
    char *p = (char *) obj + cls->param_offset;
    mrc_params_print(p, cls->param_descr, mrc_obj_name(obj), obj->comm);
  } else {
    mpi_printf(obj->comm, 
	       "\n---------------------------------------------------- class -- %s",
	       mrc_obj_name(obj));
  } 
  if (cls->view) {
    cls->view(obj);
  } else {
    mpi_printf(obj->comm, "\n");
  }

  if (obj->ops) {
    if (obj->ops->param_descr) {
      char *p = (char *) obj->subctx + obj->ops->param_offset;
      mrc_params_print(p, obj->ops->param_descr, obj->ops->name, obj->comm);
    } else {
      mpi_printf(obj->comm, 
		 "----------------------------------------------------- type -- %s\n\n",
		 obj->ops->name);
    } 

    if (obj->ops->view) {
      obj->ops->view(obj);
    }
  }
}

void
mrc_obj_view(struct mrc_obj *obj)
{
  mrc_obj_view_this(obj);

  struct mrc_obj *child;
  __list_for_each_entry(child, &obj->children_list, child_entry, struct mrc_obj) {
    mrc_obj_view(child);
  }
}

void
mrc_obj_setup_sub(struct mrc_obj *obj)
{
  if (obj->ops && obj->ops->setup) {
    obj->ops->setup(obj);
  }
}

static void
mrc_obj_setup_this(struct mrc_obj *obj)
{
  struct mrc_class *cls = obj->cls;
  if (cls->setup) {
    cls->setup(obj);
  } else {
    mrc_obj_setup_sub(obj);
  }
  if (obj->view_flag) {
    obj->view_flag = false;
    mrc_obj_view(obj);
  }
}

void
mrc_obj_setup(struct mrc_obj *obj)
{
  mrc_obj_setup_this(obj);

  struct mrc_obj *child;
  __list_for_each_entry(child, &obj->children_list, child_entry, struct mrc_obj) {
    mrc_obj_setup(child);
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

static void
mrc_obj_read2(struct mrc_obj *obj, struct mrc_io *io)
{
  struct mrc_class *cls = obj->cls;

  if (cls->param_descr) {
    char *p = (char *) obj + cls->param_offset;
    mrc_params_read(p, cls->param_descr, mrc_obj_name(obj), io);
  }
  if (!list_empty(&cls->subclasses)) {
    char *type;
    mrc_io_read_attr_string(io, mrc_obj_name(obj), "mrc_obj_type", &type);
    mrc_obj_set_type(obj, type);
    free(type);
  }
  if (obj->ops && obj->ops->param_descr) {
    char *p = (char *) obj->subctx + obj->ops->param_offset;
    mrc_params_read(p, obj->ops->param_descr, mrc_obj_name(obj), io);
  }
  if (obj->ops && obj->ops->read) {
    obj->ops->read(obj, io);
  } else if (cls->read) {
    cls->read(obj, io);
  } else {
    mrc_obj_read_children(obj, io);
    mrc_obj_setup(obj);
  }
}

void
mrc_obj_read_children(struct mrc_obj *obj, struct mrc_io *io)
{
  struct mrc_obj *child;
  __list_for_each_entry(child, &obj->children_list, child_entry, struct mrc_obj) {
    mrc_obj_read2(child, io);
  }
}

struct mrc_obj *
mrc_obj_read(struct mrc_io *io, const char *name, struct mrc_class *cls)
{
  init_class(cls);

  struct mrc_obj *obj = mrc_io_find_obj(io, name);
  if (obj)
    return obj;

  obj = mrc_obj_create(mrc_io_comm(io), cls);
  mrc_obj_set_name(obj, name);
  mrc_io_add_obj(io, obj);

  char *s;
  mrc_io_read_attr_string(io, mrc_obj_name(obj), "mrc_obj_class", &s);
  assert(strcmp(cls->name, s) == 0);
  free(s);

  mrc_obj_read2(obj, io);

  return obj;
}

void
mrc_obj_write(struct mrc_obj *obj, struct mrc_io *io)
{
  if (mrc_io_add_obj(io, obj) == 1) // exists?
    return;

  struct mrc_class *cls = obj->cls;

  mrc_io_write_attr_string(io, mrc_obj_name(obj), "mrc_obj_class", cls->name);
  if (cls->param_descr) {
    char *p = (char *) obj + cls->param_offset;
    mrc_params_write(p, cls->param_descr, mrc_obj_name(obj), io);
  }
  if (obj->ops) {
    mrc_io_write_attr_string(io, mrc_obj_name(obj), "mrc_obj_type", obj->ops->name);
    if (obj->ops->param_descr) {
      char *p = (char *) obj->subctx + obj->ops->param_offset;
      mrc_params_write(p, obj->ops->param_descr, mrc_obj_name(obj), io);
    }
  }
  if (cls->write) {
    cls->write(obj, io);
  }
  if (obj->ops && obj->ops->write) {
    obj->ops->write(obj, io);
  }

  struct mrc_obj *child;
  __list_for_each_entry(child, &obj->children_list, child_entry, struct mrc_obj) {
    mrc_obj_write(child, io);
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

