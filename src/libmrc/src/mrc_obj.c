
#include <mrc_obj.h>
#include <mrc_params.h>
#include <mrc_io.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

static struct mrc_obj *
obj_create(MPI_Comm comm, struct mrc_class *class)
{
  assert_collective(comm);

  struct mrc_obj *obj = calloc(1, class->size);
  MPI_Comm_dup(comm, &obj->comm);

  obj->class = class;
  obj->refcount = 1;
  mrc_obj_set_name(obj, class->name);

  if (class->param_descr) {
    char *p = (char *) obj + class->param_offset;
    mrc_params_set_default(p, class->param_descr);
  }

  if (class->create) {
    class->create(obj);
  }

  // set subclass to first in list as default
  if (!obj->ops && class->subclasses && !list_empty(class->subclasses)) {
    struct mrc_obj_ops *ops =
      list_entry(class->subclasses->next, struct mrc_obj_ops, list);
    mrc_obj_set_type(obj, ops->name);
  }

  // no ref count for this reference, will be deleted on final destroy()
  list_add_tail(&obj->instance_entry, &class->instances);

  return obj;
}

static struct mrc_obj_ops *
find_subclass_ops(struct mrc_class *class, const char *subclass)
{
  if (!subclass)
    return NULL;

  struct mrc_obj_ops *ops;
  list_for_each_entry(ops, class->subclasses, list) {
    assert(ops->name);
    if (strcmp(subclass, ops->name) == 0) {
      return ops;
    }
  }

  mpi_printf(MPI_COMM_WORLD, "ERROR: unknown subclass '%s' of class '%s'\n", subclass,
	  class->name);
  mpi_printf(MPI_COMM_WORLD, "valid choices are:\n");
  list_for_each_entry(ops, class->subclasses, list) {
    mpi_printf(MPI_COMM_WORLD, "- %s\n", ops->name);
  }
  abort();
}

struct mrc_obj *
mrc_obj_create(MPI_Comm comm, struct mrc_class *class)
{
  if (!class->initialized) {
    class->initialized = true;
    INIT_LIST_HEAD(&class->instances);
    if (class->init) {
      class->init();
    }
  }

  return obj_create(comm, class);
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

  struct mrc_class *class = obj->class;

  list_del(&obj->instance_entry);

  if (obj->ops) {
    if (obj->ops->destroy) {
      obj->ops->destroy(obj);
    }
  }

  free(obj->subctx);

  if (class->destroy) {
    class->destroy(obj);
  }

  MPI_Comm_free(&obj->comm);
  if (obj->name != obj->class->name) {
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

void
mrc_obj_set_name(struct mrc_obj *obj, const char *name)
{
  if (obj->name) {
    if (strcmp(name, obj->name) == 0)
      return;

    if (obj->name != obj->class->name) {
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
    list_for_each_entry(p, &obj->class->instances, instance_entry) {
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
  if (strcmp(new_name, obj->class->name) == 0) {
    obj->name = (char *) obj->class->name;
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
  
  struct mrc_obj_ops *ops = find_subclass_ops(obj->class, subclass);
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

void
mrc_obj_set_from_options(struct mrc_obj *obj)
{
  struct mrc_class *class = obj->class;

  const char *type;
  char option[strlen(mrc_obj_name(obj)) + 6];
  sprintf(option, "%s_type", mrc_obj_name(obj));
  if (mrc_params_get_option_string(option, &type) == 0) {
    mrc_obj_set_type(obj, type);
  }
  if (class->set_from_options) {
    class->set_from_options(obj);
  }

  if (class->param_descr) {
    char *p = (char *) obj + class->param_offset;
    mrc_params_parse_nodefault(p, class->param_descr, mrc_obj_name(obj), obj->comm);
    mrc_params_parse_pfx(p, class->param_descr, mrc_obj_name(obj), obj->comm);
  }

  if (obj->ops) {
    if (obj->ops->param_descr) {
      char *p = (char *) obj->subctx + obj->ops->param_offset;
      mrc_params_parse_nodefault(p, obj->ops->param_descr, obj->ops->name, obj->comm);
      mrc_params_parse_pfx(p, obj->ops->param_descr, obj->ops->name, obj->comm);
    }
    
    if (obj->ops->set_from_options) {
      obj->ops->set_from_options(obj);
    }
  }
}

void
mrc_obj_set_param_type(struct mrc_obj *obj, const char *name,
		       int type, union param_u *uval)
{
  struct mrc_class *class = obj->class;
  if (class->param_descr) {
    char *p = (char *) obj + class->param_offset;
    if (mrc_params_set_type(p, class->param_descr, name, type, uval) == 0)
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
  struct mrc_class *class = obj->class;
  if (class->param_descr) {
    char *p = (char *) obj + class->param_offset;
    if (mrc_params_get_type(p, class->param_descr, name, type, uval) == 0)
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
mrc_obj_set_param_string(struct mrc_obj *obj, const char *name, const char *val)
{
  union param_u uval = { .u_string = val };
  mrc_obj_set_param_type(obj, name, PT_STRING, &uval);
}

void
mrc_obj_set_param_select(struct mrc_obj *obj, const char *name, int val)
{
  union param_u uval = { .u_select = val };
  mrc_obj_set_param_type(obj, name, PT_SELECT, &uval);
}

void
mrc_obj_set_param_int3(struct mrc_obj *obj, const char *name, int val[3])
{
  union param_u uval = { .u_int3 = { val[0], val[1], val[2] } };
  mrc_obj_set_param_type(obj, name, PT_INT3, &uval);
}

void
mrc_obj_set_param_float3(struct mrc_obj *obj, const char *name, float val[3])
{
  union param_u uval = { .u_float3 = { val[0], val[1], val[2] } };
  mrc_obj_set_param_type(obj, name, PT_FLOAT3, &uval);
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

void
mrc_obj_view(struct mrc_obj *obj)
{
  struct mrc_class *class = obj->class;
  if (class->param_descr) {
    char *p = (char *) obj + class->param_offset;
    mrc_params_print(p, class->param_descr, mrc_obj_name(obj), obj->comm);
  } else {
    mpi_printf(obj->comm, 
	       "\n------------------------------------------------------------- %s",
	       mrc_obj_name(obj));
  } 
  if (class->view) {
    class->view(obj);
  }

  if (obj->ops) {
    if (obj->ops->param_descr) {
      char *p = (char *) obj->subctx + obj->ops->param_offset;
      mrc_params_print(p, obj->ops->param_descr, obj->ops->name, obj->comm);
    } else {
      mpi_printf(obj->comm, 
		 "------------------------------------------------------------- %s\n\n",
		 obj->ops->name);
    } 

    if (obj->ops->view) {
      obj->ops->view(obj);
    }
  }
}

void
mrc_obj_setup_sub(struct mrc_obj *obj)
{
  if (obj->ops->setup) {
    obj->ops->setup(obj);
  }
}

void
mrc_obj_setup(struct mrc_obj *obj)
{
  struct mrc_class *class = obj->class;
  if (class->setup) {
    class->setup(obj);
  } else {
    mrc_obj_setup_sub(obj);
  }
}

struct mrc_obj *
mrc_obj_read(struct mrc_io *io, const char *name, struct mrc_class *class)
{
  // FIXME dupl
  if (!class->initialized) {
    class->initialized = true;
    if (class->init) {
      class->init();
    }
  }

  struct mrc_obj *obj = mrc_io_find_obj(io, name);
  if (obj)
    return obj;

  obj = mrc_obj_create(mrc_io_comm(io), class);
  mrc_obj_set_name(obj, name);
  mrc_io_add_obj(io, obj);

  char *s;
  mrc_io_read_attr_string(io, mrc_obj_name(obj), "mrc_obj_class", &s);
  assert(strcmp(class->name, s) == 0);
  free(s);
  if (class->param_descr) {
    char *p = (char *) obj + class->param_offset;
    mrc_params_read(p, class->param_descr, mrc_obj_name(obj), io);
  }
  if (class->subclasses) {
    char *type;
    mrc_io_read_attr_string(io, mrc_obj_name(obj), "mrc_obj_type", &type);
    mrc_obj_set_type(obj, type);
    free(type);
  }
  if (obj->ops && obj->ops->param_descr) {
    char *p = (char *) obj->subctx + obj->ops->param_offset;
    mrc_params_read(p, obj->ops->param_descr, mrc_obj_name(obj), io);
  }
  if (class->read) {
    class->read(obj, io);
  }

  mrc_obj_setup(obj);

  return obj;
}

void
mrc_obj_write(struct mrc_obj *obj, struct mrc_io *io)
{
  if (mrc_io_add_obj(io, obj) == 1) // exists?
    return;

  struct mrc_class *class = obj->class;

  mrc_io_write_attr_string(io, mrc_obj_name(obj), "mrc_obj_class", class->name);
  if (class->param_descr) {
    char *p = (char *) obj + class->param_offset;
    mrc_params_write(p, class->param_descr, mrc_obj_name(obj), io);
  }
  if (obj->ops) {
    mrc_io_write_attr_string(io, mrc_obj_name(obj), "mrc_obj_type", obj->ops->name);
    if (obj->ops->param_descr) {
      char *p = (char *) obj->subctx + obj->ops->param_offset;
      mrc_params_write(p, obj->ops->param_descr, mrc_obj_name(obj), io);
    }
  }
  if (class->write) {
    class->write(obj, io);
  }
}

