
#include <mrc_mod.h>
#include <mrc_common.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

struct mrc_mod {
  struct mrc_obj obj;

  list_t list;
  MPI_Comm mod_comm;
  struct mrc_mod_entry *this_mod_entry;
};

struct mrc_mod_entry {
  char *name;
  int nr_procs;
  void (*func)(struct mrc_mod *mod, void *arg);
  void *arg;

  int rank;
  list_t entry;
};

// ======================================================================

// ----------------------------------------------------------------------
// mrc_mod_create

static void
_mrc_mod_create(struct mrc_mod *mod)
{
  INIT_LIST_HEAD(&mod->list);
  mod->mod_comm = MPI_COMM_NULL;
}

// ----------------------------------------------------------------------
// mrc_mod_destroy

static void
_mrc_mod_destroy(struct mrc_mod *mod)
{
  if (mod->mod_comm != MPI_COMM_NULL) {
    MPI_Comm_free(&mod->mod_comm);
  }
  while (!list_empty(&mod->list)) {
    struct mrc_mod_entry *p =
      list_entry(mod->list.next, struct mrc_mod_entry, entry);
    free(p->name);
    list_del(&p->entry);
    free(p);
  }
}

// ----------------------------------------------------------------------
// mrc_mod_view

static void
_mrc_mod_view(struct mrc_mod *mod)
{
  mpi_printf(mrc_mod_comm(mod), "\n");
  struct mrc_mod_entry *p;
  __list_for_each_entry(p, &mod->list, entry, struct mrc_mod_entry) {
    mpi_printf(mrc_mod_comm(mod), "%19d | '%s': %d proc%s\n", p->rank,
	       p->name, p->nr_procs, p->nr_procs == 1 ? "" : "s");
  }
}

// ----------------------------------------------------------------------
// mrc_mod_setup

static void
_mrc_mod_setup(struct mrc_mod *mod)
{
  int rank, size;
  MPI_Comm_rank(mrc_mod_comm(mod), &rank);
  MPI_Comm_size(mrc_mod_comm(mod), &size);

  int color = 0, nr_procs = 0, nr_colors = 0;
  struct mrc_mod_entry *p;
  __list_for_each_entry(p, &mod->list, entry, struct mrc_mod_entry) {
    p->rank = nr_procs;
    nr_procs += p->nr_procs;
    nr_colors++;
    if (rank >= nr_procs)
      color++;
  }

  if (nr_procs > size) {
    mpi_printf(mrc_mod_comm(mod), "ERROR: need %d procs, only have %d!\n",
	       nr_procs, size);
    abort();
  }

  if (color == nr_colors) { // unused proc
    color = MPI_UNDEFINED;
  }
  MPI_Comm_split(MPI_COMM_WORLD, color, rank, &mod->mod_comm);
  if (color == MPI_UNDEFINED) {
    return;
  }

  int c = 0;
  __list_for_each_entry(p, &mod->list, entry, struct mrc_mod_entry) {
    if (c++ == color)
      mod->this_mod_entry = p;
  }

  //  mprintf("color %d/%d %s\n", color, nr_colors, p->name);
}

void
mrc_mod_run(struct mrc_mod *mod)
{
  struct mrc_mod_entry *p = mod->this_mod_entry;
  if (p && p->func) {
    p->func(mod, p->arg);
  }
}

// ----------------------------------------------------------------------
// mrc_mod_register

void
mrc_mod_register(struct mrc_mod *mod, const char *name, int nr_procs,
		 void (*func)(struct mrc_mod *, void *), void *arg)
{
  struct mrc_mod_entry *entry = calloc(1, sizeof(*entry));
  entry->name = strdup(name);
  entry->nr_procs = nr_procs;
  entry->func = func;
  entry->arg = arg;

  list_add_tail(&entry->entry, &mod->list);
}

MPI_Comm
mrc_mod_get_comm(struct mrc_mod *mod)
{
  assert(mod->mod_comm);
  return mod->mod_comm;
}

static struct mrc_mod_entry *
find_mrc_mod_entry(struct mrc_mod *mod, const char *name)
{
  struct mrc_mod_entry *p;
  __list_for_each_entry(p, &mod->list, entry, struct mrc_mod_entry) {
    if (strcmp(p->name, name) == 0) {
      return p;
    }
  }
  return NULL;
}

int
mrc_mod_get_first_node(struct mrc_mod *mod, const char *name)
{
  struct mrc_mod_entry *p = find_mrc_mod_entry(mod, name);
  if (!p)
    return -1;

  return p->rank;
}

int
mrc_mod_get_nr_procs(struct mrc_mod *mod, const char *name)
{
  struct mrc_mod_entry *p = find_mrc_mod_entry(mod, name);
  if (!p)
    return -1;

  return p->nr_procs;
}

bool
mrc_mod_belongs_to(struct mrc_mod *mod, const char *name)
{
  struct mrc_mod_entry *p = mod->this_mod_entry;
  
  return p && (strcmp(p->name, name) == 0);
}

// ----------------------------------------------------------------------
// mrc_class_mrc_mod

struct mrc_class_mrc_mod mrc_class_mrc_mod = {
  .name         = "mrc_mod",
  .size         = sizeof(struct mrc_mod),
  .create       = _mrc_mod_create,
  .destroy      = _mrc_mod_destroy,
  .view         = _mrc_mod_view,
  .setup        = _mrc_mod_setup,
};

