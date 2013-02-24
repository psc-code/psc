
#ifndef MRC_MOD_H
#define MRC_MOD_H

#include <mrc_obj.h>

MRC_CLASS_DECLARE(mrc_mod, struct mrc_mod);

void mrc_mod_register(struct mrc_mod *mod, const char *name, int nr_procs,
		      void (*func)(struct mrc_mod *, void *), void *arg);
void mrc_mod_run(struct mrc_mod *mod);
int  mrc_mod_get_first_node(struct mrc_mod *mod, const char *name);
int  mrc_mod_get_nr_procs(struct mrc_mod *mod, const char *name);
bool mrc_mod_belongs_to(struct mrc_mod *mod, const char *name);
MPI_Comm mrc_mod_get_comm(struct mrc_mod *mod);

#endif

