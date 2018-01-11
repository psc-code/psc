
#include "psc_push_particles_private.h"

#include <string.h>

struct sub {
  struct psc_push_particles *fwd;
};

#define sub(push) mrc_to_subobj(push, struct sub)

static void
psc_push_particles_1vb_setup(struct psc_push_particles *push)
{
  struct sub *sub = sub(push);

  const char *mprts_type = psc_mparticles_type(ppsc->particles);
  char s[10 + strlen(mprts_type)];
  if (strcmp(mprts_type, "single") == 0) {
    strcpy(s, "1vb_ps");
  } else if (strcmp(mprts_type, "cuda") == 0) {
    // FIXME, we there's a 4x4 blocksize default in psc_particles_cuda.c
    // which needs to be consistent
    strcpy(s, "1vb_4x4_cuda");
  } else {
    sprintf(s, "1vb_%s", mprts_type);
  }

  MPI_Comm comm = psc_push_particles_comm(push);
  mpi_printf(comm, "INFO: using psc_push_particles '%s'\n", s);
  
  sub->fwd = psc_push_particles_create(comm);
  psc_push_particles_set_type(sub->fwd, s);
  psc_push_particles_setup(sub->fwd);
  psc_push_particles_add_child(push, (struct mrc_obj *) sub->fwd);
}

static void
psc_push_particles_1vb_push_mprts_yz(struct psc_push_particles *push,
				     struct psc_mparticles *mprts,
				     struct psc_mfields *mflds)
{
  struct sub *sub = sub(push);
  struct psc_push_particles_ops *ops = psc_push_particles_ops(sub->fwd);
  assert(ops->push_mprts_yz);
  ops->push_mprts_yz(sub->fwd, mprts, mflds);
}

// ======================================================================
// psc_push_particles: subclass "1vb"

struct psc_push_particles_ops_1vb : psc_push_particles_ops {
  psc_push_particles_ops_1vb() {
    name                  = "1vb";
    size                  = sizeof(struct sub);
    setup                 = psc_push_particles_1vb_setup;
    push_mprts_yz         = psc_push_particles_1vb_push_mprts_yz;
  }
} psc_push_particles_1vb_ops;

