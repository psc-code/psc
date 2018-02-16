
#include "fields.hxx"

#include <mrc_profile.h>
#include <string.h>

using Fields = Fields3d<fields_t>;
using real_t = mparticles_t::real_t;

static void
psc_balance_sub_communicate_particles(struct psc_balance *bal, struct communicate_ctx *ctx,
				      struct psc_mparticles *mprts_old, struct psc_mparticles *mprts_new,
				      uint *nr_particles_by_patch_new)
{
  static int pr, pr_A, pr_B, pr_C, pr_D;
  if (!pr) {
    pr   = prof_register("comm prts", 1., 0, 0);
    pr_A = prof_register("comm prts A", 1., 0, 0);
    pr_B = prof_register("comm prts B", 1., 0, 0);
    pr_C = prof_register("comm prts C", 1., 0, 0);
    pr_D = prof_register("comm prts D", 1., 0, 0);
  }

  prof_start(pr);

  prof_start(pr_A);
  for (int p = 0; p < mprts_new->nr_patches; p++) {
    mparticles_t(mprts_new)[p].reserve(nr_particles_by_patch_new[p]);
    mparticles_t(mprts_new)[p].resize(nr_particles_by_patch_new[p]);
  }

  assert(sizeof(particle_t) % sizeof(real_t) == 0); // FIXME

  MPI_Datatype mpi_dtype = mparticles_traits<mparticles_t>::mpi_dtype();
  // recv for new local patches
  MPI_Request *recv_reqs = new MPI_Request[ctx->nr_patches_new]();
  int nr_recv_reqs = 0;

  for (int ri = 0; ri < ctx->nr_recv_ranks; ri++) {
    struct by_ri *recv = &ctx->recv_by_ri[ri];
    if (recv->rank == ctx->mpi_rank) {
      continue;
    }

    for (int pi = 0; pi < recv->nr_patches; pi++) {
      int p = recv->pi_to_patch[pi];
      mparticles_t::patch_t& prts_new = mparticles_t(mprts_new)[p];
      int nn = prts_new.size() * (sizeof(particle_t)  / sizeof(real_t));
      MPI_Irecv(&*prts_new.begin(), nn, mpi_dtype, recv->rank,
		pi, ctx->comm, &recv_reqs[nr_recv_reqs++]);
    }
  }
  prof_stop(pr_A);

  prof_start(pr_B);
  // send from old local patches
  MPI_Request *send_reqs = new MPI_Request[ctx->nr_patches_old]();
  int nr_send_reqs = 0;

  for (int ri = 0; ri < ctx->nr_send_ranks; ri++) {
    struct by_ri *send = &ctx->send_by_ri[ri];
    if (send->rank == ctx->mpi_rank) {
      continue;
    }

    for (int pi = 0; pi < send->nr_patches; pi++) {
      int p = send->pi_to_patch[pi];
      mparticles_t::patch_t& prts_old = mparticles_t(mprts_old)[p];
      int nn = prts_old.size() * (sizeof(particle_t)  / sizeof(real_t));
      //mprintf("A send -> %d tag %d (patch %d)\n", send->rank, pi, p);
      MPI_Isend(&*prts_old.begin(), nn, mpi_dtype, send->rank,
      		pi, ctx->comm, &send_reqs[nr_send_reqs++]);
    }
  }
  prof_stop(pr_B);

  prof_start(pr_C);
  // local particles
  // OPT: could keep the alloced arrays, just move pointers...
  for (int p = 0; p < ctx->nr_patches_new; p++) {
    if (ctx->recv_info[p].rank != ctx->mpi_rank) {
      continue;
    }

    mparticles_t::patch_t& prts_old = mparticles_t(mprts_old)[ctx->recv_info[p].patch];
    mparticles_t::patch_t& prts_new = mparticles_t(mprts_new)[p];
    assert(prts_old.size() == prts_new.size());
#if 1
    for (int n = 0; n < prts_new.size(); n++) {
      prts_new[n] = prts_old[n];
    }
#else
    // FIXME, this needs at least a proper interface -- if not separately alloc'ed, bad things
    // are going to happen
    free(c_new->particles);
    c_new->particles = c_old->particles;
    c_new->n_alloced = c_old->n_alloced;
    c_old->particles = NULL; // prevent from being freed
    c_old->n_alloced = 0;
#endif
  }
  prof_stop(pr_C);
  
  prof_start(pr_D);
  MPI_Waitall(nr_send_reqs, send_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(nr_recv_reqs, recv_reqs, MPI_STATUSES_IGNORE);
  delete[] send_reqs;
  delete[] recv_reqs;
  prof_stop(pr_D);

  prof_stop(pr);
}

static void
psc_balance_sub_communicate_fields(struct psc_balance *bal, struct communicate_ctx *ctx,
				   struct psc_mfields *mflds_old, struct psc_mfields *mflds_new)
{
  //HACK: Don't communicate output fields if they don't correspond to the domain
  //This is needed e.g. for the boosted output which handles its MPI communication internally
  //printf("Field: %s\n", flds->f[0].name);

  mfields_t mf_old(mflds_old), mf_new(mflds_new);
  
  if (ctx->nr_patches_old != mflds_old->nr_patches /* || strncmp(flds->f[0].name, "lab", 3) == 0 */) return;
	
  assert(ctx->nr_patches_old == mflds_old->nr_patches);
  assert(ctx->nr_patches_old > 0);
  
  MPI_Datatype mpi_dtype = fields_traits<fields_t>::mpi_dtype();

  // send from old local patches
  MPI_Request *send_reqs = new MPI_Request[ctx->nr_patches_old]();
  int *nr_patches_new_by_rank = new int[ctx->mpi_size]();
  for (int p = 0; p < ctx->nr_patches_old; p++) {
    int new_rank = ctx->send_info[p].rank;
    if (new_rank == ctx->mpi_rank || new_rank < 0) {
      send_reqs[p] = MPI_REQUEST_NULL;
    } else {
      fields_t flds_old = mf_old[p];
      Fields F_old(flds_old);
      int nn = flds_old.size();
      int *ib = flds_old.ib;
      void *addr_old = &F_old(0, ib[0], ib[1], ib[2]);
      int tag = nr_patches_new_by_rank[new_rank]++;
      MPI_Isend(addr_old, nn, mpi_dtype, new_rank, tag, ctx->comm, &send_reqs[p]);
    }
  }
  delete[] nr_patches_new_by_rank;

  // recv for new local patches
  MPI_Request *recv_reqs = new MPI_Request[ctx->nr_patches_new]();
  int *nr_patches_old_by_rank = new int[ctx->mpi_size]();
  for (int p = 0; p < ctx->nr_patches_new; p++) {
    int old_rank = ctx->recv_info[p].rank;
    if (old_rank == ctx->mpi_rank) {
      recv_reqs[p] = MPI_REQUEST_NULL;
    } else if (old_rank < 0) { //this patch did not exist before
      recv_reqs[p] = MPI_REQUEST_NULL;
      //Seed new data
    } else {
      fields_t flds_new = mf_new[p];
      Fields F_new(flds_new);
      int nn = flds_new.size();
      int *ib = flds_new.ib;
      void *addr_new = &F_new(0, ib[0], ib[1], ib[2]);
      int tag = nr_patches_old_by_rank[old_rank]++;
      MPI_Irecv(addr_new, nn, mpi_dtype, old_rank,
		tag, ctx->comm, &recv_reqs[p]);
    }
  }
  delete[] nr_patches_old_by_rank;

  static int pr;
  if (!pr) {
    pr = prof_register("bal flds local", 1., 0, 0);
  }

  prof_start(pr);
  // local fields
  // OPT: could keep the alloced arrays, just move pointers...
  for (int p = 0; p < ctx->nr_patches_new; p++) {
    if (ctx->recv_info[p].rank != ctx->mpi_rank) {
      continue;
    }

    fields_t flds_old = mf_old[ctx->recv_info[p].patch];
    fields_t flds_new = mf_new[p];
    Fields F_old(flds_old), F_new(flds_new);
    assert(flds_old.nr_comp == flds_new.nr_comp);
    assert(flds_old.size() == flds_new.size());
    int size = flds_old.size();
    int *ib = flds_new.ib;
    void *addr_new = &F_new(0, ib[0], ib[1], ib[2]);
    void *addr_old = &F_old(0, ib[0], ib[1], ib[2]);
    memcpy(addr_new, addr_old, size * sizeof(fields_t::real_t));
  }
  prof_stop(pr);

  MPI_Waitall(ctx->nr_patches_old, send_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(ctx->nr_patches_new, recv_reqs, MPI_STATUSES_IGNORE);
  delete[] send_reqs;
  delete[] recv_reqs;
}

