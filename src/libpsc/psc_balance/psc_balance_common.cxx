
#include "fields.hxx"

#include <mrc_profile.h>
#include <string.h>

using Fields = Fields3d<fields_t, DIM_XYZ>;

static void
psc_balance_sub_communicate_particles(struct psc_balance *bal, struct communicate_ctx *ctx,
				      struct psc_mparticles *mprts_old, struct psc_mparticles *mprts_new,
				      int *nr_particles_by_patch_new)
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
    mparticles_patch_resize(mprts_new, p, nr_particles_by_patch_new[p]);
  }

  assert(sizeof(particle_t) % sizeof(particle_real_t) == 0); // FIXME

  // recv for new local patches
  MPI_Request *recv_reqs = new MPI_Request[ctx->nr_patches_new];
  int nr_recv_reqs = 0;

  for (int ri = 0; ri < ctx->nr_recv_ranks; ri++) {
    struct by_ri *recv = &ctx->recv_by_ri[ri];
    if (recv->rank == ctx->mpi_rank) {
      continue;
    }

    for (int pi = 0; pi < recv->nr_patches; pi++) {
      int p = recv->pi_to_patch[pi];
      particle_range_t prts_new = particle_range_mprts(mprts_new, p);
      int nn = particle_range_size(prts_new) * (sizeof(particle_t)  / sizeof(particle_real_t));
      MPI_Irecv(particle_iter_deref(prts_new.begin), nn, MPI_PARTICLES_REAL, recv->rank,
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
      particle_range_t prts_old = particle_range_mprts(mprts_old, p);
      int nn = particle_range_size(prts_old) * (sizeof(particle_t)  / sizeof(particle_real_t));
      //mprintf("A send -> %d tag %d (patch %d)\n", send->rank, pi, p);
      MPI_Isend(particle_iter_deref(prts_old.begin), nn, MPI_PARTICLES_REAL, send->rank,
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

    particle_range_t prts_old = particle_range_mprts(mprts_old, ctx->recv_info[p].patch);
    particle_range_t prts_new = particle_range_mprts(mprts_new, p);
    assert(particle_range_size(prts_old) == particle_range_size(prts_new));
#if 1
    for (int n = 0; n < particle_range_size(prts_new); n++) {
      *particle_iter_at(prts_new.begin, n) = *particle_iter_at(prts_old.begin, n);
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
  
  if (ctx->nr_patches_old != mflds_old->nr_patches /* || strncmp(flds->f[0].name, "lab", 3) == 0 */) return;
	
  assert(ctx->nr_patches_old == mflds_old->nr_patches);
  assert(ctx->nr_patches_old > 0);
  
  // send from old local patches
  MPI_Request *send_reqs = new MPI_Request[ctx->nr_patches_old]();
  int *nr_patches_new_by_rank = new int[ctx->mpi_size]();
  for (int p = 0; p < ctx->nr_patches_old; p++) {
    int new_rank = ctx->send_info[p].rank;
    if (new_rank == ctx->mpi_rank || new_rank < 0) {
      send_reqs[p] = MPI_REQUEST_NULL;
    } else {
      fields_t flds_old = fields_t_mflds(mflds_old, p);
      Fields F_old(flds_old);
      int nn = fields_t_size(flds_old) * flds_old.nr_comp;
      int *ib = flds_old.ib;
      void *addr_old = &F_old(0, ib[0], ib[1], ib[2]);
      int tag = nr_patches_new_by_rank[new_rank]++;
      MPI_Isend(addr_old, nn, MPI_FIELDS_REAL, new_rank, tag, ctx->comm, &send_reqs[p]);
    }
  }
  delete[] nr_patches_new_by_rank;

  // recv for new local patches
  MPI_Request *recv_reqs = new MPI_Request[ctx->nr_patches_new]();
  int *nr_patches_old_by_rank = new int[ctx->mpi_size];
  for (int p = 0; p < ctx->nr_patches_new; p++) {
    int old_rank = ctx->recv_info[p].rank;
    if (old_rank == ctx->mpi_rank) {
      recv_reqs[p] = MPI_REQUEST_NULL;
    } else if (old_rank < 0) { //this patch did not exist before
      recv_reqs[p] = MPI_REQUEST_NULL;
      //Seed new data
    } else {
      fields_t flds_new = fields_t_mflds(mflds_new, p);
      Fields F_new(flds_new);
      int nn = fields_t_size(flds_new) * flds_new.nr_comp;
      int *ib = flds_new.ib;
      void *addr_new = &F_new(0, ib[0], ib[1], ib[2]);
      int tag = nr_patches_old_by_rank[old_rank]++;
      MPI_Irecv(addr_new, nn, MPI_FIELDS_REAL, old_rank,
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

    fields_t flds_old = fields_t_mflds(mflds_old, ctx->recv_info[p].patch);
    fields_t flds_new = fields_t_mflds(mflds_new, p);
    Fields F_old(flds_old), F_new(flds_new);
    assert(flds_old.nr_comp == flds_new.nr_comp);
    assert(fields_t_size(flds_old) == fields_t_size(flds_new));
    int size = fields_t_size(flds_old) * flds_old.nr_comp;
    int *ib = flds_new.ib;
    void *addr_new = &F_new(0, ib[0], ib[1], ib[2]);
    void *addr_old = &F_old(0, ib[0], ib[1], ib[2]);
    memcpy(addr_new, addr_old, size * sizeof(fields_real_t));
  }
  prof_stop(pr);

  MPI_Waitall(ctx->nr_patches_old, send_reqs, MPI_STATUSES_IGNORE);
  MPI_Waitall(ctx->nr_patches_new, recv_reqs, MPI_STATUSES_IGNORE);
  delete[] send_reqs;
  delete[] recv_reqs;
}

