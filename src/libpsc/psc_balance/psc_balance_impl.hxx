
#pragma once

#include "fields.hxx"
#include "fields_traits.hxx"
#include "balance.hxx"
#include "bnd_particles.hxx"
#include "bnd.hxx"
#include "mpi_dtype_traits.hxx"

#ifdef USE_CUDA
#include "psc_fields_cuda.h"
#include "mparticles_cuda.hxx"
#include "bnd_cuda_3_impl.hxx"
#endif

#include <mrc_profile.h>
#include <string.h>

#include <numeric>

extern double* psc_balance_comp_time_by_patch;

static double capability_default(int p) { return 1.; }

// ======================================================================

namespace psc
{
namespace balance
{

struct input
{
  std::vector<double> capability_by_proc;
  std::vector<double> load_by_patch;

  int n_procs() const { return capability_by_proc.size(); }

  double total_capability() const
  {
    return std::accumulate(capability_by_proc.begin(), capability_by_proc.end(),
                           0.);
  }

  double total_load() const
  {
    return std::accumulate(load_by_patch.begin(), load_by_patch.end(), 0.);
  }

  double load_target() const { return total_load() / total_capability(); }
};

#if 0
inline std::vector<int> best_mapping(const input& input)
{
  int n_procs = input.n_procs();
  std::vector<int> nr_patches_all_new(n_procs);

  auto& capability = input.capability_by_proc;
  auto& loads_all = input.load_by_patch;

  double load_target = input.load_target();
  mprintf("psc_balance: loads_sum %g capability_sum %g load_target %g\n",
          input.total_load(), input.total_capability(), load_target);
  int p = 0, nr_new_patches = 0;
  double load = 0.;
  double next_target = load_target * capability[0];
  for (int i = 0; i < loads_all.size(); i++) {
    load += loads_all[i];
    nr_new_patches++;
    if (p < n_procs - 1) {
      // if load limit is reached, or we have only as many patches as
      // processors left
      if (load > next_target || n_procs - p >= loads_all.size() - i) {
        double above_target = load - next_target;
        double below_target = next_target - (load - loads_all[i]);
        if (above_target > below_target && nr_new_patches > 1) {
          nr_patches_all_new[p] = nr_new_patches - 1;
          nr_new_patches = 1;
        } else {
          nr_patches_all_new[p] = nr_new_patches;
          nr_new_patches = 0;
        }
        p++;
        next_target += load_target * capability[p];
      }
    }
    // last proc takes what's left
    if (i == loads_all.size() - 1) {
      nr_patches_all_new[n_procs - 1] = nr_new_patches;
    }
  }

  return nr_patches_all_new;
}

#else

inline void best_mapping_recursive(const input& input,
                                   std::vector<int>& n_patches_by_proc,
                                   int proc_begin, int proc_end,
                                   int patch_begin, int patch_end)
{
  // we always must have at least one patch per proc
  assert(patch_end - patch_begin >= proc_end - proc_begin);

  if (proc_end == proc_begin + 1) {
    n_patches_by_proc[proc_begin] = patch_end - patch_begin;
  } else {
    int proc_middle = (proc_begin + proc_end) / 2;
    double load_total = std::accumulate(&input.load_by_patch[patch_begin],
                                        &input.load_by_patch[patch_end], 0.);
    double load_target =
      (load_total * (proc_middle - proc_begin)) / (proc_end - proc_begin);
    std::cout << "[" << proc_begin << ":" << proc_end << "] total "
              << load_total << " target " << load_target << "\n";
    double load = 0.;
    int patch_middle = patch_begin;
    for (;;) {
      double prev_load = load;
      load += input.load_by_patch[patch_middle];
      patch_middle++;
      if (load > load_target) {
        double above = load - load_target;
        double below = load_target - prev_load;
        // std::cout << "above " << above << " below " << below << "\n";

        // if we were closer to the target load one step ago, go with that
        if (below < above && patch_middle) {
          patch_middle--;
          load = prev_load;
        }
        break;
      }
    }
    std::cout << "middle " << patch_middle << " load " << load << "\n";

    // make sure there are at least as many patches as procs on either side
    if (patch_middle - patch_begin < proc_middle - proc_begin) {
      patch_middle = patch_begin + (proc_middle - proc_begin);
    }
    if (patch_end - patch_middle < proc_end - proc_middle) {
      patch_middle = patch_end - (proc_end - proc_middle);
    }

    best_mapping_recursive(input, n_patches_by_proc, proc_begin, proc_middle,
                           patch_begin, patch_middle);
    best_mapping_recursive(input, n_patches_by_proc, proc_middle, proc_end,
                           patch_middle, patch_end);
  }
}

inline std::vector<int> best_mapping(const input& input)
{
  int n_procs = input.n_procs();
  std::vector<int> n_patches_by_proc(n_procs);
  best_mapping_recursive(input, n_patches_by_proc, 0, n_procs, 0,
                         input.load_by_patch.size());

  return n_patches_by_proc;
}

#endif
inline void print_stats(const input& input,
                        const std::vector<int>& nr_patches_all_new,
                        bool verbose)
{
  auto& capability = input.capability_by_proc;
  auto& loads_all = input.load_by_patch;

  int n_procs = input.n_procs();
  double load_target = input.load_target();

  int pp = 0;
  double min_diff = 0, max_diff = 0;
  for (int p = 0; p < n_procs; p++) {
    double load = 0.;
    for (int i = 0; i < nr_patches_all_new[p]; i++) {
      load += loads_all[pp++];
      if (verbose) {
        mprintf("  pp %d load %g : %g\n", pp - 1, loads_all[pp - 1], load);
      }
    }
    double diff = load - load_target * capability[p];
    if (verbose) {
      mprintf("p %d # = %d load %g / %g : diff %g %%\n", p,
              nr_patches_all_new[p], load, load_target * capability[p],
              100. * diff / (load_target * capability[p]));
    }
    if (diff < min_diff) {
      min_diff = diff;
    }
    if (diff > max_diff) {
      max_diff = diff;
    }
  }
  mprintf("psc_balance: achieved target %g (%g %% -- %g %%)\n", load_target,
          100 * min_diff / load_target, 100 * max_diff / load_target);
}

inline void write_loads(const input& input,
                        const std::vector<int>& nr_patches_all_new,
                        int timestep)
{
  char s[20];
  snprintf(s, 20, "loads2-%06d.asc", timestep);
  FILE* f = fopen(s, "w");

  int gp = 0;
  for (int r = 0; r < nr_patches_all_new.size(); r++) {
    for (int p = 0; p < nr_patches_all_new[r]; p++) {
      fprintf(f, "%d %g %d\n", gp, input.load_by_patch[gp], r);
      gp++;
    }
  }

  fclose(f);
}

// ======================================================================
// get_loads

class get_loads
{
public:
  get_loads(double factor_fields) : factor_fields_{factor_fields} {}

  // given the number of particles for each patch (even though they don't
  // actually exist yet), figure out the load for each patch
  std::vector<double> initial(const Grid_t& grid,
                              const std::vector<uint>& n_prts_by_patch)
  {
    assert(n_prts_by_patch.size() == grid.n_patches());
    std::vector<double> loads;
    loads.reserve(n_prts_by_patch.size());

    auto& ldims = grid.ldims;
    for (auto n_prts : n_prts_by_patch) {
      loads.push_back(n_prts + factor_fields_ * ldims[0] * ldims[1] * ldims[2]);
    }

    return loads;
  }

  std::vector<double> operator()(const Grid_t& grid, MparticlesBase& mprts)
  {
    auto n_prts_by_patch = mprts.sizeByPatch();

    std::vector<double> loads;
    loads.reserve(mprts.n_patches());
    for (int p = 0; p < mprts.n_patches(); p++) {
      double load;
      if (factor_fields_ >= 0.) {
        const int* ldims = grid.ldims;
        load =
          n_prts_by_patch[p] + factor_fields_ * ldims[0] * ldims[1] * ldims[2];
        // mprintf("loads p %d %g %g ratio %g\n", p, loads[p], comp_time,
        // loads[p] / comp_time);
      } else {
        load = psc_balance_comp_time_by_patch[p];
      }
#if 0
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      load = 1 + rank;
#endif
      loads.push_back(load);
    }
    return loads;
  }

private:
  double factor_fields_;
};

inline std::vector<double> gather_loads(const Grid_t& grid,
                                        std::vector<double> loads)
{
  const MrcDomain& domain = grid.mrc_domain_;
  MPI_Comm comm = domain.comm();
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // gather nr_patches for all procs on proc 0
  std::vector<int> nr_patches_all;
  if (rank == 0) {
    nr_patches_all.resize(size);
  }
  int n_patches = loads.size();
  MPI_Gather(&n_patches, 1, MPI_INT, nr_patches_all.data(), 1, MPI_INT, 0,
             comm);

  // gather loads for all patches on proc 0
  std::vector<int> displs;
  std::vector<double> loads_all;
  if (rank == 0) {
    displs.resize(size);
    int off = 0;
    for (int i = 0; i < size; i++) {
      displs[i] = off;
      off += nr_patches_all[i];
    }
    int n_global_patches = domain.nGlobalPatches();
    loads_all.resize(n_global_patches);
  }
  MPI_Gatherv(loads.data(), n_patches, MPI_DOUBLE, loads_all.data(),
              nr_patches_all.data(), displs.data(), MPI_DOUBLE, 0, comm);

  if (rank == 0) {
#if 0
      int rl = 0;
      char s[20]; sprintf(s, "loads-%06d.asc", grid.timestep());
      FILE *f = fopen(s, "w");
      for (int p = 0; p < *p_nr_global_patches; p++) {
	if (rl < size - 1 && p >= displs[rl+1]) {
	  rl++;
	}
	fprintf(f, "%d %g %d\n", p, loads_all[p], rl);
      }
      fclose(f);
#endif
  }

  return loads_all;
}

} // namespace balance
} // namespace psc

// ======================================================================
// Communicate

struct send_info
{
  int rank;
  int patch;
};

struct recv_info
{
  int rank;
  int patch;
};

struct by_ri
{
  int rank;
  int nr_patches;
  std::vector<int> pi_to_patch;
};

class communicate_ctx
{
public:
  communicate_ctx(const MrcDomain& domain_old, const MrcDomain& domain_new)
  {
    comm_ = domain_old.comm();
    MPI_Comm_rank(comm_, &mpi_rank_);
    MPI_Comm_size(comm_, &mpi_size_);

    nr_patches_old_ = domain_old.nPatches();
    nr_patches_new_ = domain_new.nPatches();

    send_info_.resize(nr_patches_old_);
    for (int p = 0; p < nr_patches_old_; p++) {
      auto info_old = domain_old.localPatchInfo(p);
      auto info_new =
        domain_new.levelIdx3PatchInfo(info_old.level, info_old.idx3);
      send_info_[p].rank = info_new.rank;
      send_info_[p].patch = info_new.patch;
    }

    recv_info_.resize(nr_patches_new_);
    for (int p = 0; p < nr_patches_new_; p++) {
      auto info_new = domain_new.localPatchInfo(p);
      auto info_old =
        domain_old.levelIdx3PatchInfo(info_new.level, info_new.idx3);
      recv_info_[p].rank = info_old.rank;
      recv_info_[p].patch = info_old.patch;
    }

    // maps rank <-> rank index

    send_rank_to_ri_.resize(mpi_size_);
    recv_rank_to_ri_.resize(mpi_size_);
    for (int r = 0; r < mpi_size_; r++) {
      send_rank_to_ri_[r] = -1;
      recv_rank_to_ri_[r] = -1;
    }

    nr_send_ranks_ = 0;
    for (int p = 0; p < nr_patches_old_; p++) {
      int send_rank = send_info_[p].rank;
      if (send_rank >= 0) {
        if (send_rank_to_ri_[send_rank] < 0) {
          send_rank_to_ri_[send_rank] = nr_send_ranks_++;
        }
      }
    }
    nr_recv_ranks_ = 0;
    for (int p = 0; p < nr_patches_new_; p++) {
      int recv_rank = recv_info_[p].rank;
      if (recv_rank >= 0) {
        if (recv_rank_to_ri_[recv_rank] < 0) {
          recv_rank_to_ri_[recv_rank] = nr_recv_ranks_++;
        }
      }
    }

    send_by_ri_.resize(nr_send_ranks_);
    recv_by_ri_.resize(nr_recv_ranks_);

    for (int p = 0; p < nr_patches_old_; p++) {
      int send_rank = send_info_[p].rank;
      if (send_rank >= 0) {
        send_by_ri_[send_rank_to_ri_[send_rank]].rank = send_rank;
      }
    }
    for (int p = 0; p < nr_patches_new_; p++) {
      int recv_rank = recv_info_[p].rank;
      if (recv_rank >= 0) {
        recv_by_ri_[recv_rank_to_ri_[recv_rank]].rank = recv_rank;
      }
    }

    /* for (int ri = 0; ri < nr_send_ranks_; ri++) { */
    /*   mprintf("send -> %d (%d)\n", send_by_ri_[ri].rank, ri); */
    /* } */

    // count number of patches sent to each rank

    for (int p = 0; p < nr_patches_old_; p++) {
      int send_rank = send_info_[p].rank;
      if (send_rank >= 0) {
        send_by_ri_[send_rank_to_ri_[send_rank]].nr_patches++;
      }
    }

    // map send patch index by ri back to local patch number

    for (int ri = 0; ri < nr_send_ranks_; ri++) {
      send_by_ri_[ri].pi_to_patch.resize(send_by_ri_[ri].nr_patches);
      send_by_ri_[ri].nr_patches = 0;
    }

    for (int p = 0; p < nr_patches_old_; p++) {
      int send_rank = send_info_[p].rank;
      if (send_rank < 0) {
        continue;
      }
      int ri = send_rank_to_ri_[send_rank];
      int pi = send_by_ri_[ri].nr_patches++;
      send_by_ri_[ri].pi_to_patch[pi] = p;
    }

    // count number of patches received from each rank

    for (int p = 0; p < nr_patches_new_; p++) {
      int recv_rank = recv_info_[p].rank;
      if (recv_rank >= 0) {
        int ri = recv_rank_to_ri_[recv_rank];
        recv_by_ri_[ri].nr_patches++;
      }
    }

    // map received patch index by ri back to local patch number

    for (int ri = 0; ri < nr_recv_ranks_; ri++) {
      recv_by_ri_[ri].pi_to_patch.resize(recv_by_ri_[ri].nr_patches);
      recv_by_ri_[ri].nr_patches = 0;
    }

    for (int p = 0; p < nr_patches_new_; p++) {
      int recv_rank = recv_info_[p].rank;
      if (recv_rank < 0) {
        continue;
      }
      int ri = recv_rank_to_ri_[recv_rank];
      int pi = recv_by_ri_[ri].nr_patches++;
      recv_by_ri_[ri].pi_to_patch[pi] = p;
    }
  }

  std::vector<uint> new_n_prts(std::vector<uint> n_prts_by_patch_old)
  {
    static int pr;
    if (!pr) {
      pr = prof_register("comm nr prts", 1., 0, 0);
    }

    prof_start(pr);

    std::vector<uint> n_prts_by_patch_new(nr_patches_new_);
    // post receives

    std::vector<MPI_Request> recv_reqs(nr_patches_new_);
    int nr_recv_reqs = 0;

    std::vector<std::vector<int>> nr_particles_recv_by_ri(nr_recv_ranks_);
    for (int ri = 0; ri < nr_recv_ranks_; ri++) {
      struct by_ri* recv = &recv_by_ri_[ri];
      nr_particles_recv_by_ri[ri].resize(recv->nr_patches);

      if (recv->rank != mpi_rank_) {
        // mprintf("recv <- %d (len %d)\n", r, nr_patches_recv_by_ri[ri]);
        MPI_Irecv(nr_particles_recv_by_ri[ri].data(), recv->nr_patches, MPI_INT,
                  recv->rank, 10, comm_, &recv_reqs[nr_recv_reqs++]);
      }
    }

    // post sends

    std::vector<MPI_Request> send_reqs(nr_send_ranks_);
    int nr_send_reqs = 0;

    std::vector<std::vector<int>> nr_particles_send_by_ri(nr_send_ranks_);
    for (int ri = 0; ri < nr_send_ranks_; ri++) {
      struct by_ri* send = &send_by_ri_[ri];
      nr_particles_send_by_ri[ri].resize(send->nr_patches);

      for (int pi = 0; pi < send->nr_patches; pi++) {
        nr_particles_send_by_ri[ri][pi] =
          n_prts_by_patch_old[send->pi_to_patch[pi]];
      }

      if (send->rank != mpi_rank_) {
        // mprintf("send -> %d (len %d)\n", r, nr_patches_send_by_ri[ri]);
        MPI_Isend(nr_particles_send_by_ri[ri].data(), send->nr_patches, MPI_INT,
                  send->rank, 10, comm_, &send_reqs[nr_send_reqs++]);
      }
    }
    assert(nr_send_reqs <= nr_send_ranks_);

    // copy local particle numbers
    {
      int send_ri = send_rank_to_ri_[mpi_rank_];
      int recv_ri = recv_rank_to_ri_[mpi_rank_];
      if (send_ri < 0) { // no local patches to copy
        assert(recv_ri < 0);
      } else {
        assert(send_by_ri_[send_ri].nr_patches ==
               recv_by_ri_[recv_ri].nr_patches);
        for (int n = 0; n < send_by_ri_[send_ri].nr_patches; n++) {
          nr_particles_recv_by_ri[recv_ri][n] =
            nr_particles_send_by_ri[send_ri][n];
        }
      }
    }

    MPI_Waitall(nr_recv_reqs, recv_reqs.data(), MPI_STATUSES_IGNORE);

    // update from received data

    for (int ri = 0; ri < nr_recv_ranks_; ri++) {
      struct by_ri* recv = &recv_by_ri_[ri];
      for (int pi = 0; pi < recv_by_ri_[ri].nr_patches; pi++) {
        n_prts_by_patch_new[recv->pi_to_patch[pi]] =
          nr_particles_recv_by_ri[ri][pi];
      }
    }

    MPI_Waitall(nr_send_reqs, send_reqs.data(), MPI_STATUSES_IGNORE);

    // return result

    prof_stop(pr);
    return n_prts_by_patch_new;
  }

  template <typename Mparticles>
  void communicate_particles(Mparticles& mp_old, Grid_t& new_grid)
  {
    using Particle = typename Mparticles::Particle;
    using real_t = typename Mparticles::real_t;
    // static int pr, pr_A, pr_B, pr_C, pr_D;
    // if (!pr) {
    //   pr   = prof_register("comm prts", 1., 0, 0);
    //   pr_A = prof_register("comm prts A", 1., 0, 0);
    //   pr_B = prof_register("comm prts B", 1., 0, 0);
    //   pr_C = prof_register("comm prts C", 1., 0, 0);
    //   pr_D = prof_register("comm prts D", 1., 0, 0);
    // }

    // prof_start(pr);

    auto mp_new = Mparticles{new_grid};

    auto n_prts_by_patch_new = new_n_prts(mp_old.sizeByPatch());

    // prof_start(pr_A);
    mp_new.reserve_all(n_prts_by_patch_new);
    mp_new.resize_all(n_prts_by_patch_new);

    assert(sizeof(Particle) % sizeof(real_t) == 0); // FIXME

    auto mpi_dtype = MpiDtypeTraits<typename Mparticles::real_t>::value();
    // recv for new local patches
    std::vector<MPI_Request> recv_reqs(nr_patches_new_);
    int nr_recv_reqs = 0;

    for (int ri = 0; ri < nr_recv_ranks_; ri++) {
      struct by_ri* recv = &recv_by_ri_[ri];
      if (recv->rank == mpi_rank_) {
        continue;
      }

      for (int pi = 0; pi < recv->nr_patches; pi++) {
        int p = recv->pi_to_patch[pi];
        auto&& prts_new = mp_new[p];
        int nn = mp_new.size(p) * (sizeof(Particle) / sizeof(real_t));
        MPI_Irecv(&*prts_new.begin(), nn, mpi_dtype, recv->rank, pi, comm_,
                  &recv_reqs[nr_recv_reqs++]);
      }
    }
    // prof_stop(pr_A);

    // prof_start(pr_B);
    // send from old local patches
    std::vector<MPI_Request> send_reqs(nr_patches_old_);
    int nr_send_reqs = 0;

    for (int ri = 0; ri < nr_send_ranks_; ri++) {
      struct by_ri* send = &send_by_ri_[ri];
      if (send->rank == mpi_rank_) {
        continue;
      }

      for (int pi = 0; pi < send->nr_patches; pi++) {
        int p = send->pi_to_patch[pi];
        auto&& prts_old = mp_old[p];
        int nn = mp_old.size(p) * (sizeof(Particle) / sizeof(real_t));
        // mprintf("A send -> %d tag %d (patch %d)\n", send->rank, pi, p);
        MPI_Isend(&*prts_old.begin(), nn, mpi_dtype, send->rank, pi, comm_,
                  &send_reqs[nr_send_reqs++]);
      }
    }
    // prof_stop(pr_B);

    // prof_start(pr_C);
    // local particles
    // OPT: could keep the alloced arrays, just move pointers...
    for (int p = 0; p < nr_patches_new_; p++) {
      if (recv_info_[p].rank != mpi_rank_) {
        continue;
      }

      auto&& prts_old = mp_old[recv_info_[p].patch];
      auto&& prts_new = mp_new[p];
      assert(mp_old.size(recv_info_[p].patch) == mp_new.size(p));
#if 1
      for (int n = 0; n < mp_new.size(p); n++) {
        prts_new[n] = prts_old[n];
      }
#else
      // FIXME, this needs at least a proper interface -- if not separately
      // alloc'ed, bad things are going to happen
      free(c_new->particles);
      c_new->particles = c_old->particles;
      c_new->n_alloced = c_old->n_alloced;
      c_old->particles = NULL; // prevent from being freed
      c_old->n_alloced = 0;
#endif
    }
    // prof_stop(pr_C);

    // prof_start(pr_D);
    MPI_Waitall(nr_send_reqs, send_reqs.data(), MPI_STATUSES_IGNORE);
    MPI_Waitall(nr_recv_reqs, recv_reqs.data(), MPI_STATUSES_IGNORE);
    // prof_stop(pr_D);

    mp_old = std::move(mp_new);
    // prof_stop(pr);
  }

  template <typename Mfields>
  void communicate_fields(Mfields& mf_old, Grid_t& new_grid)
  {
    MPI_Datatype mpi_dtype = Mfields_traits<Mfields>::mpi_dtype();

    mpi_printf(new_grid.comm(),
               "***** Balance: balancing field, #components %d\n",
               mf_old.n_comps());

    auto mf_new = Mfields{new_grid, mf_old.n_comps(), mf_old.ibn()};

    // send from old local patches
    std::vector<MPI_Request> send_reqs(nr_patches_old_);
    std::vector<int> nr_patches_new_by_rank(mpi_size_);
    for (int p = 0; p < nr_patches_old_; p++) {
      int new_rank = send_info_[p].rank;
      if (new_rank == mpi_rank_ || new_rank < 0) {
        send_reqs[p] = MPI_REQUEST_NULL;
      } else {
        auto flds_old = mf_old[p];
        int nn = flds_old.storage().size();
        void* addr_old = flds_old.storage().data();
        int tag = nr_patches_new_by_rank[new_rank]++;
        MPI_Isend(addr_old, nn, mpi_dtype, new_rank, tag, comm_, &send_reqs[p]);
      }
    }

    // recv for new local patches
    std::vector<MPI_Request> recv_reqs(nr_patches_new_);
    std::vector<int> nr_patches_old_by_rank(mpi_size_);
    for (int p = 0; p < nr_patches_new_; p++) {
      int old_rank = recv_info_[p].rank;
      if (old_rank == mpi_rank_) {
        recv_reqs[p] = MPI_REQUEST_NULL;
      } else if (old_rank < 0) { // this patch did not exist before
        recv_reqs[p] = MPI_REQUEST_NULL;
        // Seed new data
      } else {
        auto flds_new = mf_new[p];
        int nn = flds_new.storage().size();
        void* addr_new = flds_new.storage().data();
        int tag = nr_patches_old_by_rank[old_rank]++;
        MPI_Irecv(addr_new, nn, mpi_dtype, old_rank, tag, comm_, &recv_reqs[p]);
      }
    }

    static int pr;
    if (!pr) {
      pr = prof_register("bal flds local", 1., 0, 0);
    }

    prof_start(pr);
    // local fields
    // OPT: could keep the alloced arrays, just move pointers...
    for (int p = 0; p < nr_patches_new_; p++) {
      if (recv_info_[p].rank != mpi_rank_) {
        continue;
      }

      auto flds_old = mf_old[recv_info_[p].patch];
      auto flds_new = mf_new[p];
      assert(flds_old.storage().shape() == flds_new.storage().shape());
      int size = flds_old.storage().size();
      void* addr_new = flds_new.storage().data();
      void* addr_old = flds_old.storage().data();
      memcpy(addr_new, addr_old, size * sizeof(typename Mfields::real_t));
    }
    prof_stop(pr);

    MPI_Waitall(nr_patches_old_, send_reqs.data(), MPI_STATUSES_IGNORE);
    MPI_Waitall(nr_patches_new_, recv_reqs.data(), MPI_STATUSES_IGNORE);

    mf_old = std::move(mf_new);
  }

private:
  MPI_Comm comm_;
  int mpi_rank_;
  int mpi_size_;
  int nr_patches_old_;
  int nr_patches_new_;
  std::vector<send_info> send_info_; // by old patch on this proc
  std::vector<recv_info> recv_info_; // by new patch on this proc
  std::vector<int>
    send_rank_to_ri_; // map from send target rank to contiguous "rank index"
  std::vector<int>
    recv_rank_to_ri_; // map from recv source rank to contiguous "rank index"

  int nr_send_ranks_;
  std::vector<by_ri> send_by_ri_;
  int nr_recv_ranks_;
  std::vector<by_ri> recv_by_ri_;
};

// ======================================================================
// Balance_

template <typename Mparticles, typename MfieldsState, typename Mfields>
struct Balance_
{
  using Particle = typename Mparticles::Particle;
  using real_t = typename Mparticles::real_t;

  Balance_(double factor_fields = 1., bool print_loads = false,
           bool write_loads = false)
    : factor_fields_(factor_fields),
      get_loads_{factor_fields},
      print_loads_(print_loads),
      write_loads_(write_loads)
  {}

  ~Balance_()
  {
    delete[] psc_balance_comp_time_by_patch;
    psc_balance_comp_time_by_patch = nullptr;
  }

  void initial(Grid_t*& grid, std::vector<uint>& n_prts_by_patch)
  {
    auto loads = get_loads_.initial(*grid, n_prts_by_patch);
    auto loads_all = psc::balance::gather_loads(*grid, loads);
    n_prts_by_patch = balance(grid, loads_all, nullptr, n_prts_by_patch);
  }

  void operator()(Grid_t*& grid, MparticlesBase& mprts)
  {
    static int st_time_balance;
    if (!st_time_balance) {
      st_time_balance = psc_stats_register("time balancing");
    }

    psc_stats_start(st_time_balance);
    auto loads = get_loads_(mprts.grid(), mprts);
    auto loads_all = psc::balance::gather_loads(*grid, loads);
    balance(grid, loads_all, &mprts);
    psc_stats_stop(st_time_balance);
  }

private:
  int find_best_mapping(const Grid_t& grid,
                        const std::vector<double>& loads_all)
  {
    const MrcDomain& domain = grid.mrc_domain_;

    MPI_Comm comm = domain.comm();
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int nr_patches_old = domain.nPatches();

    std::vector<int> nr_patches_all_old(size);
    MPI_Allgather(&nr_patches_old, 1, MPI_INT, nr_patches_all_old.data(), 1,
                  MPI_INT, comm);

    int nr_patches_new;

    if (rank == 0) { // do the mapping on proc 0
      std::vector<double> capability(size);
      for (int p = 0; p < size; p++) {
        capability[p] = capability_default(p);
      }

      auto nr_patches_all_new =
        psc::balance::best_mapping({capability, loads_all});
      psc::balance::print_stats({capability, loads_all}, nr_patches_all_new,
                                print_loads_);
      if (write_loads_) {
        psc::balance::write_loads({capability, loads_all}, nr_patches_all_new,
                                  grid.timestep());
      }

      if (nr_patches_all_new == nr_patches_all_old) {
        std::fill(nr_patches_all_new.begin(), nr_patches_all_new.end(),
                  -1); // unchanged mapping, no communication etc needed
      }

      MPI_Scatter(nr_patches_all_new.data(), 1, MPI_INT, &nr_patches_new, 1,
                  MPI_INT, 0, comm);
    } else {
      MPI_Scatter(nullptr, 1, MPI_INT, &nr_patches_new, 1, MPI_INT, 0, comm);
    }
    return nr_patches_new;
  }

  std::vector<uint> balance(Grid_t*& gridp, std::vector<double> loads_all,
                            MparticlesBase* mp,
                            std::vector<uint> n_prts_by_patch_old = {})
  {
    static int pr_bal_load, pr_bal_ctx, pr_bal_prts, pr_bal_flds;
    if (!pr_bal_load) {
      pr_bal_load = prof_register("bal load", 1., 0, 0);
      pr_bal_ctx = prof_register("bal ctx", 1., 0, 0);
      pr_bal_prts = prof_register("bal prts", 1., 0, 0);
      pr_bal_flds = prof_register("bal flds", 1., 0, 0);
    }

    prof_start(pr_bal_load);
    auto old_grid = gridp;
    int n_patches_new = find_best_mapping(*old_grid, loads_all);
    prof_stop(pr_bal_load);

    if (n_patches_new < 0) { // unchanged mapping, nothing tbd
      mpi_printf(old_grid->comm(), "***** Balance: decomposition unchanged\n");
      return n_prts_by_patch_old;
    }

    mpi_printf(old_grid->comm(),
               "***** Balance: new decomposition: balancing\n");

    auto new_grid = new Grid_t{old_grid->domain, old_grid->bc, old_grid->kinds,
                               old_grid->norm,   old_grid->dt, n_patches_new};
    new_grid->ibn = old_grid->ibn; // FIXME, sucky ibn handling...
    new_grid->timestep_ = old_grid->timestep_;

    delete[] psc_balance_comp_time_by_patch;
    psc_balance_comp_time_by_patch = new double[new_grid->n_patches()];

    prof_start(pr_bal_ctx);
    communicate_ctx ctx(old_grid->mrc_domain(), new_grid->mrc_domain());
    prof_stop(pr_bal_ctx);

    MEM_STATS();
    // copy particles to host, free on gpu
    Mparticles* p_mp_host = nullptr;
    if (mp && typeid(*mp) != typeid(Mparticles)) {
      auto& mp_base = *mp;
      mpi_printf(old_grid->comm(), "***** Balance: particles to host\n");
      p_mp_host = new Mparticles{mp_base.grid()};
      auto& mp_host = *p_mp_host;
      MparticlesBase::convert(mp_base, mp_host);
#ifdef USE_CUDA
      mp_base.~MparticlesBase();
#else
      mp_base.reset(*new_grid); // frees memory here already
#endif
    }

#ifdef USE_CUDA
    std::vector<std::reference_wrapper<MfieldsCuda>> mfields_cuda;
    mfields_cuda.reserve(MfieldsBase::instances.size());
#endif
    std::vector<std::reference_wrapper<Mfields>> mfields_host;
    mfields_host.reserve(MfieldsBase::instances.size());

    for (auto mf : MfieldsBase::instances) {
#ifdef USE_CUDA
      if (typeid(*mf) == typeid(MfieldsCuda)) {
        auto& mf_cuda = dynamic_cast<MfieldsCuda&>(*mf);
        // mprintf("MfieldsCuda #components %d\n", mf->_n_comps());
        mfields_cuda.emplace_back(dynamic_cast<MfieldsCuda&>(*mf));
        continue;
      }
#endif
      mfields_host.emplace_back(dynamic_cast<Mfields&>(*mf));
    }

#ifdef USE_CUDA
    // copy fields to host, free on gpu
    std::vector<Mfields> mfields_old;
    mfields_old.reserve(mfields_cuda.size());

    for (int n = 0; n < mfields_cuda.size(); n++) {
      MfieldsCuda& mf_cuda = mfields_cuda[n];
      mpi_printf(old_grid->comm(),
                 "***** Balance: copying CUDA field to host, #components %d\n",
                 mf_cuda.n_comps());
      mfields_old.emplace_back(*old_grid, mf_cuda.n_comps(), mf_cuda.ibn());
      auto& mf_old = mfields_old.back();
      MfieldsBase::convert(mf_cuda, mf_old, 0, mf_cuda.n_comps());
      mf_cuda.~MfieldsCuda();
    }
    MEM_STATS();
    BndCuda3::clear();
#endif
    MEM_STATS();

    // particles
    std::vector<uint> n_prts_by_patch_new;
    if (mp) {
      mpi_printf(old_grid->comm(), "***** Balance: balancing particles\n");
      prof_start(pr_bal_prts);

      if (typeid(*mp) == typeid(Mparticles)) {
        ctx.communicate_particles(dynamic_cast<Mparticles&>(*mp), *new_grid);
      } else {
        assert(p_mp_host);
        ctx.communicate_particles(*p_mp_host, *new_grid);
      }

      prof_stop(pr_bal_prts);
    } else {
      n_prts_by_patch_new = ctx.new_n_prts(n_prts_by_patch_old);
    }

    prof_start(pr_bal_flds);
    // state field
    mpi_printf(old_grid->comm(), "***** Balance: balancing state field\n");
    for (auto mf : MfieldsStateBase::instances) {
      // for now MfieldsStateFromMfields will have its underlying Mfields
      // rebalanced, anyway, so all we need ot do is reset the grid.
      // FIXME, that his however broken if using MfieldsState that isn't
      // MfieldsStateFromMfields...
      mf->reset(*new_grid);
    }

    // fields
    for (Mfields& mf : mfields_host) {
      ctx.communicate_fields(mf, *new_grid);
    }

#ifdef USE_CUDA
    //  communicate on host
    for (Mfields& mf : mfields_old) {
      ctx.communicate_fields(mf, *new_grid);
    }

    // move back to gpu
    for (int n = 0; n < mfields_cuda.size(); n++) {
      Mfields& mf = mfields_old[n];
      MfieldsCuda& mf_cuda = mfields_cuda[n];
      new (&mf_cuda) MfieldsCuda{*new_grid, mf.n_comps(), mf.ibn()};
      MfieldsBase::convert(mf, mf_cuda, 0, mf._n_comps());
    }
#endif
    prof_stop(pr_bal_flds);

    // mv particles back to gpu
    if (mp && typeid(*mp) != typeid(Mparticles)) {
#ifdef USE_CUDA
      new (mp) MparticlesCuda<BS444>{*new_grid};
      auto& mp_base = *mp;

      mpi_printf(old_grid->comm(), "***** Balance: particles to device\n");
      assert(p_mp_host);
      auto& mp_host = *p_mp_host;
      MparticlesBase::convert(mp_host, mp_base);
      delete p_mp_host;
#endif
    }

    // update psc etc
    delete gridp;
    gridp = new_grid;
    psc_balance_generation_cnt++;

    return n_prts_by_patch_new;
  }

private:
  double factor_fields_;
  psc::balance::get_loads get_loads_;
  bool print_loads_;
  bool write_loads_;
};
