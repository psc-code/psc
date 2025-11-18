
#include "psc.h"
#include <mrc_profile.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "grid.hxx"
#include "output_particles.hxx"

#include "psc_particles_single.h"
#include "../libpsc/vpic/mparticles_vpic.hxx"
#ifdef USE_CUDA
#include "mparticles_cuda.hxx"
#endif

struct hdf5_prt
{
  hdf5_prt() = default;

  template <typename Particle, typename Patch>
  hdf5_prt(const Particle& prt, const Patch& patch)
  {
    x = prt.x()[0] + patch.xb[0];
    y = prt.x()[1] + patch.xb[1];
    z = prt.x()[2] + patch.xb[2];
    px = prt.u()[0];
    py = prt.u()[1];
    pz = prt.u()[2];
    q = prt.q();
    m = prt.m();
    w = prt.w();
    id = prt.id();
    tag = prt.tag();
  }

  float x, y, z;
  float px, py, pz;
  float q, m, w;
  psc::particle::Id id;
  psc::particle::Tag tag;
};

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

namespace hdf5
{

// ----------------------------------------------------------------------
// ToHdf5Type helper class

template <typename T>
struct H5Type;

template <>
struct H5Type<int>
{
  static hid_t value() { return H5T_NATIVE_INT; }
};

template <>
struct H5Type<unsigned long>
{
  static hid_t value() { return H5T_NATIVE_ULONG; }
};

template <>
struct H5Type<unsigned long long>
{
  static hid_t value() { return H5T_NATIVE_ULLONG; }
};

} // namespace hdf5

// ======================================================================
// Hdf5ParticleType

class Hdf5ParticleType
{
public:
  Hdf5ParticleType()
  {
    id_ = H5Tcreate(H5T_COMPOUND, sizeof(struct hdf5_prt));
    H5Tinsert(id_, "x", HOFFSET(struct hdf5_prt, x), H5T_NATIVE_FLOAT);
    H5Tinsert(id_, "y", HOFFSET(struct hdf5_prt, y), H5T_NATIVE_FLOAT);
    H5Tinsert(id_, "z", HOFFSET(struct hdf5_prt, z), H5T_NATIVE_FLOAT);
    H5Tinsert(id_, "px", HOFFSET(struct hdf5_prt, px), H5T_NATIVE_FLOAT);
    H5Tinsert(id_, "py", HOFFSET(struct hdf5_prt, py), H5T_NATIVE_FLOAT);
    H5Tinsert(id_, "pz", HOFFSET(struct hdf5_prt, pz), H5T_NATIVE_FLOAT);
    H5Tinsert(id_, "q", HOFFSET(struct hdf5_prt, q), H5T_NATIVE_FLOAT);
    H5Tinsert(id_, "m", HOFFSET(struct hdf5_prt, m), H5T_NATIVE_FLOAT);
    H5Tinsert(id_, "w", HOFFSET(struct hdf5_prt, w), H5T_NATIVE_FLOAT);
    H5Tinsert(id_, "id", HOFFSET(struct hdf5_prt, id),
              hdf5::H5Type<psc::particle::Id>::value());
    H5Tinsert(id_, "tag", HOFFSET(struct hdf5_prt, tag),
              hdf5::H5Type<psc::particle::Tag>::value());
  }

  ~Hdf5ParticleType() { H5Tclose(id_); }

  operator hid_t() const { return id_; }

private:
  hid_t id_;
};

class OutputParticlesWriterHDF5
{
public:
  using particle_type = hdf5_prt;
  using particles_type = std::vector<particle_type>;

  explicit OutputParticlesWriterHDF5(const OutputParticlesParams& params,
                                     const Grid_t::Kinds& kinds, MPI_Comm comm)
    : params_{params}, wdims_{params.hi - params.lo}, kinds_{kinds}, comm_{comm}
  {}

  void operator()(const gt::gtensor<size_t, 4>& gidx_begin,
                  const gt::gtensor<size_t, 4>& gidx_end, size_t n_off,
                  size_t n_total, const particles_type& arr, const Grid_t& grid)
  {
    int ierr;

    static int pr_C, pr_D, pr_E;
    if (!pr_C) {
      pr_C = prof_register("outp: write prep", 1., 0, 0);
      pr_D = prof_register("outp: write idx", 1., 0, 0);
      pr_E = prof_register("outp: write prts", 1., 0, 0);
    }

    prof_start(pr_C);

    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);

    MPI_Info mpi_info;
    MPI_Info_create(&mpi_info);

#ifdef H5_HAVE_PARALLEL
    if (params_.romio_cb_write) {
      MPI_Info_set(mpi_info, (char*)"romio_cb_write",
                   (char*)params_.romio_cb_write);
    }
    if (params_.romio_ds_write) {
      MPI_Info_set(mpi_info, (char*)"romio_ds_write",
                   (char*)params_.romio_ds_write);
    }
    H5Pset_fapl_mpio(plist, comm_, mpi_info);
#else
    mprintf("ERROR: particle output requires parallel hdf5\n");
    abort();
#endif

    // FIXME there must be a safer way to do this.

    // Currently, if timestep has too many digits, the filename string gets
    // truncated to e.g. prt.1234567890.h (with .h instead of .h5)

    // This is probably a problem with other outputs, too
    int slen = strlen(params_.data_dir) + strlen(params_.basename) + 15;
    char filename[slen];
    if (grid.timestep() > 999999999) {
      LOG_WARN("%s step index exceeds digit limit in file name. Data still "
               "written, but file name is truncated.\n",
               params_.basename);
    }
    snprintf(filename, slen, "%s/%s.%09d.h5", params_.data_dir,
             params_.basename, grid.timestep());

    hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    H5_CHK(file);
    H5Pclose(plist);
    MPI_Info_free(&mpi_info);

    hid_t group =
      H5Gcreate(file, "particles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5_CHK(group);
    hid_t groupp =
      H5Gcreate(group, "p0", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5_CHK(groupp);

    ierr = H5LTset_attribute_int(group, ".", "lo", params_.lo, 3);
    CE;
    ierr = H5LTset_attribute_int(group, ".", "hi", params_.hi, 3);
    CE;

    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
    H5_CHK(dxpl);
#ifdef H5_HAVE_PARALLEL
    if (params_.use_independent_io) {
      ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT);
      CE;
    } else {
      ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
      CE;
    }
#endif
    prof_stop(pr_C);

    write_time(grid.time(), file, dxpl);

    prof_start(pr_D);
    write_idx(gidx_begin, gidx_end, group, dxpl);
    prof_stop(pr_D);

    prof_start(pr_E);
    write_particles(n_off, n_total, arr, groupp, dxpl);
    prof_stop(pr_E);

    ierr = H5Pclose(dxpl);
    CE;
    ierr = H5Gclose(groupp);
    CE;
    ierr = H5Gclose(group);
    CE;
    ierr = H5Fclose(file);
    CE;
  }

private:
  void write_time(double time, hid_t group, hid_t dxpl)
  {
    herr_t ierr;

    hid_t filespace = H5Screate(H5S_SCALAR);
    H5_CHK(filespace);
    hid_t memspace;

    int mpi_rank;
    MPI_Comm_rank(comm_, &mpi_rank);
    if (mpi_rank == 0) {
      memspace = H5Screate(H5S_SCALAR);
      H5_CHK(memspace);
    } else {
      memspace = H5Screate(H5S_NULL);
      H5Sselect_none(memspace);
      H5Sselect_none(filespace);
    }

    hid_t dset = H5Dcreate(group, "time", H5T_NATIVE_DOUBLE, filespace,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5_CHK(dset);
    ierr = H5Dwrite(dset, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl, &time);
    CE;
    ierr = H5Dclose(dset);
    CE;

    ierr = H5Sclose(filespace);
    CE;
    ierr = H5Sclose(memspace);
    CE;
  }

  void write_idx(const gt::gtensor<size_t, 4>& gidx_begin,
                 const gt::gtensor<size_t, 4>& gidx_end, hid_t group,
                 hid_t dxpl)
  {
    herr_t ierr;

    hsize_t fdims[4];
    fdims[0] = kinds_.size();
    fdims[1] = wdims_[2];
    fdims[2] = wdims_[1];
    fdims[3] = wdims_[0];
    hid_t filespace = H5Screate_simple(4, fdims, NULL);
    H5_CHK(filespace);
    hid_t memspace;

    assert(sizeof(size_t) == sizeof(hsize_t));
    int mpi_rank;
    MPI_Comm_rank(comm_, &mpi_rank);
    if (mpi_rank == 0) {
      memspace = H5Screate_simple(4, fdims, NULL);
      H5_CHK(memspace);
    } else {
      memspace = H5Screate(H5S_NULL);
      H5Sselect_none(memspace);
      H5Sselect_none(filespace);
    }

    hid_t dset = H5Dcreate(group, "idx_begin", H5T_NATIVE_HSIZE, filespace,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5_CHK(dset);
    ierr = H5Dwrite(dset, H5T_NATIVE_HSIZE, memspace, filespace, dxpl,
                    gidx_begin.data());
    CE;
    ierr = H5Dclose(dset);
    CE;

    dset = H5Dcreate(group, "idx_end", H5T_NATIVE_HSIZE, filespace, H5P_DEFAULT,
                     H5P_DEFAULT, H5P_DEFAULT);
    H5_CHK(dset);
    ierr = H5Dwrite(dset, H5T_NATIVE_HSIZE, memspace, filespace, dxpl,
                    gidx_end.data());
    CE;
    ierr = H5Dclose(dset);
    CE;

    ierr = H5Sclose(filespace);
    CE;
    ierr = H5Sclose(memspace);
    CE;
  }

  void write_particles(size_t n_off, size_t n_total, const particles_type& arr,
                       hid_t group, hid_t dxpl)
  {
    herr_t ierr;

    hsize_t mdims[1] = {arr.size()};
    hsize_t fdims[1] = {n_total};
    hsize_t foff[1] = {n_off};
    hid_t memspace = H5Screate_simple(1, mdims, NULL);
    H5_CHK(memspace);
    hid_t filespace = H5Screate_simple(1, fdims, NULL);
    H5_CHK(filespace);
    ierr =
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, foff, NULL, mdims, NULL);
    CE;

    hid_t dset = H5Dcreate(group, "1d", prt_type_, filespace, H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);
    H5_CHK(dset);
    ierr = H5Dwrite(dset, prt_type_, memspace, filespace, dxpl, arr.data());
    CE;

    ierr = H5Dclose(dset);
    CE;
    ierr = H5Sclose(filespace);
    CE;
    ierr = H5Sclose(memspace);
    CE;
  }

private:
  OutputParticlesParams params_;
  Int3 wdims_;
  Grid_t::Kinds kinds_;
  MPI_Comm comm_;
  Hdf5ParticleType prt_type_;
};

template <int EVERY>
class ParticleSelectorEveryNth
{
public:
  template <typename Particle>
  bool operator()(const Particle& prt)
  {
    // return true for every `every_`th particle
    return (cnt_++ % EVERY) == 0;
  }

private:
  int cnt_ = 0;
};

class ParticleSelectorAll
{
public:
  template <typename Particle>
  bool operator()(const Particle& prt)
  {
    return true;
  }
};

namespace detail
{

// ======================================================================
// OutputParticlesHdf5

template <typename Mparticles, typename ParticleSelector>
struct OutputParticlesHdf5
{
  using writer_type = OutputParticlesWriterHDF5;
  using writer_particles_type = writer_type::particles_type;
  using Particles = typename Mparticles::Patch;

  OutputParticlesHdf5(const Grid_t& grid, const OutputParticlesParams& params)
    : lo_{params.lo},
      hi_{params.hi},
      wdims_{params.hi - params.lo},
      kinds_{grid.kinds},
      comm_{grid.comm()}
  {}

  // ----------------------------------------------------------------------
  // get_sort_index
  // FIXME, lots of stuff here is pretty much duplicated from countsort2

  static inline int sort_index(const int* ldims, int nr_kinds, Int3 idx3,
                               int kind)
  {
    return ((idx3[2] * ldims[1] + idx3[1]) * ldims[0] + idx3[0]) * nr_kinds +
           kind;
  };

  static inline int get_sort_index(Particles& prts, int n)
  {
    const Grid_t& grid = prts.grid();
    const int* ldims = grid.ldims;
    const auto& prt = prts[n];

    Int3 pos;
    for (int d = 0; d < 3; d++) {
      pos[d] = prts.cellPosition(prt.x[d], d);
      // FIXME, this is hoping that reason is that we were just on the right
      // bnd...
      if (pos[d] == ldims[d])
        pos[d]--;
      assert(pos[d] >= 0 && pos[d] < ldims[d]);
    }
    assert(prt.kind < grid.kinds.size());

    return sort_index(ldims, grid.kinds.size(), pos, prt.kind);
  }

  // ----------------------------------------------------------------------
  // count_sort

  static void count_sort(Mparticles& mprts, std::vector<std::vector<int>>& off,
                         std::vector<std::vector<int>>& map)
  {
    int nr_kinds = mprts.grid().kinds.size();

    ParticleSelector selector;

    for (int p = 0; p < mprts.n_patches(); p++) {
      const int* ldims = mprts.grid().ldims;
      int nr_indices = ldims[0] * ldims[1] * ldims[2] * nr_kinds;
      off[p].resize(nr_indices + 1);
      auto&& prts = mprts[p];
      unsigned int n_prts = prts.size();
      std::vector<int> particle_indices;
      particle_indices.reserve(n_prts);

      // counting sort to get map
      for (int n = 0; n < n_prts; n++) {
        if (selector(prts[n])) {
          particle_indices.push_back(n);
          int si = get_sort_index(prts, n);
          off[p][si]++;
        }
      }
      // prefix sum to get offsets
      int o = 0;
      std::vector<int> off2(nr_indices + 1);
      for (int si = 0; si <= nr_indices; si++) {
        int cnt = off[p][si];
        off[p][si] = o; // this will be saved for later
        off2[si] = o;   // this one will overwritten when making the map
        o += cnt;
      }

      // sort a map only, not the actual particles
      map[p].resize(n_prts);
      for (auto n : particle_indices) {
        int si = get_sort_index(prts, n);
        map[p][off2[si]++] = n;
      }
    }
  }

  // ----------------------------------------------------------------------
  //

  int find_patch_bounds(const int ldims[3], const int off[3], int ilo[3],
                        int ihi[3], int ld[3])
  {
    for (int d = 0; d < 3; d++) {
      ilo[d] = std::max(lo_[d], off[d]) - off[d];
      ihi[d] = std::min(hi_[d], off[d] + ldims[d]) - off[d];
      ilo[d] = std::min(ldims[d], ilo[d]);
      ihi[d] = std::max(0, ihi[d]);
      ld[d] = ihi[d] - ilo[d];
    }
    return kinds_.size() * ld[0] * ld[1] * ld[2];
  }

  // ----------------------------------------------------------------------
  // make_local_particle_array

  writer_particles_type make_local_particle_array(
    Mparticles& mprts, const std::vector<std::vector<int>>& off,
    const std::vector<std::vector<int>>& map,
    std::vector<gt::gtensor<size_t, 5>>& idx, size_t* p_n_off,
    size_t* p_n_total)
  {
    const auto& grid = mprts.grid();
    int nr_kinds = grid.kinds.size();

    // count all particles to be written locally
    size_t n_write = 0;
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto& patch = grid.patches[p];
      int ilo[3], ihi[3], ld[3];
      int sz = find_patch_bounds(grid.ldims, patch.off, ilo, ihi, ld);
      for (int jz = ilo[2]; jz < ihi[2]; jz++) {
        for (int jy = ilo[1]; jy < ihi[1]; jy++) {
          for (int jx = ilo[0]; jx < ihi[0]; jx++) {
            for (int kind = 0; kind < nr_kinds; kind++) {
              int si = sort_index(grid.ldims, nr_kinds, {jx, jy, jz}, kind);
              n_write += off[p][si + 1] - off[p][si];
            }
          }
        }
      }
    }

    assert(sizeof(size_t) == sizeof(unsigned long));
    size_t n_total, n_off = 0;
    MPI_Allreduce(&n_write, &n_total, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm_);
    MPI_Exscan(&n_write, &n_off, 1, MPI_UNSIGNED_LONG, MPI_SUM, comm_);

    writer_particles_type arr(n_write);

    // copy particles to be written into temp array
    int nn = 0;
    {
      auto accessor = mprts.accessor();
      for (int p = 0; p < mprts.n_patches(); p++) {
        auto prts = accessor[p];
        const auto& patch = grid.patches[p];
        int ilo[3], ihi[3], ld[3];
        int sz = find_patch_bounds(grid.ldims, patch.off, ilo, ihi, ld);
        idx[p] = gt::empty<size_t>({ld[0], ld[1], ld[2], nr_kinds, 2});

        for (int jz = ilo[2]; jz < ihi[2]; jz++) {
          for (int jy = ilo[1]; jy < ihi[1]; jy++) {
            for (int jx = ilo[0]; jx < ihi[0]; jx++) {
              for (int kind = 0; kind < nr_kinds; kind++) {
                int si = sort_index(grid.ldims, nr_kinds, {jx, jy, jz}, kind);
                idx[p](jx - ilo[0], jy - ilo[1], jz - ilo[2], kind, 0) =
                  nn + n_off;
                idx[p](jx - ilo[0], jy - ilo[1], jz - ilo[2], kind, 1) =
                  nn + n_off + off[p][si + 1] - off[p][si];
                for (int n = off[p][si]; n < off[p][si + 1]; n++, nn++) {
                  arr[nn] = writer_type::particle_type(prts[map[p][n]], patch);
                }
              }
            }
          }
        }
      }
    }
    *p_n_off = n_off;
    *p_n_total = n_total;
    return arr;
  }

  // ----------------------------------------------------------------------
  // operator()

  void operator()(Mparticles& mprts, OutputParticlesWriterHDF5& writer)
  {
    const Grid_t& grid = mprts.grid();
    herr_t ierr;

    static int pr_A, pr_B, pr_C;
    if (!pr_A) {
      pr_A = prof_register("outp: local", 1., 0, 0);
      pr_B = prof_register("outp: comm", 1., 0, 0);
      pr_C = prof_register("outp: write prep", 1., 0, 0);
    }

    prof_start(pr_A);
    int mpi_rank, mpi_size;
    MPI_Comm_rank(comm_, &mpi_rank);
    MPI_Comm_size(comm_, &mpi_size);

    struct mrc_patch_info info;
    int nr_kinds = mprts.grid().kinds.size();

    std::vector<std::vector<int>> off(mprts.n_patches());
    std::vector<std::vector<int>> map(mprts.n_patches());
    std::vector<gt::gtensor<size_t, 5>> idx(mprts.n_patches());

    count_sort(mprts, off, map);

    gt::gtensor<size_t, 4> gidx_begin, gidx_end;
    if (mpi_rank == 0) {
      // alloc global idx array
      gidx_begin =
        gt::empty<size_t>({wdims_[0], wdims_[1], wdims_[2], nr_kinds});
      gidx_end = gt::empty<size_t>({wdims_[0], wdims_[1], wdims_[2], nr_kinds});
    }

    // find local particle and idx arrays
    size_t n_off, n_total;
    auto arr =
      make_local_particle_array(mprts, off, map, idx, &n_off, &n_total);
    prof_stop(pr_A);

    prof_start(pr_B);
    if (mpi_rank == 0) {
      int nr_global_patches = grid.nGlobalPatches();

      std::vector<int> remote_sz(mpi_size);
      for (int p = 0; p < nr_global_patches; p++) {
        info = grid.globalPatchInfo(p);
        if (info.rank == mpi_rank) { // skip local patches
          continue;
        }
        int ilo[3], ihi[3], ld[3];
        int sz = find_patch_bounds(info.ldims, info.off, ilo, ihi, ld);
        remote_sz[info.rank] += sz;
      }

      // build global idx array, local part
      for (int p = 0; p < mprts.n_patches(); p++) {
        info = grid.localPatchInfo(p);
        int ilo[3], ihi[3], ld[3];
        int sz = find_patch_bounds(info.ldims, info.off, ilo, ihi, ld);

        for (int jz = ilo[2]; jz < ihi[2]; jz++) {
          for (int jy = ilo[1]; jy < ihi[1]; jy++) {
            for (int jx = ilo[0]; jx < ihi[0]; jx++) {
              for (int kind = 0; kind < nr_kinds; kind++) {
                int ix = jx + info.off[0] - lo_[0];
                int iy = jy + info.off[1] - lo_[1];
                int iz = jz + info.off[2] - lo_[2];
                gidx_begin(ix, iy, iz, kind) = idx[p](jx, jy, jz, kind, 0);
                gidx_end(ix, iy, iz, kind) = idx[p](jx, jy, jz, kind, 1);
              }
            }
          }
        }
      }

      std::vector<std::vector<size_t>> recv_buf(mpi_size);
      std::vector<size_t*> recv_buf_p(mpi_size);
      std::vector<MPI_Request> recv_req(mpi_size, MPI_REQUEST_NULL);
      for (int r = 0; r < mpi_size; r++) {
        if (r == mpi_rank) { // skip this proc
          continue;
        }
        recv_buf[r].resize(2 * remote_sz[r]);
        recv_buf_p[r] = recv_buf[r].data();
        MPI_Irecv(recv_buf[r].data(), recv_buf[r].size(), MPI_UNSIGNED_LONG, r,
                  0x4000, comm_, &recv_req[r]);
      }
      MPI_Waitall(mpi_size, recv_req.data(), MPI_STATUSES_IGNORE);

      // build global idx array, remote part
      for (int p = 0; p < nr_global_patches; p++) {
        info = grid.globalPatchInfo(p);
        if (info.rank == mpi_rank) { // skip local patches
          continue;
        }
        int ilo[3], ihi[3], ld[3];
        int sz = find_patch_bounds(info.ldims, info.off, ilo, ihi, ld);
        for (int jz = ilo[2]; jz < ihi[2]; jz++) {
          for (int jy = ilo[1]; jy < ihi[1]; jy++) {
            for (int jx = ilo[0]; jx < ihi[0]; jx++) {
              for (int kind = 0; kind < nr_kinds; kind++) {
                int ix = jx + info.off[0] - lo_[0];
                int iy = jy + info.off[1] - lo_[1];
                int iz = jz + info.off[2] - lo_[2];
                int jj =
                  ((kind * ld[2] + jz - ilo[2]) * ld[1] + jy - ilo[1]) * ld[0] +
                  jx - ilo[0];
                gidx_begin(ix, iy, iz, kind) = recv_buf_p[info.rank][jj];
                gidx_end(ix, iy, iz, kind) = recv_buf_p[info.rank][jj + sz];
              }
            }
          }
        }
        recv_buf_p[info.rank] += 2 * sz;
      }
    } else { // mpi_rank != 0: send to mpi_rank 0
      // FIXME, alloc idx[] as one array in the first place
      int local_sz = 0;
      for (int p = 0; p < mprts.n_patches(); p++) {
        const auto& patch = grid.patches[p];
        int ilo[3], ihi[3], ld[3];
        int sz = find_patch_bounds(grid.ldims, patch.off, ilo, ihi, ld);
        local_sz += sz;
      }
      std::vector<size_t> l_idx(2 * local_sz);
      size_t* l_idx_p = l_idx.data();
      for (int p = 0; p < mprts.n_patches(); p++) {
        const auto& patch = grid.patches[p];
        int ilo[3], ihi[3], ld[3];
        int sz = find_patch_bounds(grid.ldims, patch.off, ilo, ihi, ld);
        memcpy(l_idx_p, idx[p].data(), 2 * sz * sizeof(*l_idx_p));
        l_idx_p += 2 * sz;
      }

      MPI_Send(l_idx.data(), l_idx.size(), MPI_UNSIGNED_LONG, 0, 0x4000, comm_);
    }
    prof_stop(pr_B);

    writer(gidx_begin, gidx_end, n_off, n_total, arr, grid);
  }

private:
  Int3 lo_, hi_, wdims_;
  Grid_t::Kinds kinds_;
  MPI_Comm comm_;
};

} // namespace detail

template <typename ParticleSelector = ParticleSelectorAll>
class OutputParticlesHdf5 : OutputParticlesBase
{
  static OutputParticlesParams adjust_params(
    const OutputParticlesParams& params_in, const Grid_t& grid)
  {
    OutputParticlesParams params = params_in;
    for (int d = 0; d < 3; d++) {
      if (params.hi[d] == 0) {
        params.hi[d] = grid.domain.gdims[d];
      }
      assert(params.lo[d] >= 0);
      assert(params.hi[d] <= grid.domain.gdims[d]);
    }
    return params;
  }

public:
  OutputParticlesHdf5(const Grid_t& grid, const OutputParticlesParams& params)
    : params_{adjust_params(params, grid)},
      writer_{params_, grid.kinds, grid.comm()}
  {}

  template <typename Mparticles>
  void operator()(Mparticles& mprts)
  {
    const auto& grid = mprts.grid();

    if (params_.every_step <= 0 || grid.timestep() % params_.every_step != 0) {
      return;
    }

    detail::OutputParticlesHdf5<Mparticles, ParticleSelector> impl{grid,
                                                                   params_};
    impl(mprts, writer_);
  }

  // FIXME, handles MparticlesVpic by conversion for now
  template <typename Particles>
  void operator()(MparticlesVpic_<Particles>& mprts)
  {
    const auto& grid = mprts.grid();

    if (params_.every_step <= 0 || grid.timestep() % params_.every_step != 0) {
      return;
    }

    auto& mprts_single = mprts.template get_as<MparticlesSingle>();
    (*this)(mprts_single);
    mprts.put_as(mprts_single, MP_DONT_COPY);
  }

#ifdef USE_CUDA
  // FIXME, handles MparticlesCuda by conversion for now
  template <typename BS>
  void operator()(MparticlesCuda<BS>& mprts)
  {
    const auto& grid = mprts.grid();

    if (params_.every_step <= 0 || grid.timestep() % params_.every_step != 0) {
      return;
    }

    auto& mprts_single = mprts.template get_as<MparticlesSingle>();
    (*this)(mprts_single);
    mprts.put_as(mprts_single, MP_DONT_COPY);
  }
#endif

private:
  const OutputParticlesParams params_;
  OutputParticlesWriterHDF5 writer_;
};
