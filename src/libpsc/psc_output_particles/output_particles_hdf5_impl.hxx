
#include "psc.h"
#include <mrc_profile.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "output_particles.hxx"

// needs to be changed accordingly

struct hdf5_prt {
  float x, y, z;
  float px, py, pz;
  float q, m, w;
};

#define H5_CHK(ierr) assert(ierr >= 0)
#define CE assert(ierr == 0)

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define to_psc_output_particles_hdf5(out) \
  mrc_to_subobj(out, struct psc_output_particles_hdf5)

template<typename Mparticles>
struct psc_output_particles_hdf5 : OutputParticlesParams, OutputParticlesBase
{
  using particle_t = typename Mparticles::particle_t;
  using Particles = typename Mparticles::patch_t;
  using real_t = typename Mparticles::real_t;

  psc_output_particles_hdf5(const OutputParticlesParams& params)
    : OutputParticlesParams(params),
      comm_{psc_comm(ppsc)}
  {
    hid_t id = H5Tcreate(H5T_COMPOUND, sizeof(struct hdf5_prt));
    H5Tinsert(id, "x" , HOFFSET(struct hdf5_prt, x) , H5T_NATIVE_FLOAT);
    H5Tinsert(id, "y" , HOFFSET(struct hdf5_prt, y) , H5T_NATIVE_FLOAT);
    H5Tinsert(id, "z" , HOFFSET(struct hdf5_prt, z) , H5T_NATIVE_FLOAT);
    H5Tinsert(id, "px", HOFFSET(struct hdf5_prt, px), H5T_NATIVE_FLOAT);
    H5Tinsert(id, "py", HOFFSET(struct hdf5_prt, py), H5T_NATIVE_FLOAT);
    H5Tinsert(id, "pz", HOFFSET(struct hdf5_prt, pz), H5T_NATIVE_FLOAT);
    H5Tinsert(id, "q" , HOFFSET(struct hdf5_prt, q) , H5T_NATIVE_FLOAT);
    H5Tinsert(id, "m" , HOFFSET(struct hdf5_prt, m) , H5T_NATIVE_FLOAT);
    H5Tinsert(id, "w" , HOFFSET(struct hdf5_prt, w) , H5T_NATIVE_FLOAT);
    prt_type = id;
    
    // set hi to gdims by default (if not set differently before)
    // and calculate wdims (global dims of region we're writing)
    for (int d = 0; d < 3; d++) {
      assert(lo[d] >= 0);
      if (hi[d] == 0) {
	hi[d] = ppsc->grid().domain.gdims[d];
      }
      assert(hi[d] <= ppsc->grid().domain.gdims[d]);
      wdims[d] = hi[d] - lo[d];
    }
  }

  ~psc_output_particles_hdf5()
  {
    H5Tclose(prt_type);
  }
  
  // ----------------------------------------------------------------------
  // get_cell_index
  // FIXME, lots of stuff here is pretty much duplicated from countsort2

  static inline int
  cell_index_3_to_1(const int *ldims, int j0, int j1, int j2)
  {
    return ((j2) * ldims[1] + j1) * ldims[0] + j0;
  }

  static inline int
  get_sort_index(Particles& prts, particle_t *part)
  {
    const Grid_t& grid = ppsc->grid();
    const int *ldims = grid.ldims;
  
    int j0 = prts.cellPosition(part->x[0], 0);
    int j1 = prts.cellPosition(part->x[1], 1);
    int j2 = prts.cellPosition(part->x[2], 2);
    // FIXME, this is hoping that reason is that we were just on the right bnd...
    if (j0 == ldims[0]) j0--;
    if (j1 == ldims[1]) j1--;
    if (j2 == ldims[2]) j2--;
    assert(j0 >= 0 && j0 < ldims[0]);
    assert(j1 >= 0 && j1 < ldims[1]);
    assert(j2 >= 0 && j2 < ldims[2]);

    int kind = part->kind;
    assert(kind < grid.kinds.size());
 
    return cell_index_3_to_1(ldims, j0, j1, j2) * grid.kinds.size() + kind;
  }

  // ----------------------------------------------------------------------
  // count_sort

  static void count_sort(Mparticles& mprts, int **off, int **map)
  {
    int nr_kinds = mprts.grid().kinds.size();

    for (int p = 0; p < mprts.n_patches(); p++) {
      const int *ldims = ppsc->grid().ldims;
      int nr_indices = ldims[0] * ldims[1] * ldims[2] * nr_kinds;
      off[p] = (int *) calloc(nr_indices + 1, sizeof(*off[p]));
      Particles& prts = mprts[p];
      unsigned int n_prts = prts.size();

      // counting sort to get map 
      for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter) {
	particle_t *part = &*prt_iter;
	int si = get_sort_index(prts, part);
	off[p][si]++;
      }
      // prefix sum to get offsets
      int o = 0;
      int *off2 = (int *) malloc((nr_indices + 1) * sizeof(*off2));
      for (int si = 0; si <= nr_indices; si++) {
	int cnt = off[p][si];
	off[p][si] = o; // this will be saved for later
	off2[si] = o; // this one will overwritten when making the map
	o += cnt;
      }

      // sort a map only, not the actual particles
      map[p] = (int *) malloc(n_prts * sizeof(*map[p]));
      int n = 0;
      for (auto prt_iter = prts.begin(); prt_iter != prts.end(); ++prt_iter, n++) {
	particle_t *part = &*prt_iter;
	int si = get_sort_index(prts, part);
	map[p][off2[si]++] = n;
      }
      free(off2);
    }  
  }

  // ----------------------------------------------------------------------
  // 

  int find_patch_bounds(const int ldims[3], const int off[3],
			int ilo[3], int ihi[3], int ld[3])
  {
    for (int d = 0; d < 3; d++) {
      ilo[d] = MIN(ldims[d], MAX(0, lo[d] - off[d]));
      ihi[d] = MAX(0, MIN(ldims[d], hi[d] - off[d]));
      ilo[d] = MAX(0, lo[d] - off[d]);
      ihi[d] = MIN(ldims[d], hi[d] - off[d]);
      ld[d] = ihi[d] - ilo[d];
    }
    return ppsc->grid().kinds.size() * ld[0] * ld[1] * ld[2];
  }

  // ----------------------------------------------------------------------
  // make_local_particle_array

  hdf5_prt* make_local_particle_array(Mparticles& mprts, int **off, int **map,
				      size_t **idx, size_t *p_n_write, size_t *p_n_off, size_t *p_n_total)
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
	      int si = cell_index_3_to_1(grid.ldims, jx, jy, jz) * nr_kinds + kind;
	      n_write += off[p][si+1] - off[p][si];
	    }
	  }
	}
      }
    }

    assert(sizeof(size_t) == sizeof(unsigned long));
    size_t n_total, n_off = 0;
    MPI_Allreduce(&n_write, &n_total, 1, MPI_LONG, MPI_SUM, comm_);
    MPI_Exscan(&n_write, &n_off, 1, MPI_LONG, MPI_SUM, comm_);

    struct hdf5_prt *arr = (struct hdf5_prt *) malloc(n_write * sizeof(*arr));

    // copy particles to be written into temp array
    int nn = 0;
    for (int p = 0; p < mprts.n_patches(); p++) {
      auto& prts = mprts[p];
      const auto& patch = ppsc->grid().patches[p];
      int ilo[3], ihi[3], ld[3];
      int sz = find_patch_bounds(grid.ldims, patch.off, ilo, ihi, ld);
      idx[p] = (size_t *) malloc(2 * sz * sizeof(*idx));

      for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	  for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	    for (int kind = 0; kind < nr_kinds; kind++) {
	      int si = cell_index_3_to_1(grid.ldims, jx, jy, jz) * nr_kinds + kind;
	      int jj = ((kind * ld[2] + jz - ilo[2])
			* ld[1] + jy - ilo[1]) * ld[0] + jx - ilo[0];
	      idx[p][jj     ] = nn + n_off;
	      idx[p][jj + sz] = nn + n_off + off[p][si+1] - off[p][si];
	      for (int n = off[p][si]; n < off[p][si+1]; n++, nn++) {
		particle_t *part = &prts[map[p][n]];
		arr[nn].x  = part->x[0] + patch.xb[0];
		arr[nn].y  = part->x[1] + patch.xb[1];
		arr[nn].z  = part->x[2] + patch.xb[2];
		arr[nn].px = part->p[0];
		arr[nn].py = part->p[1];
		arr[nn].pz = part->p[2];
		arr[nn].q  = prts.prt_qni(*part);
		arr[nn].m  = prts.prt_mni(*part);
		arr[nn].w  = prts.prt_wni(*part);
	      }
	    }
	  }
	}
      }
    }

    *p_n_write = n_write;
    *p_n_off = n_off;
    *p_n_total = n_total;
    return arr;
  }

  void write_particles(size_t n_write, size_t n_off, size_t n_total,
		       hdf5_prt *arr, hid_t group, hid_t dxpl)
  {
    herr_t ierr;

    hsize_t mdims[1] = { n_write };
    hsize_t fdims[1] = { n_total };
    hsize_t foff[1] = { n_off };
    hid_t memspace = H5Screate_simple(1, mdims, NULL); H5_CHK(memspace);
    hid_t filespace = H5Screate_simple(1, fdims, NULL); H5_CHK(filespace);
    ierr = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, foff, NULL,
			       mdims, NULL); CE;
  
    hid_t dset = H5Dcreate(group, "1d", prt_type, filespace, H5P_DEFAULT,
			   H5P_DEFAULT, H5P_DEFAULT); H5_CHK(dset);
    ierr = H5Dwrite(dset, prt_type, memspace, filespace, dxpl, arr); CE;
  
    ierr = H5Dclose(dset); CE;
    ierr = H5Sclose(filespace); CE;
    ierr = H5Sclose(memspace); CE;
  }

  void write_idx(size_t *gidx_begin, size_t *gidx_end, hid_t group, hid_t dxpl)
  {
    herr_t ierr;

    hsize_t fdims[4];
    fdims[0] = ppsc->grid().kinds.size();
    fdims[1] = wdims[2];
    fdims[2] = wdims[1];
    fdims[3] = wdims[0];
    hid_t filespace = H5Screate_simple(4, fdims, NULL); H5_CHK(filespace);
    hid_t memspace;

    assert(sizeof(size_t) == sizeof(hsize_t));
    int rank;
    MPI_Comm_rank(comm_, &rank);
    if (rank == 0) {
      memspace = H5Screate_simple(4, fdims, NULL); H5_CHK(memspace);
    } else {
      memspace = H5Screate(H5S_NULL);
      H5Sselect_none(memspace);
      H5Sselect_none(filespace);
    }

    hid_t dset = H5Dcreate(group, "idx_begin", H5T_NATIVE_HSIZE, filespace,
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(dset);
    ierr = H5Dwrite(dset, H5T_NATIVE_HSIZE, memspace, filespace, dxpl, gidx_begin); CE;
    ierr = H5Dclose(dset); CE;

    dset = H5Dcreate(group, "idx_end", H5T_NATIVE_HSIZE, filespace,
		     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(dset);
    ierr = H5Dwrite(dset, H5T_NATIVE_HSIZE, memspace, filespace, dxpl, gidx_end); CE;
    ierr = H5Dclose(dset); CE;

    ierr = H5Sclose(filespace); CE;
    ierr = H5Sclose(memspace); CE;
  }

  // ----------------------------------------------------------------------
  // operator()

  void operator()(Mparticles& mprts)
  {
    const auto& grid = ppsc->grid();
    herr_t ierr;

    static int pr_A, pr_B, pr_C, pr_D, pr_E;
    if (!pr_A) {
      pr_A = prof_register("outp: local", 1., 0, 0);
      pr_B = prof_register("outp: comm", 1., 0, 0);
      pr_C = prof_register("outp: write prep", 1., 0, 0);
      pr_D = prof_register("outp: write idx", 1., 0, 0);
      pr_E = prof_register("outp: write prts", 1., 0, 0);
    }

    prof_start(pr_A);
    int rank, size;
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &size);

    struct mrc_patch_info info;
    int nr_kinds = mprts.grid().kinds.size();

    int **off = (int **) malloc(mprts.n_patches() * sizeof(*off));
    int **map = (int **) malloc(mprts.n_patches() * sizeof(*off));

    count_sort(mprts, map, off);

    size_t **idx = (size_t **) malloc(mprts.n_patches() * sizeof(*idx));

    size_t *gidx_begin = NULL, *gidx_end = NULL;
    if (rank == 0) {
      // alloc global idx array
      gidx_begin = (size_t *) malloc(nr_kinds * wdims[0]*wdims[1]*wdims[2] * sizeof(*gidx_begin));
      gidx_end = (size_t *) malloc(nr_kinds * wdims[0]*wdims[1]*wdims[2] * sizeof(*gidx_end));
      for (int jz = 0; jz < wdims[2]; jz++) {
	for (int jy = 0; jy < wdims[1]; jy++) {
	  for (int jx = 0; jx < wdims[0]; jx++) {
	    for (int kind = 0; kind < nr_kinds; kind++) {
	      int ii = ((kind * wdims[2] + jz) * wdims[1] + jy) * wdims[0] + jx;
	      gidx_begin[ii] = -1;
	      gidx_end[ii] = -1;
	    }
	  }
	}
      }
    }

    // find local particle and idx arrays
    size_t n_write, n_off, n_total;
    hdf5_prt *arr = make_local_particle_array(mprts, map, off, idx, &n_write, &n_off, &n_total);
    prof_stop(pr_A);

    prof_start(pr_B);
    if (rank == 0) {
      int nr_global_patches;
      ppsc->mrc_domain_.get_nr_global_patches(&nr_global_patches);

      int *remote_sz = (int *) calloc(size, sizeof(*remote_sz));
      for (int p = 0; p < nr_global_patches; p++) {
	ppsc->mrc_domain_.get_global_patch_info(p, &info);
	if (info.rank == rank) { // skip local patches
	  continue;
	}
	int ilo[3], ihi[3], ld[3];
	int sz = find_patch_bounds(info.ldims, info.off, ilo, ihi, ld);
	remote_sz[info.rank] += sz;
      }

      // build global idx array, local part
      for (int p = 0; p < mprts.n_patches(); p++) {
	ppsc->mrc_domain_.get_local_patch_info(p, &info);
	int ilo[3], ihi[3], ld[3];
	int sz = find_patch_bounds(info.ldims, info.off, ilo, ihi, ld);
      
	for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	  for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	    for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	      for (int kind = 0; kind < nr_kinds; kind++) {
		int ix = jx + info.off[0] - lo[0];
		int iy = jy + info.off[1] - lo[1];
		int iz = jz + info.off[2] - lo[2];
		int jj = ((kind * ld[2] + jz - ilo[2])
			  * ld[1] + jy - ilo[1]) * ld[0] + jx - ilo[0];
		int ii = ((kind * wdims[2] + iz) * wdims[1] + iy) * wdims[0] + ix;
		gidx_begin[ii] = idx[p][jj];
		gidx_end[ii]   = idx[p][jj + sz];
	      }
	    }
	  }
	}
      }

      size_t **recv_buf = (size_t **) calloc(size, sizeof(*recv_buf));
      size_t **recv_buf_p = (size_t **) calloc(size, sizeof(*recv_buf_p));
      MPI_Request *recv_req = (MPI_Request *) calloc(size, sizeof(*recv_req));
      recv_req[rank] = MPI_REQUEST_NULL;
      for (int r = 0; r < size; r++) {
	if (r == rank) { // skip this proc
	  continue;
	}
	recv_buf[r] = (size_t *) malloc(2 * remote_sz[r] * sizeof(*recv_buf[r]));
	recv_buf_p[r] = recv_buf[r];
	MPI_Irecv(recv_buf[r], 2 * remote_sz[r], MPI_LONG, r, 0x4000, comm_,
		  &recv_req[r]);
      }
      MPI_Waitall(size, recv_req, MPI_STATUSES_IGNORE);
      free(recv_req);

      // build global idx array, remote part
      for (int p = 0; p < nr_global_patches; p++) {
	ppsc->mrc_domain_.get_global_patch_info(p, &info);
	if (info.rank == rank) { // skip local patches
	  continue;
	}
	int ilo[3], ihi[3], ld[3];
	int sz = find_patch_bounds(info.ldims, info.off, ilo, ihi, ld);
	for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	  for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	    for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	      for (int kind = 0; kind < nr_kinds; kind++) {
		int ix = jx + info.off[0] - lo[0];
		int iy = jy + info.off[1] - lo[1];
		int iz = jz + info.off[2] - lo[2];
		int jj = ((kind * ld[2] + jz - ilo[2])
			  * ld[1] + jy - ilo[1]) * ld[0] + jx - ilo[0];
		int ii = ((kind * wdims[2] + iz) * wdims[1] + iy) * wdims[0] + ix;
		gidx_begin[ii] = recv_buf_p[info.rank][jj];
		gidx_end[ii]   = recv_buf_p[info.rank][jj + sz];
	      }
	    }
	  }
	}
	recv_buf_p[info.rank] += 2 * sz;
      }

      for (int r = 0; r < size; r++) {
	free(recv_buf[r]);
      }
      free(recv_buf);
      free(recv_buf_p);
      free(remote_sz);
    } else { // rank != 0: send to rank 0
      // FIXME, alloc idx[] as one array in the first place
      int local_sz = 0;
      for (int p = 0; p < mprts.n_patches(); p++) {
	const auto& patch = grid.patches[p];
	int ilo[3], ihi[3], ld[3];
	int sz = find_patch_bounds(grid.ldims, patch.off, ilo, ihi, ld);
	local_sz += sz;
      }
      size_t *l_idx = (size_t *) malloc(2 * local_sz * sizeof(*l_idx));
      size_t *l_idx_p = (size_t *) l_idx;
      for (int p = 0; p < mprts.n_patches(); p++) {
	const auto& patch = grid.patches[p];
	int ilo[3], ihi[3], ld[3];
	int sz = find_patch_bounds(grid.ldims, patch.off, ilo, ihi, ld);
	memcpy(l_idx_p, idx[p], 2 * sz * sizeof(*l_idx_p));
	l_idx_p += 2 * sz;
      }

      MPI_Send(l_idx, 2 * local_sz, MPI_LONG, 0, 0x4000, comm_);

      free(l_idx);
    }
    prof_stop(pr_B);

    prof_start(pr_C);
    char filename[strlen(data_dir) + strlen(basename) + 20];
    sprintf(filename, "%s/%s.%06d_p%06d.h5", data_dir,
	    basename, ppsc->timestep, 0);

    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);

    MPI_Info mpi_info;
    MPI_Info_create(&mpi_info);
#ifdef H5_HAVE_PARALLEL
    if (romio_cb_write) {
      MPI_Info_set(mpi_info, "romio_cb_write", romio_cb_write);
    }
    if (romio_ds_write) {
      MPI_Info_set(mpi_info, "romio_ds_write", romio_ds_write);
    }
    H5Pset_fapl_mpio(plist, comm_, mpi_info);
#endif

    hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist); H5_CHK(file);
    H5Pclose(plist);
    MPI_Info_free(&mpi_info);

    hid_t group = H5Gcreate(file, "particles", H5P_DEFAULT,
			    H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
    hid_t groupp = H5Gcreate(group, "p0", H5P_DEFAULT,
			     H5P_DEFAULT, H5P_DEFAULT); H5_CHK(groupp);

    ierr = H5LTset_attribute_int(group, ".", "lo", lo, 3); CE;
    ierr = H5LTset_attribute_int(group, ".", "hi", hi, 3); CE;

    hid_t dxpl = H5Pcreate(H5P_DATASET_XFER); H5_CHK(dxpl);
#ifdef H5_HAVE_PARALLEL
    if (use_independent_io) {
      ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT); CE;
    } else {
      ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE); CE;
    }
#endif
    prof_stop(pr_C);

    prof_start(pr_D);
    write_idx(gidx_begin, gidx_end, groupp, dxpl);
    prof_stop(pr_D);

    prof_start(pr_E);
    write_particles(n_write, n_off, n_total, arr, groupp, dxpl);

    ierr = H5Pclose(dxpl); CE;

    H5Gclose(groupp);
    H5Gclose(group);
    H5Fclose(file);
    prof_stop(pr_E);

    free(arr);
    free(gidx_begin);
    free(gidx_end);

    for (int p = 0; p < mprts.n_patches(); p++) {
      free(off[p]);
      free(map[p]);
      free(idx[p]);
    }
    free(off);
    free(map);
    free(idx);
  }
  
  // ----------------------------------------------------------------------
  // run

  void run(MparticlesBase& mprts_base) override
  {
    if (every_step <= 0 || ppsc->timestep % every_step != 0) {
      return;
    }

    auto& mprts = mprts_base.get_as<Mparticles>();
    (*this)(mprts);
    mprts_base.put_as(mprts, MP_DONT_COPY);
  }

  // private:
  hid_t prt_type;
  Int3 wdims; // dimensions of the subdomain we're actually writing
  MPI_Comm comm_;
};
