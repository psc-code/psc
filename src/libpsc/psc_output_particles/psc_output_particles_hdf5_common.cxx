
// when changing the following struct, the HDF5 compound data type prt_type
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

struct psc_output_particles_hdf5 {
  // parameters
  const char *data_dir;
  const char *basename;
  int every_step;
  int lo[3];
  int hi[3];
  bool use_independent_io;
  char *romio_cb_write;
  char *romio_ds_write;

  // internal
  hid_t prt_type;
  int wdims[3]; // dimensions of the subdomain we're actually writing
};

#define VAR(x) (void *)offsetof(struct psc_output_particles_hdf5, x)
static struct param psc_output_particles_hdf5_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(".")       },
  { "basename"           , VAR(basename)             , PARAM_STRING("prt")     },
  { "every_step"         , VAR(every_step)           , PARAM_INT(-1)           },
  { "lo"                 , VAR(lo)                   , PARAM_INT3(0, 0, 0)     },
  { "hi"                 , VAR(hi)                   , PARAM_INT3(0, 0, 0)     },
  { "use_independent_io" , VAR(use_independent_io)   , PARAM_BOOL(false)       },
  { "romio_cb_write"     , VAR(romio_cb_write)       , PARAM_STRING(NULL)      },
  { "romio_ds_write"     , VAR(romio_ds_write)       , PARAM_STRING(NULL)      },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_output_particles_hdf5_create

static void
psc_output_particles_hdf5_create(struct psc_output_particles *out)
{
  struct psc_output_particles_hdf5 *hdf5 = to_psc_output_particles_hdf5(out);

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
  
  hdf5->prt_type = id;
}

// ----------------------------------------------------------------------
// psc_output_particles_hdf5_setup

static void
psc_output_particles_hdf5_setup(struct psc_output_particles *out)
{
  struct psc_output_particles_hdf5 *hdf5 = to_psc_output_particles_hdf5(out);

  // set hi to gdims by default (if not set differently before)
  // and calculate wdims (global dims of region we're writing)
  for (int d = 0; d < 3; d++) {
    assert(hdf5->lo[d] >= 0);
    if (hdf5->hi[d] == 0) {
      hdf5->hi[d] = ppsc->domain.gdims[d];
    }
    assert(hdf5->hi[d] <= ppsc->domain.gdims[d]);
    hdf5->wdims[d] = hdf5->hi[d] - hdf5->lo[d];
  }
}

// ----------------------------------------------------------------------
// psc_output_particles_hdf5_destroy

static void
psc_output_particles_hdf5_destroy(struct psc_output_particles *out)
{
  struct psc_output_particles_hdf5 *hdf5 = to_psc_output_particles_hdf5(out);

  H5Tclose(hdf5->prt_type);
}

// ----------------------------------------------------------------------
// get_cell_index
// FIXME, lots of stuff here is pretty much duplicated from countsort2

static inline int
cell_index_3_to_1(int *ldims, int j0, int j1, int j2)
{
  return ((j2) * ldims[1] + j1) * ldims[0] + j0;
}

static inline int
get_sort_index(int p, particle_t *part)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_real_t dxi = 1.f / patch->dx[0];
  particle_real_t dyi = 1.f / patch->dx[1];
  particle_real_t dzi = 1.f / patch->dx[2];
  int *ldims = patch->ldims;
  
  particle_real_t u = part->xi * dxi;
  particle_real_t v = part->yi * dyi;
  particle_real_t w = part->zi * dzi;
  int j0 = fint(u);
  int j1 = fint(v);
  int j2 = fint(w);
  if (u == ldims[0]) j0--;
  if (v == ldims[1]) j1--;
  if (w == ldims[2]) j2--;
  assert(j0 >= 0 && j0 < ldims[0]);
  assert(j1 >= 0 && j1 < ldims[1]);
  if (!(j2 >= 0 && j2 < ldims[2])) {
    mprintf("j2 %d ldims %d w %g\n", j2, ldims[2], w);
  }
  assert(j2 >= 0 && j2 < ldims[2]);

  int kind = part->kind();
  assert(kind < ppsc->nr_kinds);
 
  return cell_index_3_to_1(ldims, j0, j1, j2) * ppsc->nr_kinds + kind;
}

// ----------------------------------------------------------------------
// count_sort

static void
count_sort(mparticles_t mprts, int **off, int **map)
{
  int nr_kinds = ppsc->nr_kinds;

  for (int p = 0; p < mprts.n_patches(); p++) {
    int *ldims = ppsc->patch[p].ldims;
    int nr_indices = ldims[0] * ldims[1] * ldims[2] * nr_kinds;
    off[p] = (int *) calloc(nr_indices + 1, sizeof(*off[p]));
    particle_range_t prts = mprts[p].range();
    unsigned int n_prts = prts.size();

    // counting sort to get map 
    PARTICLE_ITER_LOOP(prt_iter, prts.begin, prts.end) {
      particle_t *part = &*prt_iter;
      int si = get_sort_index(p, part);
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
    for (particle_iter_t prt_iter = prts.begin;
	 !particle_iter_equal(prt_iter, prts.end);
	 prt_iter = particle_iter_next(prt_iter), n++) {
      particle_t *part = &*prt_iter;
      int si = get_sort_index(p, part);
      map[p][off2[si]++] = n;
    }
    free(off2);
  }  
}

// ----------------------------------------------------------------------
// 

static void
find_patch_bounds(struct psc_output_particles_hdf5 *hdf5,
		  struct mrc_patch_info *info,
		  int ilo[3], int ihi[3], int ld[3], int *p_sz)
{
  int *ldims = info->ldims, *off = info->off;

  for (int d = 0; d < 3; d++) {
    ilo[d] = MIN(ldims[d], MAX(0, hdf5->lo[d] - off[d]));
    ihi[d] = MAX(0, MIN(ldims[d], hdf5->hi[d] - off[d]));
    ilo[d] = MAX(0, hdf5->lo[d] - off[d]);
    ihi[d] = MIN(ldims[d], hdf5->hi[d] - off[d]);
    ld[d] = ihi[d] - ilo[d];
  }
  *p_sz = ppsc->nr_kinds * ld[0] * ld[1] * ld[2];
}

// ----------------------------------------------------------------------
// make_local_particle_array

static struct hdf5_prt *
make_local_particle_array(struct psc_output_particles *out,
			  mparticles_t mprts, int **off, int **map,
			  size_t **idx, size_t *p_n_write, size_t *p_n_off, size_t *p_n_total)
{
  struct psc_output_particles_hdf5 *hdf5 = to_psc_output_particles_hdf5(out);
  MPI_Comm comm = psc_output_particles_comm(out);
  int nr_kinds = ppsc->nr_kinds;
  struct mrc_patch_info info;

  // count all particles to be written locally
  size_t n_write = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    mrc_domain_get_local_patch_info(ppsc->mrc_domain, p, &info);
    int ilo[3], ihi[3], ld[3], sz;
    find_patch_bounds(hdf5, &info, ilo, ihi, ld, &sz);
    for (int jz = ilo[2]; jz < ihi[2]; jz++) {
      for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	  for (int kind = 0; kind < nr_kinds; kind++) {
	    int si = cell_index_3_to_1(info.ldims, jx, jy, jz) * nr_kinds + kind;
	    n_write += off[p][si+1] - off[p][si];
	  }
	}
      }
    }
  }

  assert(sizeof(size_t) == sizeof(unsigned long));
  size_t n_total, n_off = 0;
  MPI_Allreduce(&n_write, &n_total, 1, MPI_LONG, MPI_SUM, comm);
  MPI_Exscan(&n_write, &n_off, 1, MPI_LONG, MPI_SUM, comm);

  struct hdf5_prt *arr = (struct hdf5_prt *) malloc(n_write * sizeof(*arr));

  // copy particles to be written into temp array
  int nn = 0;
  for (int p = 0; p < mprts.n_patches(); p++) {
    particle_range_t prts = mprts[p].range();
    mrc_domain_get_local_patch_info(ppsc->mrc_domain, p, &info);
    int ilo[3], ihi[3], ld[3], sz;
    find_patch_bounds(hdf5, &info, ilo, ihi, ld, &sz);
    idx[p] = (size_t *) malloc(2 * sz * sizeof(*idx));
    struct psc_patch *patch = &ppsc->patch[p];

    for (int jz = ilo[2]; jz < ihi[2]; jz++) {
      for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	  for (int kind = 0; kind < nr_kinds; kind++) {
	    int si = cell_index_3_to_1(info.ldims, jx, jy, jz) * nr_kinds + kind;
	    int jj = ((kind * ld[2] + jz - ilo[2])
		      * ld[1] + jy - ilo[1]) * ld[0] + jx - ilo[0];
	    idx[p][jj     ] = nn + n_off;
	    idx[p][jj + sz] = nn + n_off + off[p][si+1] - off[p][si];
	    for (int n = off[p][si]; n < off[p][si+1]; n++, nn++) {
	      particle_t *part = &prts.begin[map[p][n]];
	      arr[nn].x  = part->xi + patch->xb[0];
	      arr[nn].y  = part->yi + patch->xb[1];
	      arr[nn].z  = part->zi + patch->xb[2];
	      arr[nn].px = part->pxi;
	      arr[nn].py = part->pyi;
	      arr[nn].pz = part->pzi;
	      arr[nn].q  = particle_qni(part);
	      arr[nn].m  = particle_mni(part);
	      arr[nn].w  = particle_wni(part);
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

static void
write_particles(struct psc_output_particles *out, size_t n_write, size_t n_off,
		size_t n_total, struct hdf5_prt *arr, hid_t group, hid_t dxpl)
{
  struct psc_output_particles_hdf5 *hdf5 = to_psc_output_particles_hdf5(out);
  herr_t ierr;

  hsize_t mdims[1] = { n_write };
  hsize_t fdims[1] = { n_total };
  hsize_t foff[1] = { n_off };
  hid_t memspace = H5Screate_simple(1, mdims, NULL); H5_CHK(memspace);
  hid_t filespace = H5Screate_simple(1, fdims, NULL); H5_CHK(filespace);
  ierr = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, foff, NULL,
			       mdims, NULL); CE;
  
  hid_t dset = H5Dcreate(group, "1d", hdf5->prt_type, filespace, H5P_DEFAULT,
			 H5P_DEFAULT, H5P_DEFAULT); H5_CHK(dset);
  ierr = H5Dwrite(dset, hdf5->prt_type, memspace, filespace, dxpl, arr); CE;
  
  ierr = H5Dclose(dset); CE;
  ierr = H5Sclose(filespace); CE;
  ierr = H5Sclose(memspace); CE;
}

static void
write_idx(struct psc_output_particles *out, size_t *gidx_begin, size_t *gidx_end,
	  hid_t group, hid_t dxpl)
{
  struct psc_output_particles_hdf5 *hdf5 = to_psc_output_particles_hdf5(out);
  herr_t ierr;

  hsize_t fdims[4];
  fdims[0] = ppsc->nr_kinds;
  fdims[1] = hdf5->wdims[2];
  fdims[2] = hdf5->wdims[1];
  fdims[3] = hdf5->wdims[0];
  hid_t filespace = H5Screate_simple(4, fdims, NULL); H5_CHK(filespace);
  hid_t memspace;

  assert(sizeof(size_t) == sizeof(hsize_t));
  int rank;
  MPI_Comm_rank(psc_output_particles_comm(out), &rank);
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
// psc_output_particles_hdf5_run

static void
psc_output_particles_hdf5_run(struct psc_output_particles *out,
			      struct psc_mparticles *mprts_base)
{
  // OPT: this is not optimal in that it will convert particles to PARTICLE_TYPE twice,
  // though that only matters if particle type isn't PARTICLE_TYPE to start with.
  struct psc_output_particles_hdf5 *hdf5 = to_psc_output_particles_hdf5(out);
  MPI_Comm comm = psc_output_particles_comm(out);
  herr_t ierr;

  static int pr_A, pr_B, pr_C, pr_D, pr_E;
  if (!pr_A) {
    pr_A = prof_register("outp: local", 1., 0, 0);
    pr_B = prof_register("outp: comm", 1., 0, 0);
    pr_C = prof_register("outp: write prep", 1., 0, 0);
    pr_D = prof_register("outp: write idx", 1., 0, 0);
    pr_E = prof_register("outp: write prts", 1., 0, 0);
  }

  if (hdf5->every_step < 0 ||
      ppsc->timestep % hdf5->every_step != 0) {
    return;
  }

  mparticles_t mprts = mprts_base->get_as<mparticles_t>();

  prof_start(pr_A);
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  struct mrc_patch_info info;
  int nr_kinds = ppsc->nr_kinds;
  int *wdims = hdf5->wdims;

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
  struct hdf5_prt *arr =
    make_local_particle_array(out, mprts, map, off, idx, &n_write, &n_off, &n_total);
  prof_stop(pr_A);

  prof_start(pr_B);
  if (rank == 0) {
    int nr_global_patches;
    mrc_domain_get_nr_global_patches(ppsc->mrc_domain, &nr_global_patches);

    int *remote_sz = (int *) calloc(size, sizeof(*remote_sz));
    for (int p = 0; p < nr_global_patches; p++) {
      mrc_domain_get_global_patch_info(ppsc->mrc_domain, p, &info);
      if (info.rank == rank) { // skip local patches
	continue;
      }
      int ilo[3], ihi[3], ld[3], sz;
      find_patch_bounds(hdf5, &info, ilo, ihi, ld, &sz);
      remote_sz[info.rank] += sz;
    }

    // build global idx array, local part
    for (int p = 0; p < mprts.n_patches(); p++) {
      mrc_domain_get_local_patch_info(ppsc->mrc_domain, p, &info);
      int ilo[3], ihi[3], ld[3], sz;
      find_patch_bounds(hdf5, &info, ilo, ihi, ld, &sz);
      
      for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	  for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	    for (int kind = 0; kind < nr_kinds; kind++) {
	      int ix = jx + info.off[0] - hdf5->lo[0];
	      int iy = jy + info.off[1] - hdf5->lo[1];
	      int iz = jz + info.off[2] - hdf5->lo[2];
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
      MPI_Irecv(recv_buf[r], 2 * remote_sz[r], MPI_LONG, r, 0x4000, comm,
		&recv_req[r]);
    }
    MPI_Waitall(size, recv_req, MPI_STATUSES_IGNORE);
    free(recv_req);

    // build global idx array, remote part
    for (int p = 0; p < nr_global_patches; p++) {
      mrc_domain_get_global_patch_info(ppsc->mrc_domain, p, &info);
      if (info.rank == rank) { // skip local patches
	continue;
      }
      int ilo[3], ihi[3], ld[3], sz;
      find_patch_bounds(hdf5, &info, ilo, ihi, ld, &sz);
      for (int jz = ilo[2]; jz < ihi[2]; jz++) {
	for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	  for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	    for (int kind = 0; kind < nr_kinds; kind++) {
	      int ix = jx + info.off[0] - hdf5->lo[0];
	      int iy = jy + info.off[1] - hdf5->lo[1];
	      int iz = jz + info.off[2] - hdf5->lo[2];
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
      mrc_domain_get_local_patch_info(ppsc->mrc_domain, p, &info);
      int ilo[3], ihi[3], ld[3], sz;
      find_patch_bounds(hdf5, &info, ilo, ihi, ld, &sz);
      local_sz += sz;
    }
    size_t *l_idx = (size_t *) malloc(2 * local_sz * sizeof(*l_idx));
    size_t *l_idx_p = (size_t *) l_idx;
    for (int p = 0; p < mprts.n_patches(); p++) {
      mrc_domain_get_local_patch_info(ppsc->mrc_domain, p, &info);
      int ilo[3], ihi[3], ld[3], sz;
      find_patch_bounds(hdf5, &info, ilo, ihi, ld, &sz);
      memcpy(l_idx_p, idx[p], 2 * sz * sizeof(*l_idx_p));
      l_idx_p += 2 * sz;
    }

    MPI_Send(l_idx, 2 * local_sz, MPI_LONG, 0, 0x4000, comm);

    free(l_idx);
  }
  prof_stop(pr_B);

  prof_start(pr_C);
  char filename[strlen(hdf5->data_dir) + strlen(hdf5->basename) + 20];
  sprintf(filename, "%s/%s.%06d_p%06d.h5", hdf5->data_dir,
	  hdf5->basename, ppsc->timestep, 0);

  hid_t plist = H5Pcreate(H5P_FILE_ACCESS);

  MPI_Info mpi_info;
  MPI_Info_create(&mpi_info);
#ifdef H5_HAVE_PARALLEL
  if (hdf5->romio_cb_write) {
    MPI_Info_set(mpi_info, "romio_cb_write", hdf5->romio_cb_write);
  }
  if (hdf5->romio_ds_write) {
    MPI_Info_set(mpi_info, "romio_ds_write", hdf5->romio_ds_write);
  }
  H5Pset_fapl_mpio(plist, comm, mpi_info);
#endif

  hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist); H5_CHK(file);
  H5Pclose(plist);
  MPI_Info_free(&mpi_info);

  hid_t group = H5Gcreate(file, "particles", H5P_DEFAULT,
			  H5P_DEFAULT, H5P_DEFAULT); H5_CHK(group);
  hid_t groupp = H5Gcreate(group, "p0", H5P_DEFAULT,
			   H5P_DEFAULT, H5P_DEFAULT); H5_CHK(groupp);

  ierr = H5LTset_attribute_int(group, ".", "lo", hdf5->lo, 3); CE;
  ierr = H5LTset_attribute_int(group, ".", "hi", hdf5->hi, 3); CE;

  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER); H5_CHK(dxpl);
#ifdef H5_HAVE_PARALLEL
  if (hdf5->use_independent_io) {
    ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_INDEPENDENT); CE;
  } else {
    ierr = H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE); CE;
  }
#endif
  prof_stop(pr_C);

  prof_start(pr_D);
  write_idx(out, gidx_begin, gidx_end, groupp, dxpl);
  prof_stop(pr_D);

  prof_start(pr_E);
  write_particles(out, n_write, n_off, n_total, arr, groupp, dxpl);

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

  mprts.put_as(mprts_base, MP_DONT_COPY);
}


