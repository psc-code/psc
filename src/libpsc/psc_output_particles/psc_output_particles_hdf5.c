
#include "psc_output_particles_private.h"

#include <mrc_params.h>

#include <hdf5.h>
#include <hdf5_hl.h>

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
// psc_output_particles_hdf5_setup

static void
psc_output_particles_hdf5_setup(struct psc_output_particles *out)
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
get_sort_index(int p, const particle_c_t *part)
{
  struct psc_patch *patch = &ppsc->patch[p];
  particle_c_real_t dxi = 1.f / ppsc->dx[0];
  particle_c_real_t dyi = 1.f / ppsc->dx[1];
  particle_c_real_t dzi = 1.f / ppsc->dx[2];
  int *ldims = patch->ldims;
  
  particle_c_real_t u = (part->xi - patch->xb[0]) * dxi;
  particle_c_real_t v = (part->yi - patch->xb[1]) * dyi;
  particle_c_real_t w = (part->zi - patch->xb[2]) * dzi;
  int j0 = particle_c_real_fint(u);
  int j1 = particle_c_real_fint(v);
  int j2 = particle_c_real_fint(w);
  assert(j0 >= 0 && j0 < ldims[0]);
  assert(j1 >= 0 && j1 < ldims[1]);
  assert(j2 >= 0 && j2 < ldims[2]);

  int kind;
  if (part->qni < 0.) {
    kind = 0; // electron
  } else if (part->qni > 0.) {
    kind = 1; // ion
  } else {
    kind = 2; // neutral
  }
  assert(kind < ppsc->prm.nr_kinds);
 
  return cell_index_3_to_1(ldims, j0, j1, j2) * ppsc->prm.nr_kinds + kind;
}

// ----------------------------------------------------------------------
// count_sort

static void
count_sort(mparticles_c_t *particles, int **off, int **map)
{
  int nr_kinds = ppsc->prm.nr_kinds;

  for (int p = 0; p < particles->nr_patches; p++) {
    int *ldims = ppsc->patch[p].ldims;
    int nr_indices = ldims[0] * ldims[1] * ldims[2] * nr_kinds;
    off[p] = calloc(nr_indices + 1, sizeof(*off[p]));
    particles_c_t *pp = psc_mparticles_get_patch_c(particles, p);

    // counting sort to get map 
    for (int n = 0; n < pp->n_part; n++) {
      particle_c_t *part = particles_c_get_one(pp, n);
      int si = get_sort_index(p, part);
      off[p][si]++;
    }
    // prefix sum to get offsets
    int o = 0;
    int *off2 = malloc((nr_indices + 1) * sizeof(*off2));
    for (int si = 0; si <= nr_indices; si++) {
      int cnt = off[p][si];
      off[p][si] = o; // this will be saved for later
      off2[si] = o; // this one will overwritten when making the map
      o += cnt;
    }

    // sort a map only, not the actual particles
    map[p] = malloc(pp->n_part * sizeof(*map[p]));
    for (int n = 0; n < pp->n_part; n++) {
      particle_c_t *part = particles_c_get_one(pp, n);
      int si = get_sort_index(p, part);
      map[p][off2[si]++] = n;
    }
    free(off2);
  }  
}

// ----------------------------------------------------------------------
// make_local_particle_array

static struct hdf5_prt *
make_local_particle_array(struct psc_output_particles *out,
			  mparticles_c_t *particles, int **off, int **map,
			  int *idx_begin, int *idx_end, int *rdims, int *p_n_write)
{
  struct psc_output_particles_hdf5 *hdf5 = to_psc_output_particles_hdf5(out);
  int nr_kinds = ppsc->prm.nr_kinds;

  // count all particles to be written locally
  int n_write = 0;
  for (int p = 0; p < particles->nr_patches; p++) {
    int *ldims = ppsc->patch[p].ldims, *loff = ppsc->patch[p].off;

    int ilo[3], ihi[3];
    for (int d = 0; d < 3; d++) {
      ilo[d] = MAX(0, hdf5->lo[d] - loff[d]);
      ihi[d] = MIN(ldims[d], hdf5->hi[d] - loff[d]);
    }

    for (int jz = ilo[2]; jz < ihi[2]; jz++) {
      for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	  for (int kind = 0; kind < nr_kinds; kind++) {
	    int si = cell_index_3_to_1(ldims, jx, jy, jz) * nr_kinds + kind;
	    n_write += off[p][si+1] - off[p][si];
	  }
	}
      }
    }
  }

  struct hdf5_prt *arr = malloc(n_write * sizeof(*arr));

  // copy particles to be written into temp array
  int nn = 0;
  for (int p = 0; p < particles->nr_patches; p++) {
    particles_c_t *pp = psc_mparticles_get_patch_c(particles, p);
    int *ldims = ppsc->patch[p].ldims, *loff = ppsc->patch[p].off;

    int ilo[3], ihi[3];
    for (int d = 0; d < 3; d++) {
      ilo[d] = MAX(0, hdf5->lo[d] - loff[d]);
      ihi[d] = MIN(ldims[d], hdf5->hi[d] - loff[d]);
    }

    for (int jz = ilo[2]; jz < ihi[2]; jz++) {
      for (int jy = ilo[1]; jy < ihi[1]; jy++) {
	for (int jx = ilo[0]; jx < ihi[0]; jx++) {
	  for (int kind = 0; kind < nr_kinds; kind++) {
	    int si = cell_index_3_to_1(ldims, jx, jy, jz) * nr_kinds + kind;
	    int ii = ((kind * rdims[2] + jz - ilo[2])
		      * rdims[1] + jy - ilo[1]) * rdims[0] + jx - ilo[0];
	    idx_begin[ii] = nn;
	    idx_end[ii] = nn + off[p][si+1] - off[p][si];
	    for (int n = off[p][si]; n < off[p][si+1]; n++, nn++) {
	      particle_c_t *part = particles_c_get_one(pp, map[p][n]);
	      arr[nn].x  = part->xi;
	      arr[nn].y  = part->yi;
	      arr[nn].z  = part->zi;
	      arr[nn].px = part->pxi;
	      arr[nn].py = part->pyi;
	      arr[nn].pz = part->pzi;
	      arr[nn].q  = part->qni;
	      arr[nn].m  = part->mni;
	      arr[nn].w  = part->wni;
	    }
	  }
	}
      }
    }
  }

  *p_n_write = n_write;
  return arr;
}

// ----------------------------------------------------------------------
// psc_output_particles_hdf5_run

static void
psc_output_particles_hdf5_run(struct psc_output_particles *out,
			       mparticles_base_t *particles_base)
{
  struct psc_output_particles_hdf5 *hdf5 = to_psc_output_particles_hdf5(out);
  int ierr;

  if (hdf5->every_step < 0 ||
      ppsc->timestep % hdf5->every_step != 0) {
    return;
  }

  mparticles_c_t *particles = psc_mparticles_get_c(particles_base, 0);
  int nr_kinds = ppsc->prm.nr_kinds;

  // set hi to gdims by default (if not set differently before)
  // and calculate rdims (global dims of region we're writing)
  int rdims[3];
  for (int d = 0; d < 3; d++) {
    assert(hdf5->lo[d] >= 0);
    if (hdf5->hi[d] == 0) {
      hdf5->hi[d] = ppsc->domain.gdims[d];
    }
    assert(hdf5->hi[d] <= ppsc->domain.gdims[d]);
    rdims[d] = hdf5->hi[d] - hdf5->lo[d];
  }

  int **off = malloc(particles->nr_patches * sizeof(*off));
  int **map = malloc(particles->nr_patches * sizeof(*off));

  count_sort(particles, map, off);

  // alloc
  int *idx_begin = malloc(nr_kinds * rdims[0] * rdims[1] * rdims[2] * sizeof(*idx_begin));
  int *idx_end   = malloc(nr_kinds * rdims[0] * rdims[1] * rdims[2] * sizeof(*idx_end));

  int n_write;
  struct hdf5_prt *arr =
    make_local_particle_array(out, particles, map, off, idx_begin, idx_end, rdims,
			      &n_write);

  int rank;
  MPI_Comm_rank(psc_output_particles_comm(out), &rank);
  char filename[strlen(hdf5->data_dir) + strlen(hdf5->basename) + 20];
  sprintf(filename, "%s/%s.%06d_p%06d.h5", hdf5->data_dir,
	  hdf5->basename, ppsc->timestep, rank);

  hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
  MPI_Info info;
  MPI_Info_create(&info);
#ifdef H5_HAVE_PARALLEL
  if (hdf5->romio_cb_write) {
    MPI_Info_set(info, "romio_cb_write", hdf5->romio_cb_write);
  }
  if (hdf5->romio_ds_write) {
    MPI_Info_set(info, "romio_ds_write", hdf5->romio_ds_write);
  }
  H5Pset_fapl_mpio(plist, psc_output_particles_comm(out), info);
#endif

  hid_t file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist); H5_CHK(file);
  H5Pclose(plist);
  MPI_Info_free(&info);

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

  hsize_t fdims[4] = { nr_kinds, rdims[2], rdims[1], rdims[0] };
  hid_t filespace = H5Screate_simple(4, fdims, NULL); H5_CHK(filespace);

  hid_t dset = H5Dcreate(groupp, "idx_begin", H5T_NATIVE_INT, filespace,
			 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(dset);
  ierr = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, filespace, dxpl, idx_begin); CE;
  ierr = H5Dclose(dset); CE;

  dset = H5Dcreate(groupp, "idx_end", H5T_NATIVE_INT, filespace,
			 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT); H5_CHK(dset);
  ierr = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, filespace, dxpl, idx_end); CE;
  ierr = H5Dclose(dset); CE;

  ierr = H5Sclose(filespace); CE;

  if (n_write > 0) {
    hsize_t fdims[1] = { n_write };
    hid_t filespace = H5Screate_simple(1, fdims, NULL); H5_CHK(filespace);
    hid_t dset = H5Dcreate(groupp, "1d", hdf5->prt_type, filespace, H5P_DEFAULT,
			   H5P_DEFAULT, H5P_DEFAULT); H5_CHK(dset);
    ierr = H5Dwrite(dset, hdf5->prt_type, H5S_ALL, filespace, dxpl, arr); CE;
    
    ierr = H5Dclose(dset); CE;
    ierr = H5Sclose(filespace); CE;
  }

  ierr = H5Pclose(dxpl); CE;

  H5Gclose(groupp);
  H5Gclose(group);
  H5Fclose(file);

  free(arr);
  free(idx_begin);
  free(idx_end);

  for (int p = 0; p < particles->nr_patches; p++) {
    free(off[p]);
    free(map[p]);
  }
  free(off);
  free(map);
  
  psc_mparticles_put_c(particles, particles_base);
}

// ======================================================================
// psc_output_particles: subclass "hdf5"

struct psc_output_particles_ops psc_output_particles_hdf5_ops = {
  .name                  = "hdf5",
  .size                  = sizeof(struct psc_output_particles_hdf5),
  .param_descr           = psc_output_particles_hdf5_descr,
  .setup                 = psc_output_particles_hdf5_setup,
  .destroy               = psc_output_particles_hdf5_destroy,
  .run                   = psc_output_particles_hdf5_run,
};
