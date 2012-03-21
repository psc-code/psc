
#include "psc_output_particles_private.h"

#include <mrc_params.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define to_psc_output_particles_hdf5(out) \
  mrc_to_subobj(out, struct psc_output_particles_hdf5)

struct psc_output_particles_hdf5 {
  const char *data_dir;
  const char *basename;
  int every_step;
  int lo[3];
  int hi[3];
};

#define VAR(x) (void *)offsetof(struct psc_output_particles_hdf5, x)
static struct param psc_output_particles_hdf5_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(".")       },
  { "basename"           , VAR(basename)             , PARAM_STRING("prt")     },
  { "every_step"         , VAR(every_step)           , PARAM_INT(-1)           },
  { "lo"                 , VAR(lo)                   , PARAM_INT3(0, 0, 0)     },
  { "hi"                 , VAR(hi)                   , PARAM_INT3(0, 0, 0)     },
  {},
};
#undef VAR

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
// psc_output_particles_hdf5_run

static void
psc_output_particles_hdf5_run(struct psc_output_particles *out,
			       mparticles_base_t *particles_base)
{
  struct psc_output_particles_hdf5 *hdf5 = to_psc_output_particles_hdf5(out);

  if (hdf5->every_step < 0 ||
      ppsc->timestep % hdf5->every_step != 0) {
    return;
  }

  mparticles_c_t *particles = psc_mparticles_get_c(particles_base, 0);

  int **off = malloc(particles->nr_patches * sizeof(*off));
  int **map = malloc(particles->nr_patches * sizeof(*off));
  for (int p = 0; p < particles->nr_patches; p++) {
    int *ldims = ppsc->patch[p].ldims;
    int nr_indices = ldims[0] * ldims[1] * ldims[2] * ppsc->prm.nr_kinds;
    off[p] = calloc(nr_indices + 1, sizeof(*off[p]));
    particles_c_t *pp = psc_mparticles_get_patch_c(particles, p);

    // counting sort to get map 
    for (int n = 0; n < pp->n_part; n++) {
      particle_c_t *part = particles_c_get_one(pp, n);
      int ci = get_sort_index(p, part);
      off[p][ci]++;
    }
    // prefix sum to get offsets
    int o = 0;
    int *off2 = malloc((nr_indices + 1) * sizeof(*off2));
    for (int ci = 0; ci <= nr_indices; ci++) {
      int cnt = off[p][ci];
      off[p][ci] = o; // this will be saved for later
      off2[ci] = o; // this one will overwritten when making the map
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

  int rank;
  MPI_Comm_rank(psc_output_particles_comm(out), &rank);
  char filename[strlen(hdf5->data_dir) + strlen(hdf5->basename) + 20];
  sprintf(filename, "%s/%s.%06d_p%06d.h5", hdf5->data_dir,
	  hdf5->basename, ppsc->timestep, rank);

  for (int d = 0; d < 3; d++) {
    assert(hdf5->lo[d] >= 0);
    if (hdf5->hi[d] == 0) {
      hdf5->hi[d] = ppsc->domain.gdims[d];
    }
    assert(hdf5->hi[d] <= ppsc->domain.gdims[d]);
  }

  FILE *file = fopen(filename, "w");
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
	  for (int kind = 0; kind < ppsc->prm.nr_kinds; kind++) {
	    int ci = cell_index_3_to_1(ldims, jx, jy, jz);
	    int si = ci * ppsc->prm.nr_kinds + kind;
	    for (int n = off[p][si]; n < off[p][si+1]; n++) {
	      particle_c_t *part = particles_c_get_one(pp, map[p][n]);
	      fprintf(file, "%d %g %g %g %g %g %g %g %g %g %d %d\n",
		      nn, part->xi, part->yi, part->zi,
		      part->pxi, part->pyi, part->pzi,
		      part->qni, part->mni, part->wni, ci, kind);
	      nn++;
	    }
	  }
	}
      }
    }
  }
  fclose(file);
  
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
  .run                   = psc_output_particles_hdf5_run,
};
