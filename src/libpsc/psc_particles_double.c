
#include "psc.h"
#include "psc_particles_as_double.h"
#include "psc_particles_inc.h"
#include "psc_particles_c.h"

// ======================================================================
// psc_mparticles: subclass "double"
  
// ----------------------------------------------------------------------
// conversion to/from "c"

static inline void
calc_vxi(particle_double_real_t vxi[3], particle_double_t *part)
{
  particle_double_real_t root =
    1.f / sqrt(1.f + sqr(part->pxi) + sqr(part->pyi) + sqr(part->pzi));
  vxi[0] = part->pxi * root;
  vxi[1] = part->pyi * root;
  vxi[2] = part->pzi * root;
}

static void
get_particle_c(particle_double_t *prt, int n, struct psc_particles *prts_c)
{
  particle_double_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  particle_c_t *prt_c = particles_c_get_one(prts_c, n);

  prt->xi      = prt_c->xi;
  prt->yi      = prt_c->yi;
  prt->zi      = prt_c->zi;
  prt->pxi     = prt_c->pxi;
  prt->pyi     = prt_c->pyi;
  prt->pzi     = prt_c->pzi;
  prt->kind    = prt_c->kind;
  prt->qni_wni = prt_c->qni * prt_c->wni;

  particle_double_real_t vxi[3];
  calc_vxi(vxi, prt);
  prt->xi += dth[0] * vxi[0];
  prt->yi += dth[1] * vxi[1];
  prt->zi += dth[2] * vxi[2];
}

static void
put_particle_c(particle_double_t *prt, int n, struct psc_particles *prts_c)
{
  particle_double_real_t dth[3] = { .5 * ppsc->dt, .5 * ppsc->dt, .5 * ppsc->dt };
  // don't shift in invariant directions
  for (int d = 0; d < 3; d++) {
    if (ppsc->domain.gdims[d] == 1) {
      dth[d] = 0.;
    }
  }
  
  particle_double_real_t vxi[3];
  calc_vxi(vxi, prt);

  particle_c_t *prt_c = particles_c_get_one(prts_c, n);

  particle_c_real_t qni = ppsc->kinds[prt->kind].q;
  particle_c_real_t mni = ppsc->kinds[prt->kind].m;
  particle_c_real_t wni = prt->qni_wni / qni;

  prt_c->xi      = prt->xi - dth[0] * vxi[0];
  prt_c->yi      = prt->yi - dth[1] * vxi[1];
  prt_c->zi      = prt->zi - dth[2] * vxi[2];
  prt_c->pxi     = prt->pxi;
  prt_c->pyi     = prt->pyi;
  prt_c->pzi     = prt->pzi;
  prt_c->kind    = prt->kind;
  prt_c->qni     = qni;
  prt_c->wni     = wni;
  prt_c->mni     = mni;
}

static void
psc_mparticles_double_copy_to_c(int p, struct psc_mparticles *mprts,
				struct psc_mparticles *mprts_c, unsigned int flags)
{
  psc_mparticles_copy_to(p, mprts, mprts_c, flags, put_particle_c);
}

static void
psc_mparticles_double_copy_from_c(int p, struct psc_mparticles *mprts,
				  struct psc_mparticles *mprts_c, unsigned int flags)
{
  psc_mparticles_copy_from(p, mprts, mprts_c, flags, get_particle_c);
}

static struct mrc_obj_method psc_particles_double_methods[] = {
  MRC_OBJ_METHOD("copy_to_c",   psc_mparticles_double_copy_to_c),
  MRC_OBJ_METHOD("copy_from_c", psc_mparticles_double_copy_from_c),
  {}
};

#include "psc_particles_common.c"

#ifdef HAVE_LIBHDF5_HL

// ----------------------------------------------------------------------
// psc_mparticles_double_read

static void
psc_mparticles_double_read(struct psc_mparticles *mprts, struct mrc_io *io)
{
  int ierr;
  long h5_file;
  mrc_io_get_h5_file(io, &h5_file);

  hid_t group = H5Gopen(h5_file, mrc_io_obj_path(io, mprts), H5P_DEFAULT); H5_CHK(group);
  // FIXME those should be superclass bits
  mprts->domain = mrc_io_read_ref(io, mprts, "domain", mrc_domain);
  mrc_domain_get_patches(mprts->domain, &mprts->nr_patches);
  mrc_io_read_int(io, mprts, "flags", (int *) &mprts->flags);

  mprts->prts = calloc(mprts->nr_patches, sizeof(*mprts->prts));
  mprts->mpatch = calloc(mprts->nr_patches, sizeof(*mprts->mpatch));

  for (int p = 0; p < mprts->nr_patches; p++) {
    char name[20]; sprintf(name, "prts%d", p);
    mprts->prts[p] = psc_particles_create(MPI_COMM_NULL);
    psc_particles_set_type(mprts->prts[p], "double");
    psc_particles_set_name(mprts->prts[p], name);
    mprts->prts[p]->p = p;

    particle_range_t prts = particle_range_mprts(mprts, p);
    char pname[10]; sprintf(pname, "p%d", p);
    hid_t pgroup = H5Gopen(group, pname, H5P_DEFAULT); H5_CHK(pgroup);
    int n_prts;
    ierr = H5LTget_attribute_int(pgroup, ".", "n_prts", &n_prts); CE;
    psc_particles_set_n_prts(mprts->prts[p], n_prts);
    psc_particles_setup(mprts->prts[p]);
    
    if (n_prts > 0) {
      ierr = H5LTread_dataset_double(pgroup, "data",
				    (double *) particle_iter_deref(prts.begin)); CE;
    }
    ierr = H5Gclose(pgroup); CE;
  }
  ierr = H5Gclose(group); CE;
  psc_mparticles_view(mprts);
}

#endif

struct psc_mparticles_ops psc_mparticles_double_ops = {
  .name                    = "double",
  .methods                 = psc_particles_double_methods,
#ifdef HAVE_LIBHDF5_HL
  .write                   = psc_mparticles_double_write,
  .read                    = psc_mparticles_double_read,
#endif
};

