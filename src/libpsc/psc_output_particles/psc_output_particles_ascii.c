
#include "psc_output_particles_private.h"
#include "psc_particles_as_double.h"

#include <mrc_params.h>
#include <string.h>

#include <string.h>

#define to_psc_output_particles_ascii(out) \
  mrc_to_subobj(out, struct psc_output_particles_ascii)

struct psc_output_particles_ascii {
  const char *data_dir;
  const char *basename;
  int every_step;
};

#define VAR(x) (void *)offsetof(struct psc_output_particles_ascii, x)
static struct param psc_output_particles_ascii_descr[] = {
  { "data_dir"           , VAR(data_dir)             , PARAM_STRING(".")       },
  { "basename"           , VAR(basename)             , PARAM_STRING("prt")     },
  { "every_step"         , VAR(every_step)           , PARAM_INT(-1)           },
  {},
};
#undef VAR

// ----------------------------------------------------------------------
// psc_output_particles_ascii_run

static void
psc_output_particles_ascii_run(struct psc_output_particles *out,
			       struct psc_mparticles *mprts_base)
{
  struct psc_output_particles_ascii *asc = to_psc_output_particles_ascii(out);

  if (asc->every_step < 0 ||
      ppsc->timestep % asc->every_step != 0) {
    return;
  }

  int rank;
  MPI_Comm_rank(psc_output_particles_comm(out), &rank);
  char filename[strlen(asc->data_dir) + strlen(asc->basename) + 19];
  sprintf(filename, "%s/%s.%06d_p%06d.asc", asc->data_dir,
	  asc->basename, ppsc->timestep, rank);

  struct psc_mparticles *mprts = psc_mparticles_get_as(mprts_base, "c", 0);

  FILE *file = fopen(filename, "w");
  for (int p = 0; p < mprts->nr_patches; p++) {
    particle_range_t prts = particle_range_mprts(mprts, p);
    for (int n = 0; n < particle_range_size(prts); n++) {
      particle_t *part = particle_iter_at(prts.begin, n);
      fprintf(file, "%d %g %g %g %g %g %g %g %d\n",
	      n, part->xi, part->yi, part->zi,
	      part->pxi, part->pyi, part->pzi,
	      part->qni_wni, part->kind);
    }
  }

  psc_mparticles_put_as(mprts, mprts_base, MP_DONT_COPY);
  
  fclose(file);
}

// ======================================================================
// psc_output_particles: subclass "ascii"

struct psc_output_particles_ops psc_output_particles_ascii_ops = {
  .name                  = "ascii",
  .size                  = sizeof(struct psc_output_particles_ascii),
  .param_descr           = psc_output_particles_ascii_descr,
  .run                   = psc_output_particles_ascii_run,
};
