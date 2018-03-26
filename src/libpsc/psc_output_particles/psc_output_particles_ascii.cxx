
#include "psc_particles_double.h"
#include "output_particles.hxx"

#include <mrc_params.h>
#include <string.h>

#include <string.h>

#define to_psc_output_particles_ascii(out) \
  mrc_to_subobj(out, struct psc_output_particles_ascii)

struct psc_output_particles_ascii : PscOutputParticlesParams, OutputParticlesBase
{
  psc_output_particles_ascii(const PscOutputParticlesParams& params)
    : PscOutputParticlesParams(params)
  {}
};

// ----------------------------------------------------------------------
// psc_output_particles_ascii_setup

static void
psc_output_particles_ascii_setup(struct psc_output_particles *out)
{
  new(out->obj.subctx) psc_output_particles_ascii{out->params};
}

// ----------------------------------------------------------------------
// psc_output_particles_ascii_destroy

static void
psc_output_particles_ascii_destroy(struct psc_output_particles *out)
{
  struct psc_output_particles_ascii *sub = to_psc_output_particles_ascii(out);
  sub->~psc_output_particles_ascii();
}

// ----------------------------------------------------------------------
// psc_output_particles_ascii_run

static void
psc_output_particles_ascii_run(struct psc_output_particles *out,
			       struct psc_mparticles *_mprts_base)
{
  auto mprts_base = PscMparticlesBase{_mprts_base};
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

  auto mprts = mprts_base.get_as<PscMparticlesDouble>();

  FILE *file = fopen(filename, "w");
  for (int p = 0; p < mprts->n_patches(); p++) {
    int n = 0;
    for (auto& prt : mprts[p]) {
      fprintf(file, "%d %g %g %g %g %g %g %g %d\n",
	      n, prt.xi, prt.yi, prt.zi,
	      prt.pxi, prt.pyi, prt.pzi,
	      prt.qni_wni_, prt.kind());
      n++;
    }
  }

  mprts.put_as(mprts_base, MP_DONT_COPY);
  
  fclose(file);
}

// ======================================================================
// psc_output_particles: subclass "ascii"

struct psc_output_particles_ops_ascii : psc_output_particles_ops {
  using Wrapper_t = OutputParticlesWrapper<psc_output_particles_ascii>;
  psc_output_particles_ops_ascii() {
    name                  = "ascii";
    size                  = Wrapper_t::size;
    setup                 = Wrapper_t::setup;
    destroy               = Wrapper_t::destroy;
    run                   = psc_output_particles_ascii_run;
  }
} psc_output_particles_ascii_ops;
