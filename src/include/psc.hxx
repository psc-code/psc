
#pragma once

#include <psc_method.h>
#include <mrc_profile.h>

#include <particles.hxx>

#include <push_particles.hxx>

void psc_method_vpic_initialize(struct psc_method *method, struct psc *psc,
				MparticlesBase& mprts_base); // FIXME

template<typename PscConfig>
struct Psc
{
  using Mparticles_t = typename PscConfig::Mparticles_t;

  // ----------------------------------------------------------------------
  // ctor

  Psc(psc* psc)
    : psc_(psc)
  {}

  // ----------------------------------------------------------------------
  // dtor

  ~Psc()
  {
    psc_destroy(psc_);
  }
  
  // ----------------------------------------------------------------------
  // initialize_stats
  
  void initialize_stats()
  {
    st_nr_particles = psc_stats_register("nr particles");
    st_time_step = psc_stats_register("time entire step");

    // generic stats categories
    st_time_particle = psc_stats_register("time particle update");
    st_time_field = psc_stats_register("time field update");
    st_time_comm = psc_stats_register("time communication");
    st_time_output = psc_stats_register("time output");
  }
  
  // ----------------------------------------------------------------------
  // initialize

  void initialize()
  {
    psc_view(psc_);
    mprts__->view();
    psc_mfields_view(psc_->flds);

    if (strcmp(psc_method_type(psc_->method), "vpic") == 0) {
      psc_method_vpic_initialize(psc_->method, psc_, *mprts__);
    } else {
      initialize_default(psc_->method, psc_, *mprts__);
    }

    mpi_printf(psc_comm(psc_), "Initialization complete.\n");
  }

  // ----------------------------------------------------------------------
  // integrate

  void integrate()
  {
    static int pr;
    if (!pr) {
      pr = prof_register("psc_step", 1., 0, 0);
    }

    mpi_printf(psc_comm(psc_), "*** Advancing\n");
    double elapsed = MPI_Wtime();

    bool first_iteration = true;
    while (psc_->timestep < psc_->prm.nmax) {
      prof_start(pr);
      psc_stats_start(st_time_step);

      if (!first_iteration &&
	  psc_->prm.write_checkpoint_every_step > 0 &&
	  psc_->timestep % psc_->prm.write_checkpoint_every_step == 0) {
	psc_write_checkpoint(psc_);
      }
      first_iteration = false;

      mpi_printf(psc_comm(psc_), "**** Step %d / %d, Code Time %g, Wall Time %g\n", psc_->timestep + 1,
		 psc_->prm.nmax, psc_->timestep * psc_->dt, MPI_Wtime() - psc_->time_start);

      prof_start(pr_time_step_no_comm);
      prof_stop(pr_time_step_no_comm); // actual measurements are done w/ restart

      step();
    
      psc_->timestep++; // FIXME, too hacky
      psc_method_output(psc_->method, psc_, *mprts__);

      psc_stats_stop(st_time_step);
      prof_stop(pr);

      psc_stats_val[st_nr_particles] = mprts__->get_n_prts();

      if (psc_->timestep % psc_->prm.stats_every == 0) {
	psc_stats_log(psc_);
	psc_print_profiling(psc_);
      }

      if (psc_->prm.wallclock_limit > 0.) {
	double wallclock_elapsed = MPI_Wtime() - psc_->time_start;
	double wallclock_elapsed_max;
	MPI_Allreduce(&wallclock_elapsed, &wallclock_elapsed_max, 1, MPI_DOUBLE, MPI_MAX,
		      MPI_COMM_WORLD);
      
	if (wallclock_elapsed_max > psc_->prm.wallclock_limit) {
	  mpi_printf(MPI_COMM_WORLD, "WARNING: Max wallclock time elapsed!\n");
	  break;
	}
      }
    }

    if (psc_->prm.write_checkpoint) {
      psc_write_checkpoint(psc_);
    }

    // FIXME, merge with existing handling of wallclock time
    elapsed = MPI_Wtime() - elapsed;

    int  s = (int)elapsed, m  = s/60, h  = m/60, d  = h/24, w = d/ 7;
    /**/ s -= m*60,        m -= h*60, h -= d*24, d -= w*7;
    mpi_printf(psc_comm(psc_), "*** Finished (%gs / %iw:%id:%ih:%im:%is elapsed)\n",
	       elapsed, w, d, h, m, s );
  }

  virtual void step() = 0;

  psc* get_psc() const { return psc_; } // FIXME, should go away

private:

  // ----------------------------------------------------------------------
  // initialize_default
  
  static void initialize_default(struct psc_method *method, struct psc *psc,
				 MparticlesBase& mprts)
  {
    auto pushp = PscPushParticlesBase{psc->push_particles};
    pushp.stagger(mprts, PscMfieldsBase{psc->flds});
    
    // initial output / stats
    psc_method_output(psc->method, psc, mprts);
    psc_stats_log(psc);
    psc_print_profiling(psc);
  }

protected:
  MparticlesBase* mprts__;
  psc* psc_;

  int st_nr_particles;
  int st_time_step;
};
