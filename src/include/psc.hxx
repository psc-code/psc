
#pragma once

#include <mrc_profile.h>
#include <psc_diag.h>

#include <particles.hxx>

#include <push_particles.hxx>
#include <checks.hxx>
#include <output_particles.hxx>
#include <output_fields_c.hxx>

 #ifdef VPIC
#include "../libpsc/vpic/vpic_iface.h"
#endif

// ======================================================================
// PscParams

struct PscParams
{
  double cfl = .75;            // CFL number used to determine time step
  int nmax;                    // Number of timesteps to run
  double wallclock_limit = 0.; // Maximum wallclock time to run
  bool write_checkpoint = false;
  int write_checkpoint_every_step = 0;

  bool detailed_profiling = false; // output profiling info for each process separately
  int stats_every = 10;    // output timing and other info every so many steps
};
  
// ======================================================================
// Psc

template<typename PscConfig>
struct Psc
{
  using Mparticles_t = typename PscConfig::Mparticles_t;
  using MfieldsState = typename PscConfig::MfieldsState;
  using Balance_t = typename PscConfig::Balance_t;
  using Sort_t = typename PscConfig::Sort_t;
  using Collision_t = typename PscConfig::Collision_t;
  using PushParticles_t = typename PscConfig::PushParticles_t;
  using PushFields_t = typename PscConfig::PushFields_t;
  using Bnd_t = typename PscConfig::Bnd_t;
  using BndFields_t = typename PscConfig::BndFields_t;
  using BndParticles_t = typename PscConfig::BndParticles_t;
  using Checks_t = typename PscConfig::Checks_t;
  using Marder_t = typename PscConfig::Marder_t;

#ifdef VPIC
  using MaterialList = typename MfieldsState::MaterialList;
  using Material = typename MaterialList::Material;
  using DiagMixin = typename PscConfig::DiagMixin;
#endif
  
  // ----------------------------------------------------------------------
  // ctor

  Psc()
    : grid_{ggrid}
  {
    time_start_ = MPI_Wtime();

    // FIXME, we should use RngPool consistently throughout
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    srandom(rank);

    diag_ = psc_diag_create(MPI_COMM_WORLD);
    psc_diag_set_from_options(diag_);

    outp_ = psc_output_particles_create(MPI_COMM_WORLD);
    psc_output_particles_set_from_options(outp_);
  }

  // ----------------------------------------------------------------------
  // define_grid

  void define_grid(Grid_t::Domain& domain, GridBc& bc, Grid_t::Kinds& kinds,
		   double dt, Grid_t::NormalizationParams& norm_params)
  {
    auto coeff = Grid_t::Normalization{norm_params};
    grid_ = Grid_t::psc_make_grid(domain, bc, kinds, coeff, dt, ibn);

#ifdef VPIC
    vgrid_ = Grid::create();
    vgrid_->setup(domain.dx, dt, norm_params.cc, norm_params.eps0);

    // define the grid
    define_periodic_grid(domain.corner, domain.corner + domain.length,
			 domain.gdims, domain.np);
    

    // set field boundary conditions
    for (int p = 0; p < grid().n_patches(); p++) {
      assert(p == 0);
      for (int d = 0; d < 3; d++) {
	bool lo = grid().atBoundaryLo(p, d);
	bool hi = grid().atBoundaryHi(p, d);

	if (lo && bc.fld_lo[d] != BND_FLD_PERIODIC) {
	  Int3 bnd = {0, 0, 0};
	  bnd[d] = -1;
	  set_domain_field_bc(bnd, bc.fld_lo[d]);
	}
	
	if (hi && bc.fld_hi[d] != BND_FLD_PERIODIC) {
	  Int3 bnd = {0, 0, 0};
	  bnd[d] = 1;
	  set_domain_field_bc(bnd, bc.fld_hi[d]);
	}
      }
    }

    // set particle boundary conditions
    for (int p = 0; p < grid().n_patches(); p++) {
      assert(p == 0);
      for (int d = 0; d < 3; d++) {
	bool lo = grid().atBoundaryLo(p, d);
	bool hi = grid().atBoundaryHi(p, d);

	if (lo && bc.prt_lo[d] != BND_PRT_PERIODIC) {
	  Int3 bnd = {0, 0, 0};
	  bnd[d] = -1;
	  set_domain_particle_bc(bnd, bc.prt_lo[d]);
	}
	
	if (hi && bc.prt_hi[d] != BND_PRT_PERIODIC) {
	  Int3 bnd = {0, 0, 0};
	  bnd[d] = 1;
	  set_domain_particle_bc(bnd, bc.prt_hi[d]);
	}
      }
    }

    grid_setup_communication();
#endif
  }
  
  // ----------------------------------------------------------------------
  // define_field_array

  void define_field_array(double damp = 0.)
  {
#ifdef VPIC
    // FIXME, mv assert innto MfieldsState ctor
    assert(!material_list_.empty());

    ibn = {1, 1, 1};
#endif

#ifdef VPIC
    mflds_.reset(new MfieldsState{grid(), vgrid_, material_list_, damp});
#else
    mflds_.reset(new MfieldsState{grid()});
#endif

#ifdef VPIC
    hydro_.reset(new MfieldsHydro{grid(), vgrid_});
    interpolator_.reset(new MfieldsInterpolator{vgrid_});
    accumulator_.reset(new MfieldsAccumulator{vgrid_});
#endif
  }
  
  // ----------------------------------------------------------------------
  // init

  void init()
  {
    sort_.reset(new Sort_t{});
    pushp_.reset(new PushParticles_t{});
    pushf_.reset(new PushFields_t{});
    bnd_.reset(new Bnd_t{grid(), ibn});
    bndf_.reset(new BndFields_t{});
    bndp_.reset(new BndParticles_t{grid()});

    psc_diag_setup(diag_);
    psc_output_particles_setup(outp_);

    initialize_stats();
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

    // FIXME not quite the right place
    pr_time_step_no_comm = prof_register("time step w/o comm", 1., 0, 0);
  }
  
  // ----------------------------------------------------------------------
  // initialize

  void initialize()
  {
#ifdef VPIC
    initialize_vpic();
#else
    initialize_default();
#endif

    // initial output / stats
    mpi_printf(grid().comm(), "Performing initial diagnostics.\n");
    diagnostics();
    print_status();

    mpi_printf(grid().comm(), "Initialization complete.\n");
  }

  // ----------------------------------------------------------------------
  // integrate

  void integrate()
  {
    static int pr;
    if (!pr) {
      pr = prof_register("psc_step", 1., 0, 0);
    }

    mpi_printf(grid().comm(), "*** Advancing\n");
    double elapsed = MPI_Wtime();

    bool first_iteration = true;
    while (grid().timestep() < p_.nmax) {
      prof_start(pr);
      psc_stats_start(st_time_step);

      if (!first_iteration &&
	  p_.write_checkpoint_every_step > 0 &&
	  grid().timestep() % p_.write_checkpoint_every_step == 0) {
	//psc_write_checkpoint(psc_);
      }
      first_iteration = false;

      mpi_printf(grid().comm(), "**** Step %d / %d, Code Time %g, Wall Time %g\n", grid().timestep() + 1,
		 p_.nmax, grid().timestep() * grid().dt, MPI_Wtime() - time_start_);

      prof_start(pr_time_step_no_comm);
      prof_stop(pr_time_step_no_comm); // actual measurements are done w/ restart

      step();
      grid_->timestep_++; // FIXME, too hacky
#ifdef VPIC
      vgrid_->step++;
      assert(vgrid_->step == grid().timestep());
#endif
      
      diagnostics();

      psc_stats_stop(st_time_step);
      prof_stop(pr);

      psc_stats_val[st_nr_particles] = mprts_->get_n_prts();

      if (grid().timestep() % p_.stats_every == 0) {
	print_status();
      }

      if (p_.wallclock_limit > 0.) {
	double wallclock_elapsed = MPI_Wtime() - time_start_;
	double wallclock_elapsed_max;
	MPI_Allreduce(&wallclock_elapsed, &wallclock_elapsed_max, 1, MPI_DOUBLE, MPI_MAX,
		      MPI_COMM_WORLD);
      
	if (wallclock_elapsed_max > p_.wallclock_limit) {
	  mpi_printf(MPI_COMM_WORLD, "WARNING: Max wallclock time elapsed!\n");
	  break;
	}
      }
    }

    if (p_.write_checkpoint) {
      //psc_write_checkpoint(psc_);
    }

    // FIXME, merge with existing handling of wallclock time
    elapsed = MPI_Wtime() - elapsed;

    int  s = (int)elapsed, m  = s/60, h  = m/60, d  = h/24, w = d/ 7;
    /**/ s -= m*60,        m -= h*60, h -= d*24, d -= w*7;
    mpi_printf(grid().comm(), "*** Finished (%gs / %iw:%id:%ih:%im:%is elapsed)\n",
	       elapsed, w, d, h, m, s );
  }

#ifdef VPIC
  // ----------------------------------------------------------------------
  // step_vpic

  void step_vpic()
  {
    static int pr_sort, pr_collision, pr_checks, pr_push_prts, pr_push_flds,
      pr_bndp, pr_bndf, pr_marder, pr_inject, pr_heating;
    if (!pr_sort) {
      pr_sort = prof_register("step_sort", 1., 0, 0);
      pr_collision = prof_register("step_collision", 1., 0, 0);
      pr_push_prts = prof_register("step_push_prts", 1., 0, 0);
      pr_push_flds = prof_register("step_push_flds", 1., 0, 0);
      pr_bndp = prof_register("step_bnd_prts", 1., 0, 0);
      pr_bndf = prof_register("step_bnd_flds", 1., 0, 0);
      pr_checks = prof_register("step_checks", 1., 0, 0);
      pr_marder = prof_register("step_marder", 1., 0, 0);
      pr_inject = prof_register("step_inject", 1., 0, 0);
      pr_heating = prof_register("step_heating", 1., 0, 0);
    }

    MPI_Comm comm = grid().comm();

    // x^{n+1/2}, p^{n}, E^{n+1/2}, B^{n+1/2}

    int timestep = grid().timestep();

    auto& mprts = *mprts_;
    auto& mflds = *mflds_;

    if (balance_interval > 0 && timestep % balance_interval == 0) {
      (*balance_)(mprts);
    }

    prof_start(pr_time_step_no_comm);
    prof_stop(pr_time_step_no_comm); // actual measurements are done w/ restart

    if (sort_interval > 0 && timestep % sort_interval == 0) {
      //mpi_printf(comm, "***** Sorting...\n");
      prof_start(pr_sort);
      (*sort_)(mprts);
      prof_stop(pr_sort);
    }
    
    if (collision_->interval() > 0 && timestep % collision_->interval() == 0) {
      mpi_printf(comm, "***** Performing collisions...\n");
      prof_start(pr_collision);
      (*collision_)(mprts);
      prof_stop(pr_collision);
    }
    
    //psc_bnd_particles_open_calc_moments(psc_->bnd_particles, psc_->particles);

    checks_->continuity_before_particle_push(mprts);

    // === particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
    prof_start(pr_push_prts);
    pushp_->push_mprts(mprts, mflds, *interpolator_, *accumulator_, particle_bc_list_, num_comm_round);
    prof_stop(pr_push_prts);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}

    // field propagation B^{n+1/2} -> B^{n+1}
    pushf_->push_H(mflds, .5);
    // x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}

    prof_start(pr_bndp);
    (*bndp_)(mprts);
    prof_stop(pr_bndp);

    // field propagation E^{n+1/2} -> E^{n+3/2}

    // fill ghosts for H
    bndf_->fill_ghosts_H(mflds);
    bnd_->fill_ghosts(mflds, HX, HX + 3);
    
    // add and fill ghost for J
    bndf_->add_ghosts_J(mflds);
    bnd_->add_ghosts(mflds, JXI, JXI + 3);
    bnd_->fill_ghosts(mflds, JXI, JXI + 3);
    
    // push E
    pushf_->push_E(mflds, 1.);

    bndf_->fill_ghosts_E(mflds);
    //if (pushf_->variant == 0) {
    bnd_->fill_ghosts(mflds, EX, EX + 3);
    //}
    // x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}

    // field propagation B^{n+1} -> B^{n+3/2}
    //if (pushf_->variant == 0) {
    bndf_->fill_ghosts_E(mflds);
    bnd_->fill_ghosts(mflds, EX, EX + 3);
    //    }
    
    // push H
    pushf_->push_H(mflds, .5);

    bndf_->fill_ghosts_H(mflds);
    //if (pushf_->variant == 0) {
    bnd_->fill_ghosts(mflds, HX, HX + 3);
    //}
    // x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}

    checks_->continuity_after_particle_push(mprts, mflds);

    // E at t^{n+3/2}, particles at t^{n+3/2}
    // B at t^{n+3/2} (Note: that is not it's natural time,
    // but div B should be == 0 at any time...)
    if (marder_interval > 0 && timestep % marder_interval == 0) {
      //mpi_printf(comm, "***** Performing Marder correction...\n");
      prof_start(pr_marder);
      (*marder_)(mflds, mprts);
      prof_stop(pr_marder);
    }
    
    checks_->gauss(mprts, mflds);

#ifdef VPIC
    pushp_->prep(mprts, mflds, *interpolator_);
#endif
  }
#endif

  // ----------------------------------------------------------------------
  // step_psc

  void step_psc()
  {
    using DIM = typename PscConfig::dim_t;

    static int pr_sort, pr_collision, pr_checks, pr_push_prts, pr_push_flds,
      pr_bndp, pr_bndf, pr_marder, pr_inject_prts;
    if (!pr_sort) {
      pr_sort = prof_register("step_sort", 1., 0, 0);
      pr_collision = prof_register("step_collision", 1., 0, 0);
      pr_push_prts = prof_register("step_push_prts", 1., 0, 0);
      pr_inject_prts = prof_register("step_inject_prts", 1., 0, 0);
      pr_push_flds = prof_register("step_push_flds", 1., 0, 0);
      pr_bndp = prof_register("step_bnd_prts", 1., 0, 0);
      pr_bndf = prof_register("step_bnd_flds", 1., 0, 0);
      pr_checks = prof_register("step_checks", 1., 0, 0);
      pr_marder = prof_register("step_marder", 1., 0, 0);
    }

    // state is at: x^{n+1/2}, p^{n}, E^{n+1/2}, B^{n+1/2}
    MPI_Comm comm = grid().comm();
    int timestep = grid().timestep();

    auto& mprts = *mprts_;
    auto& mflds = *mflds_;

    if (balance_interval > 0 && timestep % balance_interval == 0) {
      (*balance_)(mprts);
    }

    if (sort_interval > 0 && timestep % sort_interval == 0) {
      mpi_printf(comm, "***** Sorting...\n");
      prof_start(pr_sort);
      (*sort_)(mprts);
      prof_stop(pr_sort);
    }
    
    if (collision_->interval() > 0 && timestep % collision_->interval() == 0) {
      mpi_printf(comm, "***** Performing collisions...\n");
      prof_start(pr_collision);
      (*collision_)(mprts);
      prof_stop(pr_collision);
    }
    
    // === particle injection
    prof_start(pr_inject_prts);
    inject_particles();
    prof_stop(pr_inject_prts);

    if (checks_->continuity_every_step > 0 && timestep % checks_->continuity_every_step == 0) {
      mpi_printf(comm, "***** Checking continuity...\n");
      prof_start(pr_checks);
      checks_->continuity_before_particle_push(mprts);
      prof_stop(pr_checks);
    }

    // === particle propagation p^{n} -> p^{n+1}, x^{n+1/2} -> x^{n+3/2}
    prof_start(pr_push_prts);
    pushp_->push_mprts(mprts, mflds);
    prof_stop(pr_push_prts);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1/2}, j^{n+1}

    // === field propagation B^{n+1/2} -> B^{n+1}
    prof_start(pr_push_flds);
    pushf_->push_H(mflds, .5, DIM{});
    prof_stop(pr_push_flds);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+1/2}, B^{n+1}, j^{n+1}

    prof_start(pr_bndp);
    (*bndp_)(mprts);
    prof_stop(pr_bndp);

    // === field propagation E^{n+1/2} -> E^{n+3/2}
    prof_start(pr_bndf);
#if 1
    bndf_->fill_ghosts_H(mflds);
    bnd_->fill_ghosts(mflds, HX, HX + 3);
#endif
    
    bndf_->add_ghosts_J(mflds);
    bnd_->add_ghosts(mflds, JXI, JXI + 3);
    bnd_->fill_ghosts(mflds, JXI, JXI + 3);
    prof_stop(pr_bndf);
    
    prof_restart(pr_push_flds);
    pushf_->push_E(mflds, 1., DIM{});
    prof_stop(pr_push_flds);
    
#if 1
    prof_restart(pr_bndf);
    bndf_->fill_ghosts_E(mflds);
    bnd_->fill_ghosts(mflds, EX, EX + 3);
    prof_stop(pr_bndf);
#endif
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+1}
      
    // === field propagation B^{n+1} -> B^{n+3/2}
    prof_restart(pr_push_flds);
    pushf_->push_H(mflds, .5, DIM{});
    prof_stop(pr_push_flds);

#if 1
    prof_start(pr_bndf);
    bndf_->fill_ghosts_H(mflds);
    bnd_->fill_ghosts(mflds, HX, HX + 3);
    prof_stop(pr_bndf);
    // state is now: x^{n+3/2}, p^{n+1}, E^{n+3/2}, B^{n+3/2}
#endif

    if (checks_->continuity_every_step > 0 && timestep % checks_->continuity_every_step == 0) {
      prof_restart(pr_checks);
      checks_->continuity_after_particle_push(mprts, mflds);
      prof_stop(pr_checks);
    }
    
    // E at t^{n+3/2}, particles at t^{n+3/2}
    // B at t^{n+3/2} (Note: that is not its natural time,
    // but div B should be == 0 at any time...)
    if (marder_interval > 0 && timestep % marder_interval == 0) {
      mpi_printf(comm, "***** Performing Marder correction...\n");
      prof_start(pr_marder);
      (*marder_)(mflds, mprts);
      prof_stop(pr_marder);
    }
    
    if (checks_->gauss_every_step > 0 && timestep % checks_->gauss_every_step == 0) {
      prof_restart(pr_checks);
      checks_->gauss(mprts, mflds);
      prof_stop(pr_checks);
    }
    
    //psc_push_particles_prep(psc->push_particles, psc->particles, psc->flds);
  }

  virtual void step()
  {
#ifdef VPIC
    step_vpic();
#else
    step_psc();
#endif
  }

  // ----------------------------------------------------------------------
  // inject_particles

  virtual void inject_particles()
  {}
  
  // ----------------------------------------------------------------------
  // define_periodic_grid
  
  void define_periodic_grid(double xl[3], double xh[3], const int gdims[3], const int np[3])
  {
#ifdef VPIC
    // SimulationMixin::setTopology(np[0], np[1], np[2]); FIXME, needed for vpic_simulation,
    // I believe only because this info is written out in diagnostics_run
    vgrid_->partition_periodic_box(xl, xh, gdims, np);
#endif
  }
  
  // ----------------------------------------------------------------------
  // set_domain_field_bc
  
  void set_domain_field_bc(Int3 bnd, int bc)
  {
#ifdef VPIC
    int boundary = BOUNDARY(bnd[0], bnd[1], bnd[2]);
    int fbc;
    switch (bc) {
    case BND_FLD_CONDUCTING_WALL: fbc = Grid::pec_fields   ; break;
    case BND_FLD_ABSORBING:       fbc = Grid::absorb_fields; break;
    default: assert(0);
    }
    vgrid_->set_fbc(boundary, fbc);
#endif
  }

  // ----------------------------------------------------------------------
  // set_domain_particle_bc
  
  void set_domain_particle_bc(Int3 bnd, int bc)
  {
#ifdef VPIC
    int boundary = BOUNDARY(bnd[0], bnd[1], bnd[2]);
    int pbc;
    switch (bc) {
    case BND_PRT_REFLECTING: pbc = Grid::reflect_particles; break;
    case BND_PRT_ABSORBING:  pbc = Grid::absorb_particles ; break;
    default: assert(0);
    }
    vgrid_->set_pbc(boundary, pbc);
#endif
  }

  // ----------------------------------------------------------------------
  // define_material

#ifdef VPIC
  Material* define_material(const char *name,
			    double eps, double mu=1.,
			    double sigma=0., double zeta=0.)
  {
    auto m = material_list_.create(name,
		       eps,   eps,   eps,
		       mu,    mu,    mu,
		       sigma, sigma, sigma,
		       zeta,  zeta,  zeta);
    return material_list_.append(m);
  }
#else
  void define_material(const char *name,
		       double eps, double mu=1.,
		       double sigma=0., double zeta=0.)
  {}
#endif

  void grid_setup_communication()
  {
#ifdef VPIC
    assert(vgrid_->nx && vgrid_->ny && vgrid_->ny);
    
    // Pre-size communications buffers. This is done to get most memory
    // allocation over with before the simulation starts running
    // FIXME, this isn't a great place. First, we shouldn't call mp
    // functions (semi-)directly. 2nd, whether we need these buffers depends
    // on b.c., which aren't yet known.
    
    // FIXME, this really isn't a good place to do this, as it requires layer breaking knowledge of
    // which communication will need the largest buffers...
    int nx1 = vgrid_->nx+1, ny1 = vgrid_->ny+1, nz1 = vgrid_->nz+1;
    vgrid_->mp_size_recv_buffer(BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_recv_buffer(BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_recv_buffer(BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_recv_buffer(BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_recv_buffer(BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_recv_buffer(BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(typename MfieldsHydro::Element));
    
    vgrid_->mp_size_send_buffer(BOUNDARY(-1, 0, 0), ny1*nz1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_send_buffer(BOUNDARY( 1, 0, 0), ny1*nz1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_send_buffer(BOUNDARY( 0,-1, 0), nz1*nx1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_send_buffer(BOUNDARY( 0, 1, 0), nz1*nx1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_send_buffer(BOUNDARY( 0, 0,-1), nx1*ny1*sizeof(typename MfieldsHydro::Element));
    vgrid_->mp_size_send_buffer(BOUNDARY( 0, 0, 1), nx1*ny1*sizeof(typename MfieldsHydro::Element));
#endif
  }

  // ----------------------------------------------------------------------
  // courant length
  
  double courant_length(const Grid_t::Domain& domain)
  {
    double inv_sum = 0.;
    for (int d = 0; d < 3; d++) {
      if (!domain.isInvar(d)) {
	inv_sum += 1. / sqr(domain.dx[d]);
      }
    }
    if (!inv_sum) { // simulation has 0 dimensions (happens in some test?)
      inv_sum = 1.;
    }
    return sqrt(1. / inv_sum);
  }

  // ----------------------------------------------------------------------
  // create_diagnotics

  void create_diagnostics(int interval)
  {
#ifdef VPIC
    diag_mixin_.diagnostics_init(interval);
#endif
  }
  
  // ----------------------------------------------------------------------
  // setup_diagnostics

  void setup_diagnostics()
  {
#ifdef VPIC
    diag_mixin_.diagnostics_setup();
#endif
  }

  // ----------------------------------------------------------------------
  // run_diagnostics
  
  void run_diagnostics()
  {
#ifdef VPIC
    diag_mixin_.diagnostics_run(*mprts_, *mflds_, *interpolator_, *hydro_, grid().domain.np);
#endif
  }

private:

  // ----------------------------------------------------------------------
  // print_profiling

  void print_profiling()
  {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (!p_.detailed_profiling) {
      prof_print_mpi(MPI_COMM_WORLD);
    } else {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      for (int i = 0; i < size; i++) {
	if (i == rank) {
	  mprintf("profile\n");
	  prof_print();
	}
	MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
  
#ifndef VPIC
  // ----------------------------------------------------------------------
  // initialize_default
  
  void initialize_default()
  {
    //pushp_.stagger(mprts, mflds); FIXME, vpic does it

    checks_->gauss(*mprts_, *mflds_);
  }
#endif

#ifdef VPIC

  // ----------------------------------------------------------------------
  // initialize_vpic
  
  void initialize_vpic()
  {
    MPI_Comm comm = grid().comm();
    
    // FIXME, just change the uses
    auto& mprts = *mprts_;
    auto& mflds = mflds_;
    // Do some consistency checks on user initialized fields
    
    mpi_printf(comm, "Checking interdomain synchronization\n");
    double err;
    TIC err = CleanDivOps::synchronize_tang_e_norm_b(*mflds_); TOC(synchronize_tang_e_norm_b, 1);
    mpi_printf(comm, "Error = %g (arb units)\n", err);
    
    mpi_printf(comm, "Checking magnetic field divergence\n");
    TIC CleanDivOps::compute_div_b_err(*mflds_); TOC(compute_div_b_err, 1);
    TIC err = CleanDivOps::compute_rms_div_b_err(*mflds_); TOC(compute_rms_div_b_err, 1);
    mpi_printf(comm, "RMS error = %e (charge/volume)\n", err);
    TIC CleanDivOps::clean_div_b(*mflds_); TOC(clean_div_b, 1);
    
    // Load fields not initialized by the user
    
    mpi_printf(comm, "Initializing radiation damping fields\n");
    TIC AccumulateOps::compute_curl_b(*mflds_); TOC(compute_curl_b, 1);
    
    mpi_printf(comm, "Initializing bound charge density\n");
    TIC CleanDivOps::clear_rhof(*mflds_); TOC(clear_rhof, 1);
    ParticlesOps::accumulate_rho_p(*mprts_, *mflds_);
    CleanDivOps::synchronize_rho(*mflds_);
    TIC AccumulateOps::compute_rhob(*mflds_); TOC(compute_rhob, 1);
    
    // Internal sanity checks
    
    mpi_printf(comm, "Checking electric field divergence\n");
    TIC CleanDivOps::compute_div_e_err(*mflds_); TOC(compute_div_e_err, 1);
    TIC err = CleanDivOps::compute_rms_div_e_err(*mflds_); TOC(compute_rms_div_e_err, 1);
    mpi_printf(comm, "RMS error = %e (charge/volume)\n", err);
    TIC CleanDivOps::clean_div_e(*mflds_); TOC(clean_div_e, 1);
    
    mpi_printf(comm, "Rechecking interdomain synchronization\n");
    TIC err = CleanDivOps::synchronize_tang_e_norm_b(*mflds_); TOC(synchronize_tang_e_norm_b, 1);
    mpi_printf(comm, "Error = %e (arb units)\n", err);
    
    mpi_printf(comm, "Uncentering particles\n");
    if (!mprts_->empty()) {
      TIC InterpolatorOps::load(*interpolator_, *mflds_); TOC(load_interpolator, 1);
      
      for (auto& sp : mprts) {
	TIC ParticlesOps::uncenter_p(&sp, *interpolator_); TOC(uncenter_p, 1);
      }
    }
  }

#endif
  
  // ----------------------------------------------------------------------
  // diagnostics

  virtual void diagnostics()
  {
#ifdef VPIC
#if 0
    TIC user_diagnostics(); TOC(user_diagnostics, 1);
#endif
#else
    // FIXME
    psc_diag_run(diag_, *mprts_, *mflds_);
    // FIXME
    (*outf_)(*mflds_, *mprts_);
#endif
    PscOutputParticlesBase{outp_}.run(*mprts_);
  }

  // ----------------------------------------------------------------------
  // print_status

  void print_status()
  {
#ifdef VPIC
#ifdef HAVE_VPIC
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  update_profile(rank == 0);
#endif
#endif
  psc_stats_log(grid().timestep());
    print_profiling();
  }

public:
  const Grid_t& grid() { return *grid_; }

protected:
  double time_start_;

  PscParams p_;
  Grid_t*& grid_;
#ifdef VPIC
  Grid* vgrid_;
#endif

#ifdef VPIC
  MaterialList material_list_;
#endif
  std::unique_ptr<MfieldsState> mflds_;
#ifdef VPIC
  std::unique_ptr<MfieldsHydro> hydro_;
  std::unique_ptr<MfieldsInterpolator> interpolator_;
  std::unique_ptr<MfieldsAccumulator> accumulator_;
  ParticleBcList particle_bc_list_;
#endif
  std::unique_ptr<Mparticles_t> mprts_;

  std::unique_ptr<Balance_t> balance_;
  std::unique_ptr<Sort_t> sort_;
  std::unique_ptr<Collision_t> collision_;
  std::unique_ptr<PushParticles_t> pushp_;
  std::unique_ptr<PushFields_t> pushf_;
  std::unique_ptr<Bnd_t> bnd_;
  std::unique_ptr<BndFields_t> bndf_;
  std::unique_ptr<BndParticles_t> bndp_;
  std::unique_ptr<Checks_t> checks_;
  std::unique_ptr<Marder_t> marder_;
  std::unique_ptr<OutputFieldsC> outf_;

#ifdef VPIC
  DiagMixin diag_mixin_;
#endif
  psc_diag* diag_;             ///< timeseries diagnostics
  psc_output_particles* outp_; ///< particle output

  // FIXME, maybe should be private
  // need to make sure derived class sets these (? -- or just leave them off by default)
  int balance_interval;
  int sort_interval;
  int marder_interval;

  int num_comm_round = {3};
  Int3 ibn = {2, 2, 2}; // FIXME!!! need to factor in invar dims (but not in vpic...)
  
  int st_nr_particles;
  int st_time_step;
};
