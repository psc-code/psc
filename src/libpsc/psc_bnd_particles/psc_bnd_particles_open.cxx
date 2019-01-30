
#include "fields.hxx"

using Fields = Fields3d<fields_t>;
using real_t = mparticles_t::real_t;

static const int debug_every_step = 10;

static inline bool at_lo_boundary(int p, int d);
static inline bool at_hi_boundary(int p, int d);

static inline double
random1()
{
  return random() / (double) RAND_MAX;
}

static void
copy_to_mrc_fld(struct mrc_fld *m3, struct psc_mfields *mflds)
{
  mfields_t mf(mflds);
  psc_foreach_patch(ppsc, p) {
    Fields F(mf[p]);
    struct mrc_fld_patch *m3p = mrc_fld_patch_get(m3, p);
    mrc_fld_foreach(m3, ix,iy,iz, 0,0) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	MRC_M3(m3p, m, ix,iy,iz) = F(m, ix,iy,iz);
      }
    } mrc_fld_foreach_end;
    mrc_fld_patch_put(m3);
  }
}

static void
debug_dump(struct mrc_io *io, struct psc_mfields *mflds)
{
  if (ppsc->timestep % debug_every_step != 0) {
    return;
  }

  struct mrc_fld *mrc_fld = mrc_domain_m3_create(ppsc->mrc_domain);
  mrc_fld_set_name(mrc_fld, psc_mfields_name(mflds));
  mrc_fld_set_param_int(mrc_fld, "nr_ghosts", 2);
  mrc_fld_set_param_int(mrc_fld, "nr_comps", mflds->nr_fields);
  mrc_fld_setup(mrc_fld);
  for (int m = 0; m < mflds->nr_fields; m++) {
    mrc_fld_set_comp_name(mrc_fld, m, psc_mfields_comp_name(mflds, m));
  }
  copy_to_mrc_fld(mrc_fld, mflds);
  mrc_fld_write(mrc_fld, io);
  mrc_fld_destroy(mrc_fld);
}

static void
average_9_point(struct psc_mfields *mflds_av, struct psc_mfields *mflds)
{
  const int b = 1;
  mfields_t mf_av(mflds_av), mf(mflds);
  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    Fields F(mf[p]), F_av(mf_av[p]);
    
    if (at_lo_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int izb = 0, ize = ppatch->ldims[2];
	if (at_lo_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] != BND_PART_PERIODIC) {
	  izb++;
	}
	if (at_hi_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] != BND_PART_PERIODIC) {
	  ize--;
	}
	for (int iz = izb; iz < ize; iz++) {
	  F_av(m, 0,0,iz) = (F(m, 0,b+0,iz-1) +
			     F(m, 0,b+1,iz-1) +
			     F(m, 0,b+2,iz-1) +
			     F(m, 0,b+0,iz  ) +
			     F(m, 0,b+1,iz  ) +
			     F(m, 0,b+2,iz  ) +
			     F(m, 0,b+0,iz+1) +
			     F(m, 0,b+1,iz+1) +
			     F(m, 0,b+2,iz+1)) / 9.;
	  //F_av(m, 0,0,iz) = F(m, 0,b,iz);
	}
      }
    }

    if (at_hi_boundary(p, 1) && ppsc->domain.bnd_part_hi[1] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int izb = 0, ize = ppatch->ldims[2];
	if (at_lo_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] != BND_PART_PERIODIC) {
	  izb++;
	}
	if (at_hi_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] != BND_PART_PERIODIC) {
	  ize--;
	}
	int iy = ppatch->ldims[1] - 1;
	for (int iz = izb; iz < ize; iz++) {
	  F_av(m, 0,iy,iz) = (F(m, 0,iy-b-0,iz-1) +
			      F(m, 0,iy-b-1,iz-1) +
			      F(m, 0,iy-b-2,iz-1) +
			      F(m, 0,iy-b-0,iz  ) +
			      F(m, 0,iy-b-1,iz  ) +
			      F(m, 0,iy-b-2,iz  ) +
			      F(m, 0,iy-b-0,iz+1) +
			      F(m, 0,iy-b-1,iz+1) +
			      F(m, 0,iy-b-2,iz+1)) / 9.;
	  //F_av(m, 0,iy,iz) = F(m, 0,iy-b,iz);
	}
      }
    }

    if (at_lo_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int iyb = 0, iye = ppatch->ldims[1];
	if (at_lo_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] != BND_PART_PERIODIC) {
	  iyb++;
	}
	if (at_hi_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] != BND_PART_PERIODIC) {
	  iye--;
	}
	for (int iy = iyb; iy < iye; iy++) {
	  F_av(m, 0,iy,0) = (F(m, 0,iy-1,b+0) +
			     F(m, 0,iy-1,b+1) +
			     F(m, 0,iy-1,b+2) +
			     F(m, 0,iy  ,b+0) +
			     F(m, 0,iy  ,b+1) +
			     F(m, 0,iy  ,b+2) +
			     F(m, 0,iy+1,b+0) +
			     F(m, 0,iy+1,b+1) +
			     F(m, 0,iy+1,b+2)) / 9.;
	  //	  F_av(m, 0,iy,0) = F(m, 0,iy,b);
	}
      }
    }

    if (at_hi_boundary(p, 2) && ppsc->domain.bnd_part_hi[2] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int iz = ppatch->ldims[2] - 1;
	int iyb = 0, iye = ppatch->ldims[1];
	if (at_lo_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] != BND_PART_PERIODIC) {
	  iyb++;
	}
	if (at_hi_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] != BND_PART_PERIODIC) {
	  iye--;
	}
	for (int iy = iyb; iy < iye; iy++) {
	  F_av(m, 0,iy,iz) = (F(m, 0,iy-1,iz-b-0) +
				       F(m, 0,iy-1,iz-b-1) +
				       F(m, 0,iy-1,iz-b-2) +
				       F(m, 0,iy  ,iz-b-0) +
				       F(m, 0,iy  ,iz-b-1) +
				       F(m, 0,iy  ,iz-b-2) +
				       F(m, 0,iy+1,iz-b-0) +
				       F(m, 0,iy+1,iz-b-1) +
				       F(m, 0,iy+1,iz-b-2)) / 9.;
	  //	  F_av(m, 0,iy,iz) = F(m, 0,iy,iz-b);
	}
      }
    }

    if (at_lo_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] == BND_PART_OPEN &&
	at_lo_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int iy = 0;
	int iz = 0;
	F_av(m, 0,iy,iz) = (F(m, 0,iy+1,iz+1) +
				     F(m, 0,iy+1,iz+2) +
				     F(m, 0,iy+1,iz+3) +
				     F(m, 0,iy+2,iz+1) +
				     F(m, 0,iy+2,iz+2) +
				     F(m, 0,iy+2,iz+3) +
				     F(m, 0,iy+3,iz+1) +
				     F(m, 0,iy+3,iz+2) +
				     F(m, 0,iy+3,iz+3)) / 9.;
      }
    }

    if (at_lo_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] == BND_PART_OPEN &&
	at_hi_boundary(p, 2) && ppsc->domain.bnd_part_hi[2] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int iy = 0;
	int iz = ppatch->ldims[2] - 1;
	F_av(m, 0,iy,iz) = (F(m, 0,iy+1,iz-1) +
				     F(m, 0,iy+1,iz-2) +
				     F(m, 0,iy+1,iz-3) +
				     F(m, 0,iy+2,iz-1) +
				     F(m, 0,iy+2,iz-2) +
				     F(m, 0,iy+2,iz-3) +
				     F(m, 0,iy+3,iz-1) +
				     F(m, 0,iy+3,iz-2) +
				     F(m, 0,iy+3,iz-3)) / 9.;
      }
    }

    if (at_hi_boundary(p, 1) && ppsc->domain.bnd_part_hi[1] == BND_PART_OPEN &&
	at_lo_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int iy = ppatch->ldims[1] - 1;
	int iz = 0;
	F_av(m, 0,iy,iz) = (F(m, 0,iy-1,iz+1) +
				     F(m, 0,iy-1,iz+2) +
				     F(m, 0,iy-1,iz+3) +
				     F(m, 0,iy-2,iz+1) +
				     F(m, 0,iy-2,iz+2) +
				     F(m, 0,iy-2,iz+3) +
				     F(m, 0,iy-3,iz+1) +
				     F(m, 0,iy-3,iz+2) +
				     F(m, 0,iy-3,iz+3)) / 9.;
      }
    }

    if (at_hi_boundary(p, 1) && ppsc->domain.bnd_part_hi[1] == BND_PART_OPEN &&
	at_hi_boundary(p, 2) && ppsc->domain.bnd_part_hi[2] == BND_PART_OPEN) {
      for (int m = 0; m < mflds->nr_fields; m++) {
	int iy = ppatch->ldims[1] - 1;
	int iz = ppatch->ldims[2] - 1;
	F_av(m, 0,iy,iz) = (F(m, 0,iy-1,iz-1) +
				     F(m, 0,iy-1,iz-2) +
				     F(m, 0,iy-1,iz-3) +
				     F(m, 0,iy-2,iz-1) +
				     F(m, 0,iy-2,iz-2) +
				     F(m, 0,iy-2,iz-3) +
				     F(m, 0,iy-3,iz-1) +
				     F(m, 0,iy-3,iz-2) +
				     F(m, 0,iy-3,iz-3)) / 9.;
      }
    }
  }
}

static void
average_in_time(struct psc_bnd_particles *bnd,
		struct psc_mfields *mflds_av, struct psc_mfields *mflds_last)
{
  mfields_t mf_av(mflds_av), mf_last(mflds_last);
  const double R = bnd->time_relax;
  for (int p = 0; p < ppsc->nr_patches; p++) {
    struct psc_patch *ppatch = &ppsc->patch[p];
    Fields F_av(mf_av[p]), F_last(mf_last[p]);
    
    for (int m = 0; m < mflds_av->nr_fields; m++) {
      if (bnd->first_time) {
	for (int iz = 0; iz < ppatch->ldims[2]; iz++) {
	  for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	    F_last(m, 0,iy,iz) = F_av(m, 0,iy,iz);
	  }
	}
      } else {
	// at lower and upper z bnd only FIXME
	for (int iz = 0; iz < ppatch->ldims[2]; iz += ppatch->ldims[2] + 1) {
	  for (int iy = 0; iy < ppatch->ldims[1]; iy++) {
	    F_av(m, 0,iy,iz) = 
	      R * F_av(m, 0,iy,iz)
	      + (1. - R) * F_last(m, 0,iy,iz);
	    F_last(m, 0,iy,iz) = F_av(m, 0,iy,iz);
	  }
	}
	for (int iz = 1; iz < ppatch->ldims[2] - 1; iz++) {
	  for (int iy = 0; iy < ppatch->ldims[1]; iy += ppatch->ldims[1] + 1) {
	    F_av(m, 0,iy,iz) = 
	      R * F_av(m, 0,iy,iz)
	      + (1. - R) * F_last(m, 0,iy,iz);
	    F_last(m, 0,iy,iz) = F_av(m, 0,iy,iz);
	  }
	}
      }
    }
  }
}

// ----------------------------------------------------------------------
// psc_bnd_particles_open_setup

static void _mrc_unused
psc_bnd_particles_open_setup(struct psc_bnd_particles *bnd)
{
  // open b.c. setup (FIXME: only when necessary)

  bnd->first_time = true;

  // FIXME: not necessary (?)
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  srandom(rank);

  bnd->flds_bnd = psc_bnd_create(psc_bnd_particles_comm(bnd));
  psc_bnd_set_type(bnd->flds_bnd, "c");
  psc_bnd_set_psc(bnd->flds_bnd, ppsc);
  psc_bnd_setup(bnd->flds_bnd);

  bnd->item_nvt = psc_output_fields_item_create(psc_bnd_particles_comm(bnd));
  //psc_output_fields_item_set_type(bnd->item_nvt, "nvp_1st_double"); FIXME
  //psc_output_fields_item_set_psc_bnd(bnd->item_nvt, bnd->flds_bnd);
  //psc_output_fields_item_setup(bnd->item_nvt);

  bnd->mflds_nvt_av = psc_output_fields_item_create_mfields(bnd->item_nvt);
  bnd->mflds_nvt_last = psc_output_fields_item_create_mfields(bnd->item_nvt);
  // A slight FIXME: We're not tracking partial particles at corners separately per direction,
  // which could lead so systematic errors, though it seems unlikely.
  bnd->mflds_n_in = psc_mfields_create(psc_bnd_particles_comm(bnd));
  psc_mfields_set_type(bnd->mflds_n_in, "c");
  psc_mfields_set_param_int(bnd->mflds_n_in, "nr_fields", ppsc->nr_kinds);
  psc_mfields_set_param_int3(bnd->mflds_n_in, "ibn", ppsc->ibn);
  bnd->mflds_n_in->grid = &ppsc->grid();
  psc_mfields_setup(bnd->mflds_n_in);
}

// ----------------------------------------------------------------------
// psc_bnd_particles_open_unsetup

static void _mrc_unused
psc_bnd_particles_open_unsetup(struct psc_bnd_particles *bnd)
{
  psc_mfields_destroy(bnd->mflds_nvt_av);
  psc_mfields_destroy(bnd->mflds_nvt_last);
  psc_mfields_destroy(bnd->mflds_n_in);
  psc_output_fields_item_destroy(bnd->item_nvt);
  psc_bnd_destroy(bnd->flds_bnd);
}

static void _mrc_unused
psc_bnd_particles_sub_open_calc_moments(struct psc_bnd_particles *bnd,
					struct psc_mparticles *mprts_base)
{
  static int pr_A, pr_B, pr_C;
  if (!pr_A) {
    pr_A = prof_register("open_calc_moments", 1., 0, 0);
    pr_B = prof_register("open_item_run", 1., 0, 0);
    pr_C = prof_register("open_average", 1., 0, 0);
  }

  prof_start(pr_A);

  static struct mrc_io *io;
  if (!io) {
    io = mrc_io_create(psc_comm(ppsc));
    mrc_io_set_type(io, "xdmf_collective");
    mrc_io_set_param_string(io, "basename", "open");
    mrc_io_set_from_options(io);
    mrc_io_setup(io);
  }

  if (ppsc->timestep % debug_every_step == 0) {
    mrc_io_open(io, "w", ppsc->timestep, ppsc->timestep * ppsc->dt);
  }

  mparticles_t mprts = mprts_base->get_as<mparticles_t>();
  struct psc_mfields *mflds_nvt = psc_output_fields_item_create_mfields(bnd->item_nvt);

  prof_start(pr_B);
  psc_output_fields_item_run(bnd->item_nvt, ppsc->flds, mprts.mprts(), mflds_nvt);
  prof_stop(pr_B);

  //debug_dump(io, mflds_nvt);

  prof_start(pr_C);
  average_9_point(bnd->mflds_nvt_av, mflds_nvt);
  debug_dump(io, bnd->mflds_nvt_av);

  average_in_time(bnd, bnd->mflds_nvt_av, bnd->mflds_nvt_last);
  prof_stop(pr_C);

  psc_mfields_destroy(mflds_nvt);

  bnd->first_time = false;

  mprts.put_as(mprts_base, MP_DONT_COPY);

  if (ppsc->timestep % debug_every_step == 0) {
    mrc_io_close(io);
  }
  prof_stop(pr_A);
}

enum {
  XX,
  YY,
  ZZ,
  XY,
  YZ,
  ZX,
};

#ifndef NO_OPEN_BC

static void
calc_W(double W[6], double vv[6])
{
  double determ =
    (vv[XX]*vv[YY]*vv[ZZ] + vv[XY]*vv[YZ]*vv[ZX] + vv[ZX]*vv[XY]*vv[YZ] -
     vv[ZX]*vv[YY]*vv[ZX] - vv[XY]*vv[XY]*vv[ZZ] - vv[XX]*vv[YZ]*vv[YZ]);
  W[XX] = .5 * (vv[YY]*vv[ZZ] - vv[YZ]*vv[YZ]) / determ;
  W[YY] = .5 * (vv[XX]*vv[ZZ] - vv[ZX]*vv[ZX]) / determ;
  W[ZZ] = .5 * (vv[XX]*vv[YY] - vv[XY]*vv[XY]) / determ;
  W[XY] = .5 * (vv[XY]*vv[ZZ] - vv[ZX]*vv[YZ]) / determ;
  W[ZX] = .5 * (vv[XY]*vv[YZ] - vv[ZX]*vv[YY]) / determ;
  W[YZ] = .5 * (vv[XX]*vv[YZ] - vv[XY]*vv[ZX]) / determ;
}

enum {
  NVT_N,
  NVT_VX,
  NVT_VY,
  NVT_VZ,
  NVT_VXVX,
  NVT_VYVY,
  NVT_VZVZ,
  NVT_VXVY,
  NVT_VYVZ,
  NVT_VZVX,
};

static inline double
inject_particles(int p, struct psc_mparticles *mprts, fields_t flds, 
		 fields_t flds_nvt_av, int ix, int iy, int iz,
		 double ninjo, int kind, double pos[3], double dir,
		 int X, int Y, int Z)
{
  Fields F(flds), F_nvt_av(flds_nvt_av);
  mparticles_t mp(mprts);
  mparticles_t::patch_t& prts = mp[p];

  double n     =         F_nvt_av(10*kind + NVT_N     , ix,iy,iz);
  double v[3]  = {       F_nvt_av(10*kind + NVT_VX, ix,iy,iz),
		         F_nvt_av(10*kind + NVT_VY, ix,iy,iz),
			 F_nvt_av(10*kind + NVT_VZ, ix,iy,iz), };
  double vv[6] = {       F_nvt_av(10*kind + NVT_VXVX + X, ix,iy,iz),
			 F_nvt_av(10*kind + NVT_VXVX + Y, ix,iy,iz),
			 F_nvt_av(10*kind + NVT_VXVX + Z, ix,iy,iz),
			 F_nvt_av(10*kind + NVT_VXVY + X, ix,iy,iz),
		   dir * F_nvt_av(10*kind + NVT_VXVY + Y, ix,iy,iz),
		   dir * F_nvt_av(10*kind + NVT_VXVY + Z, ix,iy,iz), };
  /* n = 1.; */
  /* v[0] = 0.; v[1] = 0.; v[2] = .1; */

  v[Z] *= dir;

  const Grid_t& grid = ppsc->grid;
  
  double W[6];
  calc_W(W, vv);

  double c=1.0;
  double vsz = sqrt(2. * vv[ZZ]);
  double gs0 = exp(-sqr(v[Z]) / sqr(vsz)) - exp(-sqr(c - v[Z]) / sqr(vsz))
    +sqrt(M_PI) * v[Z] / vsz * (erf((c - v[Z]) / vsz) + erf(v[Z] / vsz));
  double ninjn = ninjo + ppsc->dt * gs0 * n
    * vsz / sqrt(M_PI) / 2. / grid.dx[Z] / ppsc->coeff.cori;

  int ninjc = (int) ninjn;
  /* mprintf("n ele %g ninjo %g ninjn %g ninjon %g ninjc %d\n", n, */
  /* 	  ninjo, ninjn, ninjn - ninjc, ninjc); */
  ninjo = ninjn - ninjc;

  auto& grid = prts.grid();
  int nvdx = 1000;
  double dvz = c / ((double) nvdx);
	  
  if (ninjc != 0) {
    double vzdin = 0.;
    double  fin[nvdx];
    for (int jj = 0; jj < nvdx; jj++){
      vzdin += dvz;
      fin[jj] = (exp(-sqr(v[Z] / vsz)) - exp(-sqr(vzdin - v[Z]) / sqr(vsz)) + 
		 sqrt(M_PI) * v[Z] / vsz * (erf((vzdin - v[Z]) / vsz) + erf(v[Z] / vsz))) / gs0;
    }
    for (int n = 0; n < ninjc; n++) {
      unsigned int i = prts.size();
      mp[p].resize(i + 1);
      particle_t *prt = &prts[i]; 
      prt->kind_ = kind;
      prt->qni_wni = grid.kinds[kind].q;

      real_t *pxi = &prt->pxi;
      real_t *xi  = &prt->xi;

      int nnm = 0;
      do {
	nnm++;
	// FIXME, shouldn't have to loop here
	do {
	  double sr = random1();
      
	  pxi[Z] = 0.;
	  for (int k = 0; k < nvdx - 1; k++) {
	    if (sr >= fin[k] && sr < fin[k+1]) {
	      pxi[Z] = dvz * (k + 1) + (sr - fin[k]) * dvz / (fin[k+1] - fin[k]);
	      break;
	    }
	  }
	} while (pxi[Z] == 0);
        
	double sr = random1();
	double yya = 0.;
	double yy0;
	int icount = 0;
	do {
	  icount++;
	  yy0 = yya;
	  yya = yy0 - (erf(yy0) - (2.*sr - 1.)) / (2./sqrt(M_PI) * exp(-sqr(yy0)));
	} while (fabs(yya - yy0) > 1.e-15 && icount != 100);
	pxi[X] = v[X] + yya * sqrt(W[YY] / (W[XX] * W[YY] - sqr(W[XY])))
	  + (pxi[Z] - v[Z]) * vv[ZX] / vv[ZZ];
   
	sr = random1();
	yya = 0.0;
	icount = 0;
	do {
	  icount++;
	  yy0 = yya;
	  yya = yy0 - (erf(yy0) - (2.*sr - 1.)) / (2./sqrt(M_PI) * exp(-sqr(yy0)));
	} while (fabs(yya - yy0) > 1.e-15 && icount != 100);
	pxi[Y] = v[Y] + 1. / W[YY] * (yya * sqrt(W[YY])
				     - (pxi[Z] - v[Z]) * W[YZ] - (pxi[X] - v[X]) * W[XY]);
		
	if (nnm > 100) { assert(0); break; }
      } while (sqr(pxi[X]) + sqr(pxi[Y]) + sqr(pxi[Z]) > 1.);

      double xr = random1();
      xi[X] = pos[X] + xr * grid.dx[X];
      double yr = random1();
      xi[Y] = pos[Y] + yr * grid.dx[Y];
      double zr = random1();
      double dz = dir * zr * pxi[Z] * ppsc->dt;
      xi[Z] = pos[Z] + dz;

      double Jz = prt->qni_wni * dz * ppsc->coeff.cori / ppsc->dt;
      if (Z == 2) {
	F(JXI + Z, ix,iy  ,iz) += (1 - yr) * Jz;
	F(JXI + Z, ix,iy+1,iz) += (    yr) * Jz;
      } else if (Z == 1) {
	F(JXI + Z, ix,iy,iz  ) += (1 - xr) * Jz;
	F(JXI + Z, ix,iy,iz+1) += (    xr) * Jz;
      } else {
	assert(0);
      }
      double gamma = 1. / sqrt(1. - (sqr(pxi[X]) + sqr(pxi[Y]) + sqr(pxi[Z])));
      if (sqr(pxi[X]) + sqr(pxi[Y]) + sqr(pxi[Z]) > 1.) {
	gamma = 1.;
	assert(0);
      }
      pxi[X] *= gamma;
      pxi[Y] *= gamma;
      pxi[Z] *= dir * gamma;
    }
  }

  return ninjo;
}

static double
inject_particles_y(int p, struct psc_mparticles *mprts, fields_t flds, 
		   fields_t flds_nvt_av, int ix, int iy, int iz,
		   double ninjo, int kind, double pos[3], double dir)
{
  return inject_particles(p, mprts, flds, flds_nvt_av, ix, iy, iz, ninjo, kind, pos, dir,
			  2, 0, 1);
}

static double
inject_particles_z(int p, struct psc_mparticles *mprts, fields_t flds, 
		   fields_t flds_nvt_av, int ix, int iy, int iz,
		   double ninjo, int kind, double pos[3], double dir)
{
  return inject_particles(p, mprts, flds, flds_nvt_av, ix, iy, iz, ninjo, kind, pos, dir,
			  0, 1, 2);
}

#endif

static void _mrc_unused
psc_bnd_particles_open_boundary(struct psc_bnd_particles *bnd, struct psc_mparticles *mprts,
				struct psc_mfields *mflds)
{
#ifndef NO_OPEN_BC
  static int pr_A;
  if (!pr_A) {
    pr_A = prof_register("open_boundary", 1., 0, 0);
  }

  prof_start(pr_A);
  const Grid_t& grid = ppsc->grid;
  int nr_kinds = ppsc->nr_kinds;
  mfields_t mf_nvt_av(bnd->mflds_nvt_av), mf_n_in(bnd->mflds_n_in), mf(mflds);

  for (int p = 0; p < ppsc->nr_patches; p++) {
    fields_t flds_nvt_av = mf_nvt_av[p];
    Fields F_n_in(mf_n_in[p]);
    fields_t flds = mf[p];

    for (int m = 0; m < nr_kinds; m++) {
      // inject at y lo
      if (at_lo_boundary(p, 1) && ppsc->domain.bnd_part_lo[1] == BND_PART_OPEN) {
	int iy = 0;
	for (int iz = 0; iz < grid.ldims[2]; iz++) {
	  double ninjo = F_n_in(m, 0,iy,iz);
	  double pos[3] = { 0., 0., iz * grid.dx[2], };
	  F_n_in(m, 0,iy,iz) =
	    inject_particles_y(p, mprts, flds, flds_nvt_av, 0,iy,iz, ninjo, m, pos, +1.);
	}
      }
      // inject at y hi
      if (at_hi_boundary(p, 1) && ppsc->domain.bnd_part_hi[1] == BND_PART_OPEN) {
	int iy = grid.ldims[1] - 1;
	for (int iz = 0; iz < grid.ldims[2]; iz++) {
	  double ninjo = F_n_in(m, 0,iy,iz);
	  double pos[3] = { 0., (iy + 1) * (1-1e-6) * grid.dx[1], iz * grid.dx[2] };
	  F_n_in(m, 0,iy,iz) =
	    inject_particles_y(p, mprts, flds, flds_nvt_av, 0,iy,iz, ninjo, m, pos, -1.);
	}
      }
      // inject at z lo
      if (at_lo_boundary(p, 2) && ppsc->domain.bnd_part_lo[2] == BND_PART_OPEN) {
	int iz = 0;
	for (int iy = 0; iy < grid.ldims[1]; iy++) {
	  double ninjo = F_n_in(m, 0,iy,iz);
	  double pos[3] = { 0., iy * grid.dx[1], 0. };
	  F_n_in(m, 0,iy,iz) =
	    inject_particles_z(p, mprts, flds, flds_nvt_av, 0,iy,iz, ninjo, m, pos, +1.);
	}
      }
      // inject at z hi
      if (at_hi_boundary(p, 2) && ppsc->domain.bnd_part_hi[2] == BND_PART_OPEN) {
	int iz = grid.ldims[2] - 1;
	for (int iy = 0; iy < grid.ldims[1]; iy++) {
	  double ninjo = F_n_in(m, 0,iy,iz);
	  double pos[3] = { 0., iy * grid.dx[1], (iz + 1) * (1-1e-6) * grid.dx[2] };
	  F_n_in(m, 0,iy,iz) =
	    inject_particles_z(p, mprts, flds, flds_nvt_av, 0,iy,iz, ninjo, m, pos, -1.);
	}
      }
    }
  }
  prof_stop(pr_A);
#endif
}


