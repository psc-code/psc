/*

 Decaying turbulence. The simulation is initialized with a set of large-scale 
 modes with phases adjusted to minimize the variance of |B|. See Roberts, PRL, 2012 for details on how to compute the initial field.

 */

#include <psc.h>
#include <psc_push_particles.h>
#include <psc_push_fields.h>
#include <psc_sort.h>
#include <psc_balance.h>
#include <psc_particles_as_single.h>

#include <mrc_params.h>
#include <mrc_profile.h>

#include <math.h>
#include <time.h>

enum {
  T_ELECTRON,
  T_ION,
  NR_T_KINDS,
};



struct psc_turb {
  // parameters
  double beta;      // plasma beta
  double mi_me;     // mass ratio
  double wpe_wce;   
  double dB_B0;     // RMS amplitude of the perturbation
  double dB_2_dB_1; // amplitude of modes with polarization 2 wrt that of polarization 1

  // parameters of the modes
  int nmx,nmy,nmz,nmodes;      // number of modes in x,y,z and the total number
  double *phases;              // mode ampltidue
  double *amp;                 // |k| for each mode
  double *absk;                // mode amplitude
  double *dbx, *dby, *dbz;     // unit vector describing polarization for each mode
  double *kx, *ky, *kz;        // kx,ky,kz for each mode

  // calculated from the above
  double B0; 
  double T;
  double dB;
};

#define to_psc_turb(psc) mrc_to_subobj(psc, struct psc_turb)

#define VAR(x) (void *)offsetof(struct psc_turb, x)

static struct param psc_turb_descr[] = {
  { "beta"     , VAR(beta)         , PARAM_DOUBLE(.5)            },
  { "dB_B0"     , VAR(dB_B0)       , PARAM_DOUBLE(.5)            },
  { "dB_2_dB_1" , VAR(dB_2_dB_1)   , PARAM_DOUBLE(.1)            },
  { "mi_me"    , VAR(mi_me)        , PARAM_DOUBLE(100.)           },
  { "wpe_wce"  , VAR(wpe_wce)      , PARAM_DOUBLE(2.)            },
  { "nmx"     , VAR(nmx)           , PARAM_INT(0)                },
  { "nmy"     , VAR(nmy)           , PARAM_INT(5)                },
  { "nmz"     , VAR(nmz)           , PARAM_INT(5)                },
  {},
};
#undef VAR

// auxiliary grid to compute the initial field 
struct phase_grid {
  double  dx, dy, dz;  // spacing
  int nx,ny,nz;        // # of points
};


// ----------------------------------------------------------------------
// psc_turb_create

static void
psc_turb_create(struct psc *psc)
{
  psc_default_dimensionless(psc);

  psc->prm.nmax = 62415*3;
  psc->prm.nicell = 500;
  psc->prm.cfl = 0.98;

  struct psc_kind kinds[NR_T_KINDS] = {
    [T_ELECTRON] = { .name = strdup("e"), .q = -1., .m = 1, },
    [T_ION]      = { .name = strdup("i"), .q =  1.,         },
  };
  psc_set_kinds(psc, NR_T_KINDS, kinds);


  psc->domain.length[0] = 0.176776695; 
  psc->domain.length[1] = 314.1592654; 
  psc->domain.length[2] = 314.1592654;

  psc->domain.gdims[0] = 1;
  psc->domain.gdims[1] = 1728;
  psc->domain.gdims[2] = 1728;


  psc->domain.bnd_fld_lo[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[0] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[1] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_lo[2] = BND_FLD_PERIODIC;
  psc->domain.bnd_fld_hi[2] = BND_FLD_PERIODIC;

  psc->domain.bnd_part_lo[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[0] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[1] = BND_PART_PERIODIC;
  psc->domain.bnd_part_lo[2] = BND_PART_PERIODIC;
  psc->domain.bnd_part_hi[2] = BND_PART_PERIODIC;


  // FIXME: can only use 1st order pushers with current conducting wall b.c.
  // try second order?
  // psc_push_particles_set_type(psc->push_particles, "1vb");
}

// ----------------------------------------------------------------------
// psc_turb_setup
//
// the parameters are now set, calculate quantities to initialize fields,
// particles

#define swap(a,b) do {double tmp=(a); (a)=(b); (b)=tmp;} while(0)

// computes magnetic field at a goven point (x,y,z) for a set of modes
// 0<=p1<nmodes_1  and 0<=p2<nmodes_2, where p1 and p2 are indices into
// polarization 1 & 2 repsectively
static void
get_b(struct psc_turb *turb, int nmodes_1, int nmodes_2,double x, double y, double z, double *bx_, double *by_, double *bz_)
{
  int l;
  *bx_ = 0.0;
  *by_ = 0.0;
  *bz_ = 0.0;
  int kx,ky,kz;
  double cosfactor;
  int nmodes = turb->nmodes;

  // polarization 1
  for (l=1; l<nmodes_1; l++) {
    kx = turb->kx[l];
    ky = turb->ky[l];
    kz = turb->kz[l];
    cosfactor = cos( (2.0*M_PI*kx)*x + (2.0*M_PI*ky)*y + (2.0*M_PI*kz)*z + turb->phases[l] );
    *bx_ += turb->dbx[l]*turb->amp[l]*cosfactor;
    *by_ += turb->dby[l]*turb->amp[l]*cosfactor;
    *bz_ += turb->dbz[l]*turb->amp[l]*cosfactor;	  
  }
  // polarization 2
  for (l=1; l<nmodes_2; l++) {
    kx = turb->kx[l];
    ky = turb->ky[l];
    kz = turb->kz[l];
    cosfactor = cos( (2.0*M_PI*kx)*x + (2.0*M_PI*ky)*y + (2.0*M_PI*kz)*z + turb->phases[l+nmodes] );
    *bx_ += turb->dbx[l+nmodes]*turb->amp[l+nmodes]*cosfactor;
    *by_ += turb->dby[l+nmodes]*turb->amp[l+nmodes]*cosfactor;
    *bz_ += turb->dbz[l+nmodes]*turb->amp[l+nmodes]*cosfactor;	  
  }  
  // add the guide field
  *bx_ = *bx_ + turb->B0;
}



/*
 computes integral of |B| over the auxiliary grid. This is the quantity we want to maximize
*/
static void
get_absk_int(struct psc_turb *turb, int nmodes_1, int nmodes_2, struct phase_grid *pg, double *absk_int)
{
  int i,j,k;
  double absk,bx,by,bz;

  *absk_int = 0.0;
  for(i=0;i<pg->nx;i++)
    for(j=0;j<pg->ny;j++)
      for(k=0;k<pg->nz;k++)
	{
	  get_b(turb, nmodes_1, nmodes_2,i*pg->dx,j*pg->dy,k*pg->dz,&bx,&by,&bz);
	  absk = sqrt(bx*bx + by*by + bz*bz);
	  *absk_int += absk;
	}  
}

// initialize the parameters 
static void
psc_turb_setup(struct psc *psc)
{
  struct psc_turb *turb = to_psc_turb(psc);

  double B0 = 1.0 / (turb->wpe_wce);
  double T = turb->beta * sqr(B0) / 2.;

  turb->B0 = B0;
  turb->T = T;
  turb->dB = turb->dB_B0*B0;

  int l,m,n,idx,pidx,pidx_max,nmodes_1,nmodes_2;
  double dB2;
  double kx,ky,kz,e1_x,e1_y,e1_z,e2_x,e2_y,e2_z,absk,absk_int,absk_int_max;
  const int ntries = 10;
  int ml[3],mh[3],nmodes;

  // count the modes
  nmodes=0;
  ml[2] = 0; 
  mh[2] = turb->nmz;
  mh[1] = turb->nmy;
  mh[0] = turb->nmx;
  mpi_printf(MPI_COMM_WORLD,"----------------- will initialize these modes\n");
  for(n=ml[2]; n<=mh[2]; n++)  {
    if (n==0) ml[1] = 0; else ml[1] = -turb->nmy;
    for(m=ml[1]; m<=mh[1]; m++) {
      if ((n==0)&&(m==0)) ml[0] = 0; else ml[0] = -turb->nmx;
      for(l=ml[0]; l<=mh[0]; l++) {
        mpi_printf(MPI_COMM_WORLD,"%3i %3i %3i\n",l,m,n);
        nmodes ++;
      }
    }
  }
  mpi_printf(MPI_COMM_WORLD,"--------------------------------------------\n");
  mpi_printf(MPI_COMM_WORLD,"Total modes: %i\n", nmodes);

  // allocate the arrays
  struct phase_grid *pg = (struct phase_grid*) malloc(sizeof( *pg));
  turb->amp = (double*) malloc(sizeof(double)*nmodes*2);
  turb->absk = (double*) malloc(sizeof(double)*nmodes);
  turb->dbx = (double*) malloc(sizeof(double)*nmodes*2);
  turb->dby = (double*) malloc(sizeof(double)*nmodes*2);
  turb->dbz = (double*) malloc(sizeof(double)*nmodes*2);
  turb->kx = (double*) malloc(sizeof(double)*nmodes);
  turb->ky = (double*) malloc(sizeof(double)*nmodes);
  turb->kz = (double*) malloc(sizeof(double)*nmodes);
  turb->phases = (double*) malloc(sizeof(double)*nmodes*2);

  turb->nmodes = nmodes;

  mpi_printf(MPI_COMM_WORLD, "psc/turb: initializaing the modes\n");

  srand (10);  // synchronize rng's
  
  idx = -1;

  dB2 = 0.0;

  // initialize phases and amplitudes
  for(n=ml[2]; n<=mh[2]; n++)  {
    if (n==0) ml[1] = 0; else ml[1] = -turb->nmy;
    for(m=ml[1]; m<=mh[1]; m++) {
      if ((n==0)&&(m==0)) ml[0] = 0; else ml[0] = -turb->nmx;
      for(l=ml[0]; l<=mh[0]; l++)
	{
	  
	  kx = (2.0*M_PI/psc->domain.length[0])*l;
	  ky = (2.0*M_PI/psc->domain.length[1])*m;
	  kz = (2.0*M_PI/psc->domain.length[2])*n;
	  absk = sqrt(kx*kx+ky*ky+kz*kz);

	  idx++;
	  turb->absk[idx] = absk;
	  turb->kx[idx] = l;
	  turb->ky[idx] = m;
	  turb->kz[idx] = n;
	  turb->amp[idx] = 0.0;
	  turb->amp[idx+nmodes] = 0.0;

	  if  ( (n==0) && (m==0) && (l==0) ) continue;

	  // find two unit vectors orthogonal to k an to each other
	  if ((m==0) && (n==0)) // purely parallel propagation
	    {
	      double rnd_phi = 2*M_PI*( rand() / (float) RAND_MAX );
	      e1_x = 0.0;
	      e1_y = cos(rnd_phi);
	      e1_z = sin(rnd_phi);

	      e2_x = 0.0;
	      e2_y = -sin(rnd_phi);
	      e2_z =  cos(rnd_phi);

	    }  
	  else
	    {
	      // choose e1 = (ex x k)

	      double e1 = sqrt(ky*ky+kz*kz); 

	      e1_x =  0;
	      e1_y = -kz/e1;
	      e1_z =  ky/e1;


	      // chose e2 = (e1 x k)
	      double e2 = sqrt( (e1*e1+kx*kx)*e1*e1 );
	      e2_x = -(ky*ky + kz*kz)/e2;
	      e2_y = ky*kx/e2;
	      e2_z = kz*kx/e2; 


	    }  
	    
	  turb->amp[idx] = 1/pow(absk,1.0);
	  turb->amp[idx+nmodes] = turb->dB_2_dB_1*turb->amp[idx];

	  dB2 = dB2 + 0.5*turb->amp[idx]*turb->amp[idx] + 0.5*turb->amp[idx+nmodes]*turb->amp[idx+nmodes]; // 0.5 comes from <cos^2>
	  turb->phases[idx]        = 2.0*M_PI*( rand() / (float) RAND_MAX );  // polarization 1
	  turb->phases[idx+nmodes] = 2.0*M_PI*( rand() / (float) RAND_MAX );  // polarization 2

	  // save the polarization information
	  turb->dbx[idx] = e1_x;
	  turb->dby[idx] = e1_y;
	  turb->dbz[idx] = e1_z;

	  turb->dbx[idx+nmodes] = e2_x;
	  turb->dby[idx+nmodes] = e2_y;
	  turb->dbz[idx+nmodes] = e2_z;	 
	}
    }
  }
  // normalize the amplitudes so that we have the prescribed RMS amplitude
  for (l=0; l<nmodes*2; l++) turb->amp[l] *= turb->dB/sqrt(dB2);  

  // sort the arrays by absk. the phases atray does not need to be sorted since the phase are random at this point

  for (l=1; l<nmodes; l++) {
    idx = l;
    while ( (turb->absk[idx] < turb->absk[idx-1]) && (idx>0) ) {
      
      // absk
      swap(turb->absk[idx-1],turb->absk[idx]);      
      // amp
      swap(turb->amp[idx-1],turb->amp[idx]);
      swap(turb->amp[idx-1+nmodes],turb->amp[idx+nmodes]);
      // dbx
      swap(turb->dbx[idx-1],turb->dbx[idx]);
      swap(turb->dbx[idx-1+nmodes],turb->dbx[idx+nmodes]);
      // dby
      swap(turb->dby[idx-1],turb->dby[idx]);
      swap(turb->dby[idx-1+nmodes],turb->dby[idx+nmodes]);
      // dbz
      swap(turb->dbz[idx-1],turb->dbz[idx]);
      swap(turb->dbz[idx-1+nmodes],turb->dbz[idx+nmodes]);
      // kx
      swap(turb->kx[idx-1],turb->kx[idx]);
      // ky
      swap(turb->ky[idx-1],turb->ky[idx]);
      // kz
      swap(turb->kz[idx-1],turb->kz[idx]);

      idx -= 1;
    }
  }

  // parameters of the auxiliary grid
  int Nf = 4 ;    // extra factor compared to the Nyquist criterion
  
  if (turb->nmx > 0) {
    pg->dx = 1.0/(2*turb->nmx)/Nf;
    pg->nx = (int) (1.0/pg->dx);
  } else {
    pg->nx = 1;
    pg->dx = 0;
  }

  if (turb->nmy > 0) {
    pg->dy = 1.0/(2*turb->nmy)/Nf;
    pg->ny = (int) (1.0/pg->dy);
  } else {
    pg->ny = 1;
    pg->dy = 0;
  }

  if (turb->nmz > 0) {
    pg->dz = 1.0/(2*turb->nmz)/Nf;
    pg->nz = (int) (1.0/pg->dz);
  } else {
    pg->nz = 1;
    pg->dz = 0;
  }


  // fix the first polarizaiton and optmize the phases of the second
  absk_int_max = 0.0;
  for (nmodes_2=2; nmodes_2<=nmodes; nmodes_2++) 
    {
      // try phases  
      pidx_max = -1;    
      for (pidx=0; pidx<ntries; pidx++)
	{
	  turb->phases[nmodes_2-1+nmodes]        = (2.0*M_PI/ntries)*pidx;
	  get_absk_int(turb,nmodes,nmodes_2,pg,&absk_int);
	  if (absk_int > absk_int_max)
	    {
	      absk_int_max = absk_int;
	      pidx_max = pidx;
	    }
	}
      assert(pidx_max >= 0);
      turb->phases[nmodes_2-1+nmodes]        = (2.0*M_PI/ntries)*pidx_max;
    }

  // fix the phases of the second ploarization and optimize the first
  absk_int_max = 0.0;
  for (nmodes_1=2; nmodes_1<=nmodes; nmodes_1++) 
    {
      // try phases   
  pidx_max = -1;   
      for (pidx=0; pidx<ntries; pidx++)
	{
	  turb->phases[nmodes_1-1]        = (2.0*M_PI/ntries)*pidx;
	  get_absk_int(turb,nmodes_1,nmodes,pg,&absk_int);
	  if (absk_int > absk_int_max)
	    {
	      absk_int_max = absk_int;
	      pidx_max = pidx;
	    }
	}
      turb->phases[nmodes_1-1]        = (2.0*M_PI/ntries)*pidx_max;
    }
  // ---- end compute modes


  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank==0) {
    // write modes to file
    FILE *f = fopen("modes.dat","w");
    fprintf(f,"#  mode description. the columns are\n");
    fprintf(f,"#  1: mode number\n");
    fprintf(f,"#  2: kx (int) \n");
    fprintf(f,"#  3: ky (int) \n");
    fprintf(f,"#  4: kz (int) \n");
    fprintf(f,"#  5: mode amplitude\n");
    fprintf(f,"#  6: dbx (polarization 1)\n");
    fprintf(f,"#  7: dby (polarization 1)\n");
    fprintf(f,"#  8: dbz (polarization 1)\n");
    fprintf(f,"#  9: phases (polarization 1)\n");
    fprintf(f,"#  10: mode amplitude\n");
    fprintf(f,"#  11: dbx (polarization 2)\n");
    fprintf(f,"#  12: dby (polarization 2)\n");
    fprintf(f,"#  13: dbz (polarization 2)\n");
    fprintf(f,"#  14: phases (polarization 2)\n");
    for (l=1; l<nmodes; l++)
      fprintf(f,"%3i %3.0f %3.0f %3.0f %6.3E %6.3E %6.3E %6.3E %6.3E %6.3E %6.3E %6.3E %6.3E %6.3E\n",l,
	      turb->kx[l],
	      turb->ky[l],
	      turb->kz[l],
	      turb->amp[l],
	      turb->dbx[l],
	      turb->dby[l],
	      turb->dbz[l],
	      turb->phases[l],
	      turb->amp[l+nmodes],
	      turb->dbx[l+nmodes],
	      turb->dby[l+nmodes],
	      turb->dbz[l+nmodes],
	      turb->phases[l+nmodes]);
    fclose(f);
    
  }

  mpi_printf(MPI_COMM_WORLD, "psc/turb: finished initializaing the phases\n");

  // set particle kind parameters
  assert(psc->nr_kinds == NR_T_KINDS);
  psc->kinds[T_ELECTRON].T = turb->T;
  psc->kinds[T_ION].T = turb->T;
  psc->kinds[T_ION].m = turb->mi_me;;

  psc_setup_super(psc);
}


// ----------------------------------------------------------------------
// psc_turb_init_field

static double
psc_turb_init_field(struct psc *psc, double x[3], int m)
{
  struct psc_turb *turb = to_psc_turb(psc);
  double bx,by,bz,xhat,yhat,zhat;
  

  if ( (m==HX) || (m==HY) || (m==HZ) ) {
    int nmodes = turb->nmodes;
    xhat =   x[0]/psc->domain.length[0];
    yhat =   x[1]/psc->domain.length[1];
    zhat =   x[2]/psc->domain.length[2];

    get_b(turb,nmodes,nmodes,xhat,yhat,zhat,&bx,&by,&bz);
  }

  switch (m) {
  case HX: return bx;
  case HY: return by;
  case HZ: return bz;
  default: return 0.;
  }
}

// ----------------------------------------------------------------------
// psc_turb_init_npt

static void
psc_turb_init_npt(struct psc *psc, int kind, double x[3],
		struct psc_particle_npt *npt)
{
  //  struct psc_turb *turb = to_psc_turb(psc);

  npt->n = 1;

  switch (kind) {
  case T_ELECTRON: // electrons
    break;
  case T_ION:     // ions
    break;
  default:
    assert(0);
  }
}


static void
psc_turb_read(struct psc *psc, struct mrc_io *io)
{
  // do nothing -- but having this function is important so that
  // psc_kh_create() doesn't get called instead
  psc_read_super(psc, io);
}


// ======================================================================
// psc_turb_ops

struct psc_ops psc_turb_ops = {
  .name             = "turb",
  .size             = sizeof(struct psc_turb),
  .read             = psc_turb_read,
  .param_descr      = psc_turb_descr,
  .create           = psc_turb_create,
  .setup            = psc_turb_setup,
  .init_field       = psc_turb_init_field,
  .init_npt         = psc_turb_init_npt,
};

// ======================================================================
// main

int
main(int argc, char **argv)
{
  return psc_main(&argc, &argv, &psc_turb_ops);
}
