
#include <mrc_io.h>
#include <mrc_params.h>
#include <mrc_profile.h>

#include <assert.h>
#include <string.h>
#include <math.h>

struct sph_harm_params {
  int l;
  int m;
  int n;
  float c0;
  char *ic;
};

#define VAR(x) (void *)offsetof(struct sph_harm_params, x)

static struct param sph_harm_params_descr[] = {
  { "l"            , VAR(l)               , PARAM_INT(0)            },
  { "m"            , VAR(m)               , PARAM_INT(0)            },
  { "n"            , VAR(n)               , PARAM_INT(0)            },
  { "c0"           , VAR(c0)              , PARAM_FLOAT(0.)         },
  { "ic"           , VAR(ic)              , PARAM_STRING("one")     },
  {},
};

#undef VAR

static float
factorial(int n)
{
  int res = 1;
  for (int i = 1; i <=n; i++) {
    res *= i;
  }
  return res;
}

static float
P(int l, float z)
{
  switch (l) {
  case 1: return z;
  case 3: return 1./2.  * (5. * pow(z, 3.) - 3. * z);
  case 5: return 1./8.  * (63. * pow(z, 5.) - 70. * pow(z, 3.) + 15. * z);
  case 7: return 1./16. * (429. * pow(z, 7.) - 693. * pow(z, 5.) + 315. * pow(z, 3.) - 35. * z);
  case 9: return 1./128. * (12155 * pow(z, 9.) - 25740. * pow(z, 7.) + 18018. * pow(z, 5.) - 4620. * pow(z, 3.) + 315. * z);
  case 11: return 1./256. * (88179 * pow(z, 11.) - 230945 * pow(z, 9.) + 218790. * pow(z, 7.) - 90090. * pow(z, 5.) + 15015. * pow(z, 3.) - 693. * z);
  default: assert(0);
  }
}

static float
Y(int l, int m, float th, float phi)
{
  if (l == 0 && m == 0) {
    return sqrt(1. / (4.*M_PI));

  } else if (l == 1 && m == 0) {
    return sqrt(3. / (4.*M_PI)) * cos(th);
  } else if (l == 1 && m == -1) {
    return sqrt(3. / (8.*M_PI)) * sin(th) * sin(phi);
  } else if (l == 1 && m == 1) {
    return sqrt(3. / (8.*M_PI)) * sin(th) * cos(phi);

  } else if (l == 2 && m == 0) {
    return sqrt(5. / (16.*M_PI)) * (3. * pow(cos(th), 2.) - 1.);
  } else if (l == 2 && m == -1) {
    return sqrt(15. / (8.*M_PI)) * sin(th) * cos(th) * sin(phi);
  } else if (l == 2 && m ==  1) {
    return sqrt(15. / (8.*M_PI)) * sin(th) * cos(th) * cos(phi);
  } else if (l == 2 && m == -2) {
    return sqrt(15. / (32.*M_PI)) * pow(sin(th), 2.) * sin(2*phi);
  } else if (l == 2 && m == 2) {
    return sqrt(15. / (32.*M_PI)) * pow(sin(th), 2.) * cos(2*phi);

  } else if (l == 3 && m == 0) {
    return sqrt(7. / (16.*M_PI)) * (5. * pow(cos(th), 3.) - 3. * cos(th));
  } else if (l == 3 && m == 1) {
    return -sqrt(21. / (64.*M_PI)) * sin(th) * (5. * pow(cos(th), 2.) - 1.) * cos(phi);
  } else if (l == 3 && m == 2) {
    return sqrt(105. / (32.*M_PI)) * pow(sin(th), 2.) * cos(th) * cos(2*phi);
  } else if (l == 3 && m == 3) {
    return -sqrt(35. / (64.*M_PI)) * pow(sin(th), 3.) * cos(3*phi);

  } else if (l == 4 && m == 0) {
    return sqrt(9. / (256.*M_PI)) * (35. * pow(cos(th), 4.) - 30. * pow(cos(th), 2.) + 3.);

  } else {
    printf("not implemented: l = %d m = %d\n", l, m);
    assert(0);
  }
}

static float
laguerre(int n, int alpha, float x)
{
  if (n == 0) {
    return 1.;
  } else if (n == 1) {
    return -x + alpha + 1;
  } else if (n == 2) {
    return x*x - (alpha + 2) * x + (alpha + 2) * (alpha + 1) / 2.;
  } else {
    assert(0);
  }
  float L0, L1 = 0, res = 1.;
  for (int i = 1; i <= n; i++) {
    L0 = L1; L1 = res;
    res = ((2*i - 1 + alpha- x) * L1 - (i-1 + alpha) * L0) / i;
  }
  return res;
}

static void
ini_one(struct mrc_fld *f, struct sph_harm_params *par)
{
  struct mrc_crds *crds = mrc_domain_get_crds(f->_domain);
  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    float x = MRC_CRDX(crds, ix), y = MRC_CRDY(crds, iy), z = MRC_CRDZ(crds, iz);
    float r = sqrt(x*x + y*y + z*z), th = acos(z/r), phi = atan2(y/r, x/r);
    MRC_F3(f, 0, ix,iy,iz) = Y(par->l, par->m, th, phi);
  } mrc_fld_foreach_end;
}

static void
ini_semi(struct mrc_fld *f, struct sph_harm_params *par)
{
  struct mrc_crds *crds = mrc_domain_get_crds(f->_domain);
  for (int l = 0; l < (par->l + 1) / 2; l++) {
    float c = pow(-1., l) * factorial(2*l) / pow(pow(2., l) * factorial(l), 2.) * (4*l+3.)/(2*l+2.);
    printf("l = %d c = %g c' = %g\n", 2*l + 1, c, c * sqrt(2. / (4*l + 3)));
  }
  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    float x = MRC_CRDX(crds, ix), y = MRC_CRDY(crds, iy), z = MRC_CRDZ(crds, iz);
    float r = sqrt(x*x + y*y + z*z), th = acos(z/r);
    float val;
    if (r < 1.) {
      val = par->c0;
    } else {
      val = par->c0 / r;
    }
    for (int l = 0; l < (par->l + 1) / 2; l++) {
      float c = pow(-1., l) * factorial(2*l) / pow(pow(2., l) * factorial(l), 2.) *
	(4*l+3.)/(2*l+2.);
      if (r < 1.) {
	val += c * P(2*l + 1, cos(th)) * pow(r, 2*l + 1);
      } else {
	val += c * P(2*l + 1, cos(th)) * pow(r, -(2*l + 1) - 1);
      }
    }
    MRC_F3(f, 0, ix,iy,iz) = val;
  } mrc_fld_foreach_end;
}

static void
ini_hydrogen(struct mrc_fld *f, struct sph_harm_params *par)
{
  struct mrc_crds *crds = mrc_domain_get_crds(f->_domain);
  mrc_fld_foreach(f, ix,iy,iz, 0, 0) {
    float x = MRC_CRDX(crds, ix), y = MRC_CRDY(crds, iy), z = MRC_CRDZ(crds, iz);
    float r = sqrt(x*x + y*y + z*z), th = acos(z/r), phi = atan2(y/r, x/r);
    MRC_F3(f, 0, ix,iy,iz) = exp(-r) *
      laguerre(par->n - par->l - 1, 2*par->l + 1, r) * Y(par->l, par->m, th, phi);
  } mrc_fld_foreach_end;
}

static void
calc_grad(struct mrc_fld *f)
{
  struct mrc_crds *crds = mrc_domain_get_crds(f->_domain);
  mrc_fld_foreach(f, ix,iy,iz, -1, -1) {
    MRC_F3(f, 1, ix,iy,iz) = -
      (MRC_F3(f, 0, ix+1,iy,iz) - MRC_F3(f, 0, ix-1,iy,iz)) / 
      (MRC_CRDX(crds, ix+1) - MRC_CRDX(crds, ix-1));
    MRC_F3(f, 2, ix,iy,iz) = -
      (MRC_F3(f, 0, ix,iy+1,iz) - MRC_F3(f, 0, ix,iy-1,iz)) / 
      (MRC_CRDY(crds, iy+1) - MRC_CRDY(crds, iy-1));
    MRC_F3(f, 3, ix,iy,iz) = -
      (MRC_F3(f, 0, ix,iy,iz+1) - MRC_F3(f, 0, ix,iy,iz-1)) / 
      (MRC_CRDZ(crds, iz+1) - MRC_CRDZ(crds, iz-1));
  } mrc_fld_foreach_end;
}

// ----------------------------------------------------------------------
// main

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct sph_harm_params par;
  mrc_params_parse(&par, sph_harm_params_descr, "sph_harm", MPI_COMM_WORLD);
  mrc_params_print(&par, sph_harm_params_descr, "sph_harm", MPI_COMM_WORLD);

  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(domain, "simple");
  mrc_domain_set_from_options(domain);
  mrc_domain_view(domain);
  mrc_domain_setup(domain);

  struct mrc_fld *fld = mrc_domain_fld_create(domain, SW_0, "phi:ex:ey:ez");
  mrc_fld_setup(fld);
  if (strcmp(par.ic, "one") == 0) {
    ini_one(fld, &par);
  } else if (strcmp(par.ic, "semi") == 0) {
    ini_semi(fld, &par);
  } else if (strcmp(par.ic, "hydrogen") == 0) {
    ini_hydrogen(fld, &par);
  } else {
    assert(0);
  }
  calc_grad(fld);

  struct mrc_io *io = mrc_io_create(MPI_COMM_WORLD);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_fld_write(fld, io);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrc_fld_destroy(fld);
  mrc_domain_destroy(domain);
  MPI_Finalize();
  return 0;
}

