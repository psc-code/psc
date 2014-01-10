
#include <mrc_params.h>
#include <mrc_domain.h>
#include <mrc_fld.h>
#include <mrc_io.h>
#include <mrctest.h>
#include <mrc_trafo.h>

#include <stdio.h>
#include <string.h>
#include <assert.h>

static void
check_trafos(struct mrc_trafo *trafo1, struct mrc_trafo *trafo2)
{
  
  mrc_fld_foreach_patch(trafo1->_jac, p) {
    mrc_fld_foreach(trafo1->_jac, jx, jy, jz, 0, 0) {
      assert(TRAFO_CRD0(trafo1, jx, jy, jz, p) == TRAFO_CRD0(trafo2, jx, jy, jz, p));
      assert(TRAFO_CRD1(trafo1, jx, jy, jz, p) == TRAFO_CRD1(trafo2, jx, jy, jz, p));
      assert(TRAFO_CRD2(trafo1, jx, jy, jz, p) == TRAFO_CRD2(trafo2, jx, jy, jz, p));
      assert(TRAFO_JAC(trafo1, jx, jy, jz, p) == TRAFO_JAC(trafo2, jx, jy, jz, p));
      for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	  assert(TRAFO_EL(trafo1, i, j, jx, jy, jz, p) == TRAFO_EL(trafo2, i, j, jx, jy, jz, p));
	  assert(TRAFO_EU(trafo1, i, j, jx, jy, jz, p) == TRAFO_EU(trafo2, i, j, jx, jy, jz, p));
	  assert(TRAFO_GLL(trafo1, i, j, jx, jy, jz, p) == TRAFO_GLL(trafo2, i, j, jx, jy, jz, p));
	  assert(TRAFO_GUU(trafo1, i, j, jx, jy, jz, p) == TRAFO_GUU(trafo2, i, j, jx, jy, jz, p));
	  for (int k = 0; k < 3; k++) {
	    assert(TRAFO_GAM(trafo1, i, j, k, jx, jy, jz, p) == TRAFO_GAM(trafo2, i, j, k, jx, jy, jz, p));
	  }
	}
      }
      
    } mrc_fld_foreach_end;
  }
}

static void
test_read_write(struct mrc_domain *domain)
{
  struct mrc_io *io = mrc_io_create(mrc_domain_comm(domain));

  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "w", 0, 0.);
  mrc_io_write_path(io, "/domain", "domain", domain);
  mrc_io_close(io);
  mrc_io_destroy(io);

  io = mrc_io_create(mrc_domain_comm(domain));
  mrc_io_set_from_options(io);
  mrc_io_setup(io);
  mrc_io_open(io, "r", 0, 0.);
  struct mrc_domain *domain2 = mrc_io_read_path(io, "/domain", "domain", mrc_domain);
  mrc_io_close(io);
  mrc_io_destroy(io);

  mrctest_crds_compare(mrc_domain_get_crds(domain),
		       mrc_domain_get_crds(domain2));

  typedef struct mrc_trafo* (*dgt_t)(struct mrc_domain *, struct mrc_trafo **);
  dgt_t domain_get_trafo = (dgt_t) mrc_domain_get_method(domain, "get_trafo");
  struct mrc_trafo *trafo1, *trafo2;
  domain_get_trafo(domain, &trafo1);
  domain_get_trafo(domain2, &trafo2);
  check_trafos(trafo1, trafo2);

  mrc_domain_destroy(domain2);
}

int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  libmrc_params_init(argc, argv);

  struct mrc_domain *domain = mrc_domain_create(MPI_COMM_WORLD);
  mrc_domain_set_type(domain, "mb");
  mrc_domain_set_param_string(domain, "block_factory", "simple3d");
  struct mrc_crds *crds = mrc_domain_get_crds(domain);
  int testcase = 1;
  mrc_params_get_option_int("case", &testcase);

  switch (testcase) {
  case 1:
    mrc_crds_set_type(crds, "uniform");
    mrc_domain_set_from_options(domain);
    mrc_domain_setup(domain);
    test_read_write(domain);
    break;
  case 2: ;
    // THIS CASE WILL ALWAYS FAIL! TRAFO GETS ASSIGNED AT CREATION
    // NO WAY TO REGENERATE TRAFO!
    // Though I could totally do that, if I felt like it.
    mrc_crds_set_type(crds, "rectilinear");
    mrc_crds_set_param_int(crds, "sw", 2);
    mrc_domain_set_from_options(domain);
    mrc_domain_setup(domain);
    mrc_domain_view(domain);
    mrctest_set_crds_rectilinear_1(domain);
    test_read_write(domain);
    break;
  }
  mrc_domain_destroy(domain);

  MPI_Finalize();
}
