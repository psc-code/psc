
#include "grid.hxx"
#include "fields.hxx"
#include "cuda_mfields.h"

class TestAccel
{
  enum { // FIXME, duplicated
    JXI, JYI, JZI,
    EX , EY , EZ ,
    HX , HY , HZ ,
    N_FIELDS,
  };

public:
  TestAccel()
  {
    double L = 1e10;

    grid_.dt = 1.;

    Grid_t::Patch patch{};
    grid_.ldims = { 1, 1, 1 };
    grid_.dx = { L, L, L };
    
    patch.xb = { 0., 0., 0. };
    patch.xe = { L,  L,  L  };
    grid_.patches.push_back(patch);

    cmflds_ = new cuda_mfields(grid_, N_FIELDS, { 1, 1, 1 });
    fields_single_t flds = cmflds_->get_host_fields();
    Fields3d<fields_single_t> F(flds);
    F(EX, 0,0,0) = 1;
    F(EX, 0,1,0) = 1;
    F(EX, 0,0,1) = 1;
    F(EX, 0,1,1) = 1;
    
    F(EY, 0,0,0) = 2;
    F(EY, 0,0,1) = 2;
    F(EY, 1,0,0) = 2;
    F(EY, 1,0,1) = 2;
    
    F(EZ, 0,0,0) = 3;
    F(EZ, 1,0,0) = 3;
    F(EZ, 0,1,0) = 3;
    F(EZ, 1,1,0) = 3;

    cmflds_->copy_to_device(0, flds, 0, N_FIELDS);
    cmflds_->dump("accel.fld.json");
    flds.dtor();
  };

  ~TestAccel()
  {
    delete cmflds_;
  }

  void run()
  {
  };

private:
  Grid_t grid_;
  cuda_mfields* cmflds_;
};

// ----------------------------------------------------------------------
// main

int
main(void)
{
  TestAccel test;
  test.run();
}
