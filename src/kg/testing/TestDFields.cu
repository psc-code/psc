
#include <gtest/gtest.h>

#include "cuda_mfields.h"

#include <thrust/host_vector.h>

__global__
void set_raw(float *data)
{
  data[0] = 0;
  data[1] = 100;
  data[2] = 200;
  data[3] = 10;
  data[4] = 110;
  data[5] = 210;
}

__global__
void set_dfields(DFields d_flds)
{
  d_flds(0, 0, 0, 0) = 0;
  d_flds(0, 1, 0, 0) = 100;
  d_flds(0, 2, 0, 0) = 200;
  d_flds(0, 0, 1, 0) = 10;
  d_flds(0, 1, 1, 0) = 110;
  d_flds(0, 2, 1, 0) = 210;
}

TEST(DFields, Ctor)
{
  psc::device_vector<float> d_storage(6);
  
  auto d_flds = DFields{{{0, 0, 0}, {3, 2, 1}}, 1, thrust::raw_pointer_cast(d_storage.data())};
  //set_raw<<<1, 1>>>(d_flds.data());
  set_dfields<<<1, 1>>>(d_flds);
  
  thrust::host_vector<float> h_storage = d_storage;

#if 0
  auto h_flds = DFields{{{0, 0, 0}, {3, 2, 1}}, 1, h_storage.data()};
  for (int k = 0; k < 1; k++) {
    for (int j = 0; j < 2; j++) {
      for (int i = 0; i < 3; i++) {
	std::cout << h_flds(0, i,j,k) << " ";
      }
      std::cout << "\n";
    }
  }
#endif

  EXPECT_EQ(h_storage, std::vector<float>({0, 100, 200, 10, 110, 210}));
}