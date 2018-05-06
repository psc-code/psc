
#pragma once

// ======================================================================
// cuda_collision

template<typename cuda_mparticles>
struct cuda_collision
{
  cuda_collision(int interval, double nu)
    : interval_{interval}, nu_{nu}
  {}
  
  void operator()(cuda_mparticles* cmprts)
  {}

private:
  int interval_;
  double nu_;
};

