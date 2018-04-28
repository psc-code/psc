
#pragma once

// ----------------------------------------------------------------------
// cuda_heating_run_foil

struct cuda_heating_foil
{
  // params
  float zl;
  float zh;
  float xc;
  float yc;
  float rH;
  float T;
  float Mi;
  int kind;

  // state (FIXME, shouldn't be part of the interface)
  float fac;
  float heating_dt;
};

void cuda_heating_setup_foil(struct cuda_heating_foil *foil);

template<typename BS>
void cuda_heating_run_foil(cuda_mparticles<BS>* cmprts);

