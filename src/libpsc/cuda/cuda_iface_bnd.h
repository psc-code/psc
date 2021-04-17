
#ifndef CUDA_IFACE_BND_H
#define CUDA_IFACE_BND_H

class MfieldsCuda;

// ----------------------------------------------------------------------
// routines for actual domain boundaries

void cuda_conducting_wall_H_lo_hi_y(MfieldsCuda& mflds, int p);
void cuda_conducting_wall_E_lo_hi_y(MfieldsCuda& mflds, int p);
void cuda_conducting_wall_J_lo_hi_y(MfieldsCuda& mflds, int p);
void cuda_conducting_wall_H_lo_y(MfieldsCuda& mflds, int p);
void cuda_conducting_wall_H_hi_y(MfieldsCuda& mflds, int p);
void cuda_conducting_wall_E_lo_y(MfieldsCuda& mflds, int p);
void cuda_conducting_wall_E_hi_y(MfieldsCuda& mflds, int p);
void cuda_conducting_wall_J_lo_y(MfieldsCuda& mflds, int p);
void cuda_conducting_wall_J_hi_y(MfieldsCuda& mflds, int p);

#endif
