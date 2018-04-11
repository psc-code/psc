
#pragma once

// ----------------------------------------------------------------------
// cuda_push_mprts

template<typename BS>
struct cuda_mparticles;

template<typename BS>
void cuda_push_mprts_yz(cuda_mparticles<BS>* cmprts, struct cuda_mfields *cmflds,
			const int bs[3], bool ip_ec, bool deposit_vb_3d, bool currmem_global);

void cuda_push_mprts_xyz(cuda_mparticles<BS144>*cmprts, struct cuda_mfields *cmflds);

