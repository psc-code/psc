
c     ======================================================================
      module OUT_params
c     ======================================================================
c     general I/O params, like paths, the sub-box, etc.

      implicit none

c r1n: min output range along x
c r1x: max output range along x
c r2n: min output range along y
c r2x: max output range along y
c r3n: min output range along z
c r3x: max output range along z

      integer :: r1n, r1x, r2n, r2x, r3n, r3x

c data_out: data out directory
c data_chk: data checkpointing directory

      character*(200) data_out
      character*(200) data_chk
      character*(200) data_processed

c nnp: number of timesteps for one full cycle
c np:  number of timesteps for time-averaging

      integer :: nnp, np

c shift_z: density function shift in dimensionless units
c shift_c: time step counter for moving the frame in the code
c shift_l: step counter for moving the frame in the code

      real(kind=8) :: shift_z
      integer :: shift_c, shift_l

      end module OUT_params

c     ----------------------------------------------------------------------
      subroutine OUT_params_set(nnp_, np_)
c     ----------------------------------------------------------------------
      use OUT_params
      use PSC_domain

      implicit none

      integer :: nnp_, np_

c I/O CONTROL PARAMETERS

      ! output range, change here
      r1n=i1n
      r1x=i1x
      r2n=i2n
      r2x=i2x
      r3n=i3n
      r3x=i3x

      r1n=max(i1n,r1n)
      r1x=min(i1x,r1x)
      r2n=max(i2n,r2n)
      r2x=min(i2x,r2x)
      r3n=max(i3n,r3n)
      r3x=min(i3x,r3x)

c data_out: path of data output directory
c data_chk: path of checkpointing directory

      data_out='.'
      data_chk='.'
      data_processed='.'
      
c shift_z: shifts the density function in the code

      shift_z=0.0

c shift_c: time step counter for shifting the frame in the code

      shift_c=0
      shift_l=0

c time averaging

      nnp = nnp_
      np = np_

      end subroutine OUT_params_set
