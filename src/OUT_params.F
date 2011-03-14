
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

      ! FIXME, if needed, override defaults in OUT_fields_alloc!
      integer :: r1n, r1x, r2n, r2x, r3n, r3x

      end module OUT_params
