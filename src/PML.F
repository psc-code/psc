
c     ======================================================================
      module fields_pml
c     ======================================================================
      implicit none

! for the pml the declaration of d and h field is required
! dvx-hz: additional em field arrays
! dvx-dvz is chosen because of the existence of dx-dz for spatial range

      real(kind=8),allocatable,dimension(:,:,:) :: dvx,dvy,dvz   ! added by ab
      real(kind=8),allocatable,dimension(:,:,:) :: hx,hy,hz   ! added by ab

! electric and magnetic permittivity tensors

      real(kind=4),allocatable,dimension(:,:,:) :: eps, mu  ! added by ab

      end module fields_pml

c     ----------------------------------------------------------------------
      subroutine PML_fields_alloc
c     ----------------------------------------------------------------------

      use VLA_variables
      use fields_pml

      implicit none

      ! additional d field  - added by ab

      allocate(dvx(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(dvy(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(dvz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))


      ! additional h field - added by ab

      allocate(hx(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(hy(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(hz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))

      ! electric and magnetic permittivity tensors - added by ab

      allocate(eps(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(mu(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))

      call OUT_alloc_fields_pml

      dvx=0.0d0
      dvy=0.0d0
      dvz=0.0d0
      hx=0.0d0
      hy=0.0d0
      hz=0.0d0
      eps = 1.0
      mu = 1.0

      end subroutine PML_fields_alloc

c     ----------------------------------------------------------------------
      subroutine PML_fields_setup
c     ----------------------------------------------------------------------

      use VLA_variables
      use fields
      use fields_pml

      implicit none

      ! initialisation of additional fields - added by ab
      ! because of normalization these are equal to e and b fields

      dvx=ex
      dvy=ey
      dvz=ez
      hx=bx
      hy=by
      hz=bz
      
      ! initialisation of permittivity tensors, always set to 1.0
      ! due to scaling these are dimensionless
      
      eps = 1.0
      mu = 1.0

      call OUT_setup_fields_pml
      call OUT_poyn_setup

      end subroutine PML_fields_setup

c     ======================================================================
      module coeff_pml
c     ======================================================================
      implicit none

! cushion: thickness of buffer between pml region and free space
! pml: polynomial order of pml

      integer cushion, pml

! kappax_max, sigmax_max: attenuation coefficients in x direction
! kappay_max, sigmay_max: attenuation coefficients in y direction
! kappaz_max, sigmaz_max: attenuation coefficients in z direction

      real(kind=8) kappax_max,sigmax_max
      real(kind=8) kappay_max,sigmay_max
      real(kind=8) kappaz_max,sigmaz_max

! kappax, sigmax: coefficient arrays in x
! kappay, sigmay: coefficient arrays in y
! kappaz, sigmaz: coefficient arrays in z

      real(kind=8),allocatable,dimension(:) :: kappax,sigmax
      real(kind=8),allocatable,dimension(:) :: kappay,sigmay
      real(kind=8),allocatable,dimension(:) :: kappaz,sigmaz

! coefficient at integer position
      
      real(kind=8),allocatable,dimension(:) :: cxp, cxm
      real(kind=8),allocatable,dimension(:) :: fbx, fcx, fdx, fex
      real(kind=8),allocatable,dimension(:) :: cyp, cym
      real(kind=8),allocatable,dimension(:) :: fby, fcy, fdy, fey
      real(kind=8),allocatable,dimension(:) :: czp, czm
      real(kind=8),allocatable,dimension(:) :: fbz, fcz, fdz, fez

! coefficient at position moved by half space

      real(kind=8),allocatable,dimension(:) :: bxp, bxm
      real(kind=8),allocatable,dimension(:) :: gbx, gcx, gdx, gex
      real(kind=8),allocatable,dimension(:) :: byp, bym
      real(kind=8),allocatable,dimension(:) :: gby, gcy, gdy, gey
      real(kind=8),allocatable,dimension(:) :: bzp, bzm
      real(kind=8),allocatable,dimension(:) :: gbz, gcz, gdz, gez

      end module coeff_pml

c     ----------------------------------------------------------------------
      subroutine PML_coeff_init_param
c     ----------------------------------------------------------------------

      use VLA_variables
      use coeff_pml

      implicit none
      
      thick = 10
      cushion = int(thick/3)
      size = thick+cushion
      pml = 3
       
! PML PARAMETERS - added by ab

! deltax: thickness of pml in x direction
! deltay: thickness of pml in y direction
! deltaz: thickness of pml in z direction

      deltax = thick*dx
      deltay = thick*dy
      deltaz = thick*dz

! attenuation factors
      kappax_max = 1.34
      sigmax_max = 13.8*eps0/deltax

      kappay_max = 1.34
      sigmay_max = 13.8*eps0/deltay

      kappaz_max = 1.34
      sigmaz_max = 13.8*eps0/deltaz

      end subroutine PML_coeff_init_param

c     ----------------------------------------------------------------------
      subroutine PML_coeff_setup
c     ----------------------------------------------------------------------
c     deriving pml attenuation coefficients

      use VLA_variables
      use coeff_pml

      implicit none

! not moved by half step
      
      ! x

      do i1 = i1mn-rd1, i1mx+rd1
         kappax(i1) = 1.0
         sigmax(i1) = 0.0
         cxp(i1) = 2*eps0*kappax(i1)+sigmax(i1)*dt
         cxm(i1) = 2*eps0*kappax(i1)-sigmax(i1)*dt
         fbx(i1) = 2*eps0*kappax(i1)
         fcx(i1) = cxm(i1)/cxp(i1)
         fdx(i1) = 2*eps0*dt/cxp(i1)
         fex(i1) = 1.0/cxp(i1)
      end do

      ! y
      
       do i2 = i2mn-rd2, i2mx+rd2
         kappay(i2) = 1.0
         sigmay(i2) = 0.0
         cyp(i2) = 2*eps0*kappay(i2)+sigmay(i2)*dt
         cym(i2) = 2*eps0*kappay(i2)-sigmay(i2)*dt
         fby(i2) = 2*eps0*kappay(i2)
         fcy(i2) = cym(i2)/cyp(i2)
         fdy(i2) = 2*eps0*dt/cyp(i2)
         fey(i2) = 1.0/cyp(i2)
      end do

      ! z

      do i3 = i3mn-rd3, i3mx+rd3
         kappaz(i3) = 1.0
         sigmaz(i3) = 0.0
         czp(i3) = 2*eps0*kappaz(i3)+sigmaz(i3)*dt
         czm(i3) = 2*eps0*kappaz(i3)-sigmaz(i3)*dt
         fbz(i3) = 2*eps0*kappaz(i3)
         fcz(i3) = czm(i3)/czp(i3)
         fdz(i3) = 2*eps0*dt/czp(i3)
         fez(i3) = 1.0/czp(i3)
      end do


! moved by half step

      ! x

      do i1 = i1mn-rd1, i1mx+rd1
         kappax(i1) = 1.0
         sigmax(i1) = 0.0
         bxp(i1) = 2*eps0*kappax(i1)+sigmax(i1)*dt
         bxm(i1) = 2*eps0*kappax(i1)-sigmax(i1)*dt
         gbx(i1) = 2*eps0*kappax(i1)
         gcx(i1) = bxm(i1)/bxp(i1)
         gdx(i1) = 2*eps0*dt/bxp(i1)
         gex(i1) = 1.0/bxp(i1)
      end do

      ! y

      do i2 = i2mn-rd2, i2mx+rd2
         kappay(i2) = 1.0
         sigmay(i2) = 0.0
         byp(i2) = 2*eps0*kappay(i2)+sigmay(i2)*dt
         bym(i2) = 2*eps0*kappay(i2)-sigmay(i2)*dt
         gby(i2) = 2*eps0*kappay(i2)
         gcy(i2) = bym(i2)/byp(i2)
         gdy(i2) = 2*eps0*dt/byp(i2)
         gey(i2) = 1.0/byp(i2)
      end do

      ! z

      do i3 = i3mn-rd3, i3mx+rd3
         kappaz(i3) = 1.0
         sigmaz(i3) = 0.0
         bzp(i3) = 2*eps0*kappaz(i3)+sigmaz(i3)*dt
         bzm(i3) = 2*eps0*kappaz(i3)-sigmaz(i3)*dt
         gbz(i3) = 2*eps0*kappaz(i3)
         gcz(i3) = bzm(i3)/bzp(i3)
         gdz(i3) = 2*eps0*dt/bzp(i3)
         gez(i3) = 1.0/bzp(i3)
      end do

      end subroutine PML_coeff_setup
