
c     ======================================================================
      module fields_pml
c     ======================================================================
      implicit none

! for the pml the declaration of d and h field is required
! dvx-hz: additional em field arrays
! dvx-dvz is chosen because of the existence of dx-dz for spatial range

      real(kind=8),allocatable,dimension(:,:,:) :: dvx,dvy,dvz   ! added by ab
      real(kind=8),allocatable,dimension(:,:,:) :: bx,by,bz   ! added by ab

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


      ! additional b field - added by ab

      allocate(bx(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(by(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(bz(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))

      ! electric and magnetic permittivity tensors - added by ab

      allocate(eps(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))
      allocate(mu(i1mn-rd1:i1mx+rd1,i2mn-rd2:i2mx+rd2,
     &            i3mn-rd3:i3mx+rd3))

      dvx=0.0d0
      dvy=0.0d0
      dvz=0.0d0
      bx=0.0d0
      by=0.0d0
      bz=0.0d0
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
      bx=hx
      by=hy
      bz=hz
      
      ! initialisation of permittivity tensors, always set to 1.0
      ! due to scaling these are dimensionless
      
      eps = 1.0
      mu = 1.0

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
      subroutine PML_coeff_alloc
c     ----------------------------------------------------------------------

      use VLA_variables
      use coeff_pml

      implicit none

! array allocation
! size of pml should be equal to size of em field arrays

! x
      allocate(kappax(i1mn-rd1:i1mx+rd1))
      allocate(sigmax(i1mn-rd1:i1mx+rd1))
      
      allocate(cxp(i1mn-rd1:i1mx+rd1))
      allocate(cxm(i1mn-rd1:i1mx+rd1))
     
      allocate(fbx(i1mn-rd1:i1mx+rd1))
      allocate(fcx(i1mn-rd1:i1mx+rd1))
      allocate(fdx(i1mn-rd1:i1mx+rd1))
      allocate(fex(i1mn-rd1:i1mx+rd1))

      allocate(bxp(i1mn-rd1:i1mx+rd1))
      allocate(bxm(i1mn-rd1:i1mx+rd1))
      
      allocate(gbx(i1mn-rd1:i1mx+rd1))
      allocate(gcx(i1mn-rd1:i1mx+rd1))
      allocate(gdx(i1mn-rd1:i1mx+rd1))
      allocate(gex(i1mn-rd1:i1mx+rd1))

! y
      
      allocate(kappay(i2mn-rd2:i2mx+rd2))
      allocate(sigmay(i2mn-rd2:i2mx+rd2))
      
      allocate(cyp(i2mn-rd2:i2mx+rd2))
      allocate(cym(i2mn-rd2:i2mx+rd2))
      
      allocate(fby(i2mn-rd2:i2mx+rd2))
      allocate(fcy(i2mn-rd2:i2mx+rd2))
      allocate(fdy(i2mn-rd2:i2mx+rd2))
      allocate(fey(i2mn-rd2:i2mx+rd2))

      allocate(byp(i2mn-rd2:i2mx+rd2))
      allocate(bym(i2mn-rd2:i2mx+rd2))
      
      allocate(gby(i2mn-rd2:i2mx+rd2))
      allocate(gcy(i2mn-rd2:i2mx+rd2))
      allocate(gdy(i2mn-rd2:i2mx+rd2))
      allocate(gey(i2mn-rd2:i2mx+rd2))

! z
      
      allocate(kappaz(i3mn-rd3:i3mx+rd3))
      allocate(sigmaz(i3mn-rd3:i3mx+rd3))
      
      allocate(czp(i3mn-rd3:i3mx+rd3))
      allocate(czm(i3mn-rd3:i3mx+rd3))
      
      allocate(fbz(i3mn-rd3:i3mx+rd3))
      allocate(fcz(i3mn-rd3:i3mx+rd3))
      allocate(fdz(i3mn-rd3:i3mx+rd3))
      allocate(fez(i3mn-rd3:i3mx+rd3))

      allocate(bzp(i3mn-rd3:i3mx+rd3))
      allocate(bzm(i3mn-rd3:i3mx+rd3))
      
      allocate(gbz(i3mn-rd3:i3mx+rd3))
      allocate(gcz(i3mn-rd3:i3mx+rd3))
      allocate(gdz(i3mn-rd3:i3mx+rd3))
      allocate(gez(i3mn-rd3:i3mx+rd3))

      end subroutine PML_coeff_alloc

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
      subroutine SET_param_pml(thick_, cushion_, size_, pml_)
c     ----------------------------------------------------------------------
c     wrapper for C to be able to set pml parameters

      use VLA_variables
      use PIC_variables
      use coeff_pml

      implicit none

      integer :: thick_, cushion_, size_, pml_

      thick = thick_
      cushion = cushion_
      size = size_
      pml = pml_

      end subroutine SET_param_pml

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
