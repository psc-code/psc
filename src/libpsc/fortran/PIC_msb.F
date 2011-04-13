c THIS SUBROUTINE DETERMINES A 3D SOLUTION OF MAXWELLS EQUATIONS.
c THE UPDATED E AND H-FIELDS ARE AT t=(n+0.5)*dt.


      subroutine PIC_msb_h(ex, ey, ez, hx, hy, hz, jxi, jyi, jzi)

      use VLA_variables
      use PSC_domain
      use PSC_patch
      use OUT_poyn

      implicit none

      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: ex, ey, ez
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: hx, hy, hz
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: jxi, jyi, jzi

      real(kind=8) :: cnx,cny,cnz

      real(kind=8) :: jx,jy,jz

      real(kind=8) :: exh,eyh,ezh,hxh,hyh,hzh
      real(kind=8) :: exk,eyk,ezk,hxk,hyk,hzk

      integer :: i1, i2, i3

c initialization

      cnx=0.5*dt/dx
      cny=0.5*dt/dy
      cnz=0.5*dt/dz


c B-field propagation E^(n+0.5), B^(n+0.5), j^(n+1.0), m^(n+0.5)
c -> E^(n+0.5), B^(n+1.0), j^(n+1.0), m^(n+0.5)


      hx2B=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               hxh=hx(i1,i2,i3)
               hxk=hx(i1,i2,i3)
     &             -cny*(ez(i1,i2+1,i3)-ez(i1,i2,i3))
     &             +cnz*(ey(i1,i2,i3+1)-ey(i1,i2,i3))
               hx2B=hx2B+dx*dy*dz*hxh*hxk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               hx(i1,i2,i3)=hx(i1,i2,i3)
     &                      -cny*(ez(i1,i2+1,i3)-ez(i1,i2,i3))
     &                      +cnz*(ey(i1,i2,i3+1)-ey(i1,i2,i3))
            enddo
         enddo
      enddo

      hy2B=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               hyh=hy(i1,i2,i3)
               hyk=hy(i1,i2,i3)
     &             -cnz*(ex(i1,i2,i3+1)-ex(i1,i2,i3))
     &             +cnx*(ez(i1+1,i2,i3)-ez(i1,i2,i3))
               hy2B=hy2B+dx*dy*dz*hyh*hyk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               hy(i1,i2,i3)=hy(i1,i2,i3)
     &                      -cnz*(ex(i1,i2,i3+1)-ex(i1,i2,i3))
     &                      +cnx*(ez(i1+1,i2,i3)-ez(i1,i2,i3))
            enddo
         enddo
      enddo

      hz2B=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               hzh=hz(i1,i2,i3)
               hzk=hz(i1,i2,i3)
     &             -cnx*(ey(i1+1,i2,i3)-ey(i1,i2,i3))
     &             +cny*(ex(i1,i2+1,i3)-ex(i1,i2,i3))
               hz2B=hz2B+dx*dy*dz*hzh*hzk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               hz(i1,i2,i3)=hz(i1,i2,i3)
     &                      -cnx*(ey(i1+1,i2,i3)-ey(i1,i2,i3))
     &                      +cny*(ex(i1,i2+1,i3)-ex(i1,i2,i3))
            enddo
         enddo
      enddo

      end subroutine PIC_msb_h


      subroutine PIC_msb_e(ex, ey, ez, hx, hy, hz, jxi, jyi, jzi)

      use VLA_variables
      use PSC_domain
      use PSC_patch
      use OUT_poyn

      implicit none

      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: ex, ey, ez
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: hx, hy, hz
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: jxi, jyi, jzi

      real(kind=8) :: cnx,cny,cnz

      real(kind=8) :: jx,jy,jz

      real(kind=8) :: exh,eyh,ezh,hxh,hyh,hzh
      real(kind=8) :: exk,eyk,ezk,hxk,hyk,hzk

      integer :: i1, i2, i3

c initialization

      cnx=0.5*dt/dx
      cny=0.5*dt/dy
      cnz=0.5*dt/dz

c E-field propagation E^(n+0.5), B^(n+1.0), j^(n+1.0) 
c -> E^(n+1.0), B^(n+1.0), j^(n+1.0)


      ex2B=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               jx=jxi(i1,i2,i3)
               exh=ex(i1,i2,i3)
               exk=ex(i1,i2,i3)
     &             +cny*(hz(i1,i2,i3)-hz(i1,i2-1,i3))
     &             -cnz*(hy(i1,i2,i3)-hy(i1,i2,i3-1))
     &             -0.5*dt*jx
               ex2B=ex2B+dx*dy*dz*exh*exk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               jx=jxi(i1,i2,i3)
               ex(i1,i2,i3)=ex(i1,i2,i3)
     &                       +cny*(hz(i1,i2,i3)-hz(i1,i2-1,i3))
     &                       -cnz*(hy(i1,i2,i3)-hy(i1,i2,i3-1))
     &                       -0.5*dt*jx
            enddo
         enddo
      enddo

      ey2B=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               jy=jyi(i1,i2,i3)
               eyh=ey(i1,i2,i3)
               eyk=ey(i1,i2,i3)
     &             +cnz*(hx(i1,i2,i3)-hx(i1,i2,i3-1))
     &             -cnx*(hz(i1,i2,i3)-hz(i1-1,i2,i3))
     &             -0.5*dt*jy
               ey2B=ey2B+dx*dy*dz*eyh*eyk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               jy=jyi(i1,i2,i3)
               ey(i1,i2,i3)=ey(i1,i2,i3)
     &                       +cnz*(hx(i1,i2,i3)-hx(i1,i2,i3-1))
     &                       -cnx*(hz(i1,i2,i3)-hz(i1-1,i2,i3))
     &                       -0.5*dt*jy
            enddo
         enddo
      enddo

      ez2B=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               jz=jzi(i1,i2,i3)
               ezh=ez(i1,i2,i3)
               ezk=ez(i1,i2,i3)
     &             +cnx*(hy(i1,i2,i3)-hy(i1-1,i2,i3))
     &             -cny*(hx(i1,i2,i3)-hx(i1,i2-1,i3))
     &             -0.5*dt*jz
               ez2B=ez2B+dx*dy*dz*ezh*ezk
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               jz=jzi(i1,i2,i3)
               ez(i1,i2,i3)=ez(i1,i2,i3)
     &                       +cnx*(hy(i1,i2,i3)-hy(i1-1,i2,i3))
     &                       -cny*(hx(i1,i2,i3)-hx(i1,i2-1,i3))
     &                       -0.5*dt*jz
            enddo
         enddo
      enddo


c energy conservation // FIXME, probably broken since ghost points haven't been set


      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               jx=jxi(i1,i2,i3)
               jy=jyi(i1,i2,i3)
               jz=jzi(i1,i2,i3)
               je=je+0.5*dx*dy*dz*(jx*ex(i1,i2,i3)
     &                            +jy*ey(i1,i2,i3)
     &                            +jz*ez(i1,i2,i3))
            enddo
         enddo
      enddo


      end subroutine PIC_msb_e

!     ------------------------------------------------------------
      subroutine PIC_fill_ghosts_h_b(hx, hy, hz, ex, ey, ez,
     &     jxi, jyi, jzi)
!     ------------------------------------------------------------
      use VLA_variables
      use PSC_domain
      use PSC_patch

      implicit none

      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: ex, ey, ez
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: hx, hy, hz
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: jxi, jyi, jzi

      integer :: i1, i2, i3
      real(kind=8) :: jx,jy,jz
      real(kind=8) :: lx,ly,lz
      real(kind=8) :: t,x,y,z

      real(kind=8) :: C_p_pulse_x1
      real(kind=8) :: C_s_pulse_x1
      real(kind=8) :: C_p_pulse_y1
      real(kind=8) :: C_s_pulse_y1
      real(kind=8) :: C_p_pulse_z1
      real(kind=8) :: C_s_pulse_z1
      real(kind=8) :: C_p_pulse_x2
      real(kind=8) :: C_s_pulse_x2
      real(kind=8) :: C_p_pulse_y2
      real(kind=8) :: C_s_pulse_y2
      real(kind=8) :: C_p_pulse_z2
      real(kind=8) :: C_s_pulse_z2

      lx=dt/dx
      ly=dt/dy
      lz=dt/dz

c open boundary at x1 plus incident laser pulse

      if (boundary_fld_x1==0 .and. i1mn==0) then       

      t=n*dt
      i1=i1mn-1
      x=i1*dx

      do i3=i3mn-1,i3mx+1
         z=i3*dz
         do i2=i2mn-1,i2mx+1
            y=i2*dy
            jy=jyi(i1,i2,i3)
            hz(i1-1,i2,i3)=(4.0*C_s_pulse_x1(x,y+0.5*dy,z,t)
     &                      -2.0*ey(i1,i2,i3)
     &                      -(1.0-lx)*hz(i1,i2,i3)
     &                      -lz*(hx(i1,i2,i3)-hx(i1,i2,i3-1))
     &                      -dt*jy)/(1.0+lx)
            jz=jzi(i1,i2,i3)
            hy(i1-1,i2,i3)=(-4.0*C_p_pulse_x1(x,y,z+0.5*dz,t)
     &                      +2.0*ez(i1,i2,i3)
     &                      -(1.0-lx)*hy(i1,i2,i3)
     &                      +ly*(hx(i1,i2,i3)-hx(i1,i2-1,i3))
     &                      +dt*jz)/(1.0+lx)
         enddo
      enddo
      endif

c open boundary at y1 plus incident laser pulse

      if (boundary_fld_y1==0 .and. i2mn==0) then       

      t=n*dt
      i2=i2mn ! ghost cell we're setting: (i2mn - 1) + 1/2
      y=i2*dy

      do i3=i3mn-1,i3mx+1
         z=i3*dz
         do i1=i1mn-1,i1mx+1
            x=i1*dx
            jz=jzi(i1,i2,i3)
            hx(i1,i2-1,i3)=(4.0*C_s_pulse_y1(x,y,z+0.5*dz,t)
     &                      -2.0*ez(i1,i2,i3)
     &                      -(1.0-ly)*hx(i1,i2,i3)
     &                      -lx*(hy(i1,i2,i3)-hy(i1-1,i2,i3))
     &                      +dt*jz)/(1.0+ly)
            jx=jxi(i1,i2,i3)
            hz(i1,i2-1,i3)=(-4.0*C_p_pulse_y1(x+0.5*dx,y,z,t)
     &                      +2.0*ex(i1,i2,i3)
     &                      -(1.0-ly)*hz(i1,i2,i3)
     &                      -lz*(hy(i1,i2,i3)-hy(i1,i2,i3-1))
     &                      +dt*jx)/(1.0+ly)
         enddo
      enddo
      endif

c open boundary at z1 plus incident laser pulse

      if (boundary_fld_z1==0 .and. i3mn==0) then       

      t=n*dt
      i3=i3mn-1 ! ghost cell we're setting: (i3mn - 1 - 1) + 1/2
      z=i3*dz

      do i2=i2mn-1,i2mx+1
         y=i2*dy
         do i1=i1mn-1,i1mx+1
            x=i1*dx
            jx=jxi(i1,i2,i3)
            hy(i1,i2,i3-1)=(+4.0*C_s_pulse_z1(x+0.5*dx,y,z,t)
     &                      -2.0*ex(i1,i2,i3)
     &                      -(1.0-lz)*hy(i1,i2,i3)
     &                      -ly*(hz(i1,i2,i3)-hz(i1,i2-1,i3))
     &                      +dt*jx)/(1.0+lz)
            jy=jyi(i1,i2,i3)
            hx(i1,i2,i3-1)=(-4.0*C_p_pulse_z1(x,y+0.5*dy,z,t)
     &                      +2.0*ey(i1,i2,i3)
     &                      -(1.0-lz)*hx(i1,i2,i3)
     &                      -lx*(hz(i1,i2,i3)-hz(i1-1,i2,i3))
     &                      -dt*jy)/(1.0+lz)
         enddo
      enddo
      endif

c open boundary at x2 plus incident laser pulse

      if (boundary_fld_x2==0 .and. i1mx==i1x) then       

      t=n*dt
      i1=i1mx
      x=(i1+1)*dx

      do i3=i3mn-1,i3mx+1
         z=i3*dz
         do i2=i2mn-1,i2mx+1
            y=i2*dy
            jy=jyi(i1+1,i2,i3)
            hz(i1+1,i2,i3)=(-4.0*C_s_pulse_x2(x,y+0.5*dy,z,t)
     &                      +2.0*ey(i1+1,i2,i3)
     &                      -(1.0-lx)*hz(i1,i2,i3)
     &                      -lz*(hx(i1+1,i2,i3)-hx(i1+1,i2,i3-1))
     &                      +dt*jy)/(1.0+lx)
            jz=jzi(i1+1,i2,i3)
            hy(i1+1,i2,i3)=(4.0*C_p_pulse_x2(x,y,z+0.5*dz,t)
     &                      -2.0*ez(i1+1,i2,i3)
     &                      -(1.0-lx)*hy(i1,i2,i3)
     &                      +ly*(hx(i1+1,i2,i3)-hx(i1+1,i2-1,i3))
     &                      -dt*jz)/(1.0+lx)
         enddo
      enddo
      endif

c open boundary at y2 plus incident laser pulse

      if (boundary_fld_y2==0 .and. i2mx==i2x) then       

      t=n*dt
      i2=i2mx
      y=(i2+1)*dy

      do i3=i3mn-1,i3mx+1
         z=i3*dz
         do i1=i1mn-1,i1mx+1
            x=i1*dx
            jz=jzi(i1,i2+1,i3)
            hx(i1,i2+1,i3)=(-4.0*C_s_pulse_y2(x,y,z+0.5*dz,t)
     &                    +2.0*ez(i1,i2+1,i3)
     &                    -(1.0-ly)*hx(i1,i2,i3)
     &                    +lx*(hy(i1,i2+1,i3)-hy(i1-1,i2+1,i3))
     &                    -dt*jz)/(1.0+ly)
            jx=jxi(i1,i2+1,i3)
            hz(i1,i2+1,i3)=(4.0*C_p_pulse_y2(x+0.5*dx,y,z,t)
     &                    -2.0*ex(i1,i2+1,i3)
     &                    -(1.0-ly)*hz(i1,i2,i3)
     &                    +lz*(hy(i1,i2+1,i3)-hy(i1,i2+1,i3-1))
     &                    -dt*jx)/(1.0+ly)
         enddo
      enddo
      endif

c open boundary at z2 plus incident laser pulse

      if (boundary_fld_z2==0 .and. i3mx==i3x) then       

      t=n*dt
      i3=i3mx ! ghost cell we're setting: (i3mx + 1) + 1/2
      z=(i3+1)*dz

      do i2=i2mn-1,i2mx+1
         y=i2*dy
         do i1=i1mn-1,i1mx+1
            x=i1*dx
            jx=jxi(i1,i2,i3+1)
            hy(i1,i2,i3+1)=(-4.0*C_s_pulse_z2(x+0.5*dx,y,z,t)
     &                    +2.0*ex(i1,i2,i3+1)
     &                    -(1.0-lz)*hy(i1,i2,i3)
     &                    +ly*(hz(i1,i2,i3)-hz(i1,i2-1,i3+1))
     &                    -dt*jx)/(1.0+lz)
            jy=jyi(i1,i2,i3+1)
            hx(i1,i2,i3+1)=(+4.0*C_p_pulse_z2(x,y+0.5*dy,z,t)
     &                    -2.0*ey(i1,i2,i3+1)
     &                    -(1.0-lz)*hx(i1,i2,i3)
     &                    +lx*(hz(i1,i2,i3+1)-hz(i1-1,i2,i3+1))
     &                    +dt*jy)/(1.0+lz)
         enddo
      enddo
      endif

      end subroutine PIC_fill_ghosts_h_b
