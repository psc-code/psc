
c THIS SUBROUTINE DETERMINES A 3D SOLUTION OF MAXWELLS EQUATIONS.
c THE UPDATED E,D,H AND B-FIELDS ARE AT t=(n+1.0)*dt.


      subroutine PIC_pml_msb(ex, ey, ez, hx, hy, hz, dvx, dvy, dvz,
     &     bx, by, bz, jxi, jyi, jzi, eps, mu)

      use VLA_variables
      use PSC_domain
      use PSC_patch
      use coeff_pml
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
     &                       i3mn-rd3:i3mx+rd3) :: dvx, dvy, dvz
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: bx, by, bz
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: jxi, jyi, jzi
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: eps, mu

      real(kind=8) :: cnx,cny,cnz
      real(kind=8) :: lx,ly,lz
      real(kind=8) :: jx,jy,jz

      real(kind=8) :: C_p_pulse_z1
      real(kind=8) :: C_s_pulse_z1
      real(kind=8) :: C_p_pulse_z2
      real(kind=8) :: C_s_pulse_z2
      real(kind=8) :: C_p_pulse_y1
      real(kind=8) :: C_s_pulse_y1
      real(kind=8) :: C_p_pulse_y2
      real(kind=8) :: C_s_pulse_y2
      real(kind=8) :: C_p_pulse_x1
      real(kind=8) :: C_s_pulse_x1
      real(kind=8) :: C_p_pulse_x2
      real(kind=8) :: C_s_pulse_x2
      real(kind=8) :: t,x,y,z

      real(kind=8) :: dns,dls,ens,els    ! added by ab
      real(kind=8) :: bns,bls,hns,hls    ! added by ab
      real(kind=8) :: dstore,bstore      ! added by ab

      integer :: i1,i2,i3

c initialization


      lx=dt/dx
      ly=dt/dy
      lz=dt/dz

      cnx=0.5*lx
      cny=0.5*ly
      cnz=0.5*lz

      dns=0.0
      dls=0.0
      ens=0.0
      els=0.0
      bns=0.0
      bls=0.0
      hns=0.0
      hls=0.0
      dstore=0.0
      bstore=0.0


c B-field propagation E^(n+0.5), B^(n+0.5), j^(n+1.0), m^(n+0.5)
c -> E^(n+0.5), B^(n+1.0), j^(n+1.0), m^(n+0.5)


      hx2B=0.0
      bx2B=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               hls=hx(i1,i2,i3)
               bls=bx(i1,i2,i3)
               bns = (2*gcy(i2)*bx(i1,i2,i3) 
     &              - gdy(i2)*((ez(i1,i2+1,i3)-ez(i1,i2,i3))/dy 
     &              - (ey(i1,i2,i3+1)-ey(i1,i2,i3))/dz))/(1+gcy(i2))
               hns = 2*(gcz(i3)*hx(i1,i2,i3) 
     &              + gez(i3)*(fbx(i1)*bns-cxm(i1)*bls)
     &              /mu(i1,i2,i3))/(1+gcz(i3))
               bx2B=bx2B+dx*dy*dz*bls*bns
               hx2B=hx2B+dx*dy*dz*hls*hns
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               bstore = bx(i1,i2,i3)
               bx(i1,i2,i3) = (2*gcy(i2)*bx(i1,i2,i3) 
     &              - gdy(i2)*((ez(i1,i2+1,i3)-ez(i1,i2,i3))/dy 
     &              - (ey(i1,i2,i3+1)-ey(i1,i2,i3))/dz))/(1+gcy(i2))          
               hx(i1,i2,i3) = 2*(gcz(i3)*hx(i1,i2,i3) 
     &              + gez(i3)*(fbx(i1)*bx(i1,i2,i3)-cxm(i1)*bstore)
     &              /mu(i1,i2,i3))/(1+gcz(i3))
            enddo
         enddo
      enddo

      hy2B=0.0
      by2B=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               hls=hy(i1,i2,i3)
               bls=by(i1,i2,i3)
               bns = (2*gcz(i3)*by(i1,i2,i3) 
     &              - gdz(i3)*((ex(i1,i2,i3+1)-ex(i1,i2,i3))/dz 
     &              - (ez(i1+1,i2,i3)-ez(i1,i2,i3))/dx))/(1+gcz(i3))
               hns = 2*(gcx(i1)*hy(i1,i2,i3) 
     &              + gex(i1)*(fby(i2)*bns-cym(i2)*bls)
     &              /mu(i1,i2,i3))/(1+gcx(i1)) 
               by2B=by2B+dx*dy*dz*bls*bns
               hy2B=hy2B+dx*dy*dz*hls*hns
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               bstore = by(i1,i2,i3)
               by(i1,i2,i3) = (2*gcz(i3)*by(i1,i2,i3) 
     &              - gdz(i3)*((ex(i1,i2,i3+1)-ex(i1,i2,i3))/dz 
     &              - (ez(i1+1,i2,i3)-ez(i1,i2,i3))/dx))/(1+gcz(i3))    
               hy(i1,i2,i3) = 2*(gcx(i1)*hy(i1,i2,i3) 
     &              + gex(i1)*(fby(i2)*by(i1,i2,i3)-cym(i2)*bstore)
     &              /mu(i1,i2,i3))/(1+gcx(i1))
            enddo
         enddo
      enddo

      hz2B=0.0
      bz2B=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               hls=hz(i1,i2,i3)
               bls=bz(i1,i2,i3)
               bns = (2*gcx(i1)*bz(i1,i2,i3) 
     &              - gdx(i1)*((ey(i1+1,i2,i3)-ey(i1,i2,i3))/dx 
     &              - (ex(i1,i2+1,i3)-ex(i1,i2,i3))/dy))/(1+gcx(i1))
               hns = 2*(gcy(i2)*hz(i1,i2,i3) 
     &              + gey(i2)*(fbz(i3)*bns-czm(i3)*bls)
     &             /mu(i1,i2,i3))/(1+gcy(i2))               
               bz2B=bz2B+dx*dy*dz*bls*bns
               hz2B=hz2B+dx*dy*dz*hls*hns
            enddo
         enddo
      enddo
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               bstore = bz(i1,i2,i3)
               bz(i1,i2,i3) = (2*gcx(i1)*bz(i1,i2,i3) 
     &              - gdx(i1)*((ey(i1+1,i2,i3)-ey(i1,i2,i3))/dx 
     &              - (ex(i1,i2+1,i3)-ex(i1,i2,i3))/dy))/(1+gcx(i1))
               hz(i1,i2,i3) = 2*(gcy(i2)*hz(i1,i2,i3) 
     &              + gey(i2)*(fbz(i3)*bz(i1,i2,i3)-czm(i3)*bstore)
     &             /mu(i1,i2,i3)) /(1+gcy(i2))
            enddo
         enddo
      enddo


c open boundary at x1 plus incident laser pulse (ppol polarization ez)


      if (boundary_fld_x1==0.or.boundary_fld_x1==3) then 

      t=n*dt
      i1=i1mn-1
      x=i1*dx

      do i3=i3mn-1,i3mx+1
         z=i3*dz
         do i2=i2mn-1,i2mx+1
            y=i2*dy
            jy=jyi(i1,i2,i3)
            bz(i1-1,i2,i3)=(4.0*C_s_pulse_x1(x,y+0.5*dy,z,t)
     &                      -2.0*ey(i1,i2,i3)
     &                      -(1.0-lx)*bz(i1,i2,i3)
     &                      -lz*(bx(i1,i2,i3)-bx(i1,i2,i3-1))
     &                      -dt*jy)/(1.0+lx)
            hz(i1-1,i2,i3)=bz(i1-1,i2,i3)         ! added by ab
         enddo
      enddo

      t=n*dt
      i1=i1mn-1
      x=i1*dx

      do i3=i3mn-1,i3mx+1
         z=i3*dz
         do i2=i2mn-1,i2mx+1
            y=i2*dy
            jz=jzi(i1,i2,i3)
            by(i1-1,i2,i3)=(-4.0*C_p_pulse_x1(x,y,z+0.5*dz,t)
     &                      +2.0*ez(i1,i2,i3)
     &                      -(1.0-lx)*by(i1,i2,i3)
     &                      +ly*(bx(i1,i2,i3)-bx(i1,i2-1,i3))
     &                      +dt*jz)/(1.0+lx)
            hy(i1-1,i2,i3)=by(i1-1,i2,i3)        ! added by ab
         enddo
      enddo
      endif


c open boundary at x2 plus incident laser pulse


      if (boundary_fld_x2==0.or.boundary_fld_x2==3) then 

      t=n*dt
      i1=i1mx+1
      x=i1*dx

      do i3=i3mn-1,i3mx+1
         z=i3*dz
         do i2=i2mn-1,i2mx+1
            y=i2*dy
            jy=jyi(i1,i2,i3)
            bz(i1,i2,i3)=(-4.0*C_s_pulse_x2(x,y+0.5*dy,z,t)
     &                      +2.0*ey(i1,i2,i3)
     &                      -(1.0-lx)*bz(i1-1,i2,i3)
     &                      -lz*(bx(i1,i2,i3)-bx(i1,i2,i3-1))
     &                      +dt*jy)/(1.0+lx)
            hz(i1,i2,i3)=bz(i1,i2,i3)         ! added by ab
         enddo
      enddo

      t=n*dt
      i1=i1mx+1
      x=i1*dx

      do i3=i3mn-1,i3mx+1
         z=i3*dz
         do i2=i2mn-1,i2mx+1
            y=i2*dy
            jz=jzi(i1,i2,i3)
            by(i1,i2,i3)=(4.0*C_p_pulse_x2(x,y,z+0.5*dz,t)
     &                      -2.0*ez(i1,i2,i3)
     &                      -(1.0-lx)*by(i1-1,i2,i3)
     &                      +ly*(bx(i1,i2,i3)-bx(i1,i2-1,i3))
     &                      -dt*jz)/(1.0+lx)
            hy(i1,i2,i3)=by(i1,i2,i3)         ! added by ab
         enddo
      enddo
      endif


c open boundary at y1 plus incident laser pulse (ppol polarization ex)


      if (boundary_fld_y1==0.or.boundary_fld_y1==3) then 

      t=n*dt
      i2=i2mn-1
      y=i2*dy

      do i3=i3mn-1,i3mx+1
         z=i3*dz
         do i1=i1mn-1,i1mx+1
            x=i1*dx
            jz=jzi(i1,i2,i3)
            bx(i1,i2-1,i3)=(4.0*C_s_pulse_y1(x,y,z+0.5*dz,t)
     &                      -2.0*ez(i1,i2,i3)
     &                      -(1.0-ly)*bx(i1,i2,i3)
     &                      -lx*(by(i1,i2,i3)-by(i1-1,i2,i3))
     &                      +dt*jz)/(1.0+ly)
            hx(i1,i2-1,i3)=bx(i1,i2-1,i3)         ! added by ab
         enddo
      enddo

      t=n*dt
      i2=i2mn-1
      y=i2*dy

      do i3=i3mn-1,i3mx+1
         z=i3*dz
         do i1=i1mn-1,i1mx+1
            x=i1*dx
            jx=jxi(i1,i2,i3)
            bz(i1,i2-1,i3)=(-4.0*C_p_pulse_y1(x+0.5*dx,y,z,t)
     &                      +2.0*ex(i1,i2,i3)
     &                      -(1.0-ly)*bz(i1,i2,i3)
     &                      -lz*(by(i1,i2,i3)-by(i1,i2,i3-1))
     &                      +dt*jx)/(1.0+ly)
            hz(i1,i2-1,i3)=bz(i1,i2-1,i3)        ! added by ab
         enddo
      enddo
      endif


c open boundary at y2 plus incident laser pulse


      if (boundary_fld_y2==0.or.boundary_fld_y2==3) then 

      t=n*dt
      i2=i2mx+1
      y=i2*dy

      do i3=i3mn-1,i3mx+1
         z=i3*dz
         do i1=i1mn-1,i1mx+1
            x=i1*dx
            jz=jzi(i1,i2,i3)
            bx(i1,i2,i3)=(-4.0*C_s_pulse_y2(x,y,z+0.5*dz,t)
     &                    +2.0*ez(i1,i2,i3)
     &                    -(1.0-ly)*bx(i1,i2-1,i3)
     &                    +lx*(by(i1,i2,i3)-by(i1-1,i2,i3))
     &                    -dt*jz)/(1.0+ly)
            hx(i1,i2,i3)=bx(i1,i2,i3)              ! added by ab
         enddo
      enddo

      t=n*dt
      i2=i2mx+1
      y=i2*dy

      do i3=i3mn-1,i3mx+1
         z=i3*dz
         do i1=i1mn-1,i1mx+1
            x=i1*dx
            jx=jxi(i1,i2,i3)
            bz(i1,i2,i3)=(4.0*C_p_pulse_y2(x+0.5*dx,y,z,t)
     &                    -2.0*ex(i1,i2,i3)
     &                    -(1.0-ly)*bz(i1,i2-1,i3)
     &                    +lz*(by(i1,i2,i3)-by(i1,i2,i3-1))
     &                    -dt*jx)/(1.0+ly)
            hz(i1,i2,i3)=bz(i1,i2,i3)              ! added by ab
         enddo
      enddo
      endif


c open boundary at z1 plus incident laser pulse (ppol polarization ey)


      if (boundary_fld_z1==0.or.boundary_fld_z1==3) then 

      t=n*dt
      i3=i3mn-1
      z=i3*dz

      do i2=i2mn-1,i2mx+1
         y=i2*dy
         do i1=i1mn-1,i1mx+1
            x=i1*dx
            jx=jxi(i1,i2,i3)
            by(i1,i2,i3-1)=(+4.0*C_s_pulse_z1(x+0.5*dx,y,z,t)
     &                      -2.0*ex(i1,i2,i3)
     &                      -(1.0-lz)*by(i1,i2,i3)
     &                      -ly*(bz(i1,i2,i3)-bz(i1,i2-1,i3))
     &                      +dt*jx)/(1.0+lz)
            hy(i1,i2,i3-1)=by(i1,i2,i3-1)       ! added by ab
         enddo
      enddo

      t=n*dt
      i3=i3mn-1
      z=i3*dz

      do i2=i2mn-1,i2mx+1
         y=i2*dy
         do i1=i1mn-1,i1mx+1
            x=i1*dx
            jy=jyi(i1,i2,i3)
            bx(i1,i2,i3-1)=(-4.0*C_p_pulse_z1(x,y+0.5*dy,z,t)
     &                      +2.0*ey(i1,i2,i3)
     &                      -(1.0-lz)*bx(i1,i2,i3)
     &                      -lx*(bz(i1,i2,i3)-bz(i1-1,i2,i3))
     &                      -dt*jy)/(1.0+lz)
            hx(i1,i2,i3-1)=bx(i1,i2,i3-1)       ! added by ab
         enddo
      enddo
      endif


c open boundary at z2 plus incident laser pulse


      if (boundary_fld_z2==0.or.boundary_fld_z2==3) then 

      t=n*dt
      i3=i3mx+1
      z=i3*dz

      do i2=i2mn-1,i2mx+1
         y=i2*dy
         do i1=i1mn-1,i1mx+1
            x=i1*dx
            jx=jxi(i1,i2,i3)
            by(i1,i2,i3)=(-4.0*C_s_pulse_z2(x+0.5*dx,y,z,t)
     &                    +2.0*ex(i1,i2,i3)
     &                    -(1.0-lz)*by(i1,i2,i3-1)
     &                    +ly*(bz(i1,i2,i3)-bz(i1,i2-1,i3))
     &                    -dt*jx)/(1.0+lz)
            hy(i1,i2,i3)=by(i1,i2,i3)
         enddo
      enddo

      t=n*dt
      i3=i3mx+1
      z=i3*dz

      do i2=i2mn-1,i2mx+1
         y=i2*dy
         do i1=i1mn-1,i1mx+1
            x=i1*dx
            jy=jyi(i1,i2,i3)
            bx(i1,i2,i3)=(+4.0*C_p_pulse_z2(x,y+0.5*dy,z,t)
     &                    -2.0*ey(i1,i2,i3)
     &                    -(1.0-lz)*bx(i1,i2,i3-1)
     &                    +lx*(bz(i1,i2,i3)-bz(i1-1,i2,i3))
     &                    +dt*jy)/(1.0+lz)
            hx(i1,i2,i3)=bx(i1,i2,i3)
         enddo
      enddo
      endif


      call PIC_fex(bx)
      call PIC_fey(bx)
      call PIC_fez(bx)
      call PIC_fex(hx)
      call PIC_fey(hx)
      call PIC_fez(hx)
      call PIC_fex(by)
      call PIC_fey(by)
      call PIC_fez(by)
      call PIC_fex(hy)
      call PIC_fey(hy)
      call PIC_fez(hy)
      call PIC_fex(bz)
      call PIC_fey(bz)
      call PIC_fez(bz)
      call PIC_fex(hz)
      call PIC_fey(hz)
      call PIC_fez(hz)


c E-field propagation E^(n+0.5), B^(n+1.0), j^(n+1.0) 
c -> E^(n+1.0), B^(n+1.0), j^(n+1.0)


      dx2B=0.0
      ex2B=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               jx=jxi(i1,i2,i3)
               els=ex(i1,i2,i3)
               dls = dvx(i1,i2,i3)
               dns = (2*fcy(i2)*dvx(i1,i2,i3) 
     &              + fdy(i2)*((hz(i1,i2,i3)-hz(i1,i2-1,i3))/dy 
     &              - (hy(i1,i2,i3)-hy(i1,i2,i3-1))/dz-jx))/(1+fcy(i2))        
               ens = 2*(fcz(i3)*ex(i1,i2,i3) 
     &              + fez(i3)*(gbx(i1)*dns-bxm(i1)*dls)
     &              /eps(i1,i2,i3))/(1+fcz(i3))
               ex2B=ex2B+dx*dy*dz*els*ens
               dx2B=dx2B+dx*dy*dz*dls*dns
            enddo
         enddo
      enddo

      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               jx=jxi(i1,i2,i3)
               dstore = dvx(i1,i2,i3)
               dvx(i1,i2,i3) = (2*fcy(i2)*dvx(i1,i2,i3) 
     &              + fdy(i2)*((hz(i1,i2,i3)-hz(i1,i2-1,i3))/dy 
     &              - (hy(i1,i2,i3)-hy(i1,i2,i3-1))/dz-jx))/(1+fcy(i2))        
               ex(i1,i2,i3) = 2*(fcz(i3)*ex(i1,i2,i3) 
     &              + fez(i3)*(gbx(i1)*dvx(i1,i2,i3)-bxm(i1)*dstore)
     &              /eps(i1,i2,i3))/(1+fcz(i3))               
            enddo
         enddo
      enddo

      dy2B=0.0
      ey2B=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               jy=jyi(i1,i2,i3)
               els=ey(i1,i2,i3)
               dls = dvy(i1,i2,i3)
               dns = (2*fcz(i3)*dvy(i1,i2,i3) 
     &              + fdz(i3)*((hx(i1,i2,i3)-hx(i1,i2,i3-1))/dz 
     &              - (hz(i1,i2,i3)-hz(i1-1,i2,i3))/dx-jy))
     &              /(1+fcz(i3))        
               ens = 2*(fcx(i1)*ey(i1,i2,i3) 
     &              + fex(i1)*(gby(i2)*dns-bym(i2)*dls)
     &              /eps(i1,i2,i3))/(1+fcx(i1))
               ey2B=ey2B+dx*dy*dz*els*ens
               dy2B=dy2B+dx*dy*dz*dls*dns
            enddo
         enddo
      enddo
 
      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               jy=jyi(i1,i2,i3)
               dstore = dvy(i1,i2,i3)
               dvy(i1,i2,i3) = (2*fcz(i3)*dvy(i1,i2,i3) 
     &              + fdz(i3)*((hx(i1,i2,i3)-hx(i1,i2,i3-1))/dz 
     &              - (hz(i1,i2,i3)-hz(i1-1,i2,i3))/dx-jy))
     &              /(1+fcz(i3))        
               ey(i1,i2,i3) = 2*(fcx(i1)*ey(i1,i2,i3) 
     &              + fex(i1)*(gby(i2)*dvy(i1,i2,i3)-bym(i2)*dstore)
     &              /eps(i1,i2,i3))/(1+fcx(i1))
            enddo
         enddo
      enddo

      dz2B=0.0
      ez2B=0.0
      do i3=i3mn,i3mx
         do i2=i2mn,i2mx
            do i1=i1mn,i1mx
               jz=jzi(i1,i2,i3)
               els=ez(i1,i2,i3)
               dls = dvz(i1,i2,i3)
               dns = (2*fcx(i1)*dvz(i1,i2,i3) 
     &              + fdx(i1)*((hy(i1,i2,i3)-hy(i1-1,i2,i3))/dx 
     &              - (hx(i1,i2,i3)-hx(i1,i2-1,i3))/dy-jz))/(1+fcx(i1))
               ens = 2*(fcy(i2)*ez(i1,i2,i3) 
     &              + fey(i2)*(gbz(i3)*dns-bzm(i3)*dls)
     &              /eps(i1,i2,i3))/(1+fcy(i2))               
               ez2B=ez2B+dx*dy*dz*els*ens
               dz2B=dz2B+dx*dy*dz*dls*dns
            enddo
         enddo
      enddo

      do i3=i3mn-1,i3mx+1
         do i2=i2mn-1,i2mx+1
            do i1=i1mn-1,i1mx+1
               jz=jzi(i1,i2,i3)
               dstore = dvz(i1,i2,i3)
               dvz(i1,i2,i3) = (2*fcx(i1)*dvz(i1,i2,i3) 
     &              + fdx(i1)*((hy(i1,i2,i3)-hy(i1-1,i2,i3))/dx 
     &              - (hx(i1,i2,i3)-hx(i1,i2-1,i3))/dy-jz))/(1+fcx(i1))
               ez(i1,i2,i3) = 2*(fcy(i2)*ez(i1,i2,i3) 
     &              + fey(i2)*(gbz(i3)*dvz(i1,i2,i3)-bzm(i3)*dstore)
     &              /eps(i1,i2,i3))/(1+fcy(i2))
            enddo
         enddo
      enddo


      call PIC_fex(ex)
      call PIC_fey(ex)
      call PIC_fez(ex)
      call PIC_fex(dvx)
      call PIC_fey(dvx)
      call PIC_fez(dvx)
      call PIC_fex(ey)
      call PIC_fey(ey)
      call PIC_fez(ey)
      call PIC_fex(dvy)
      call PIC_fey(dvy)
      call PIC_fez(dvy)
      call PIC_fex(ez)
      call PIC_fey(ez)
      call PIC_fez(ez)
      call PIC_fex(dvz)
      call PIC_fey(dvz)
      call PIC_fez(dvz)


c energy conservation


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


      end subroutine PIC_pml_msb
