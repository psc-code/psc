c BEAM INJECTOR ALONG z-AXIS. PARTICLES ARE EMITTED NORMAL TO A SPERICAL SHELL.
c THEY MAY HAVE TEMPERATURE ALONG THE SURFACE NORMAL.


      subroutine PIC_bpulse

      use PSC_globals
      use PSC_patch

      implicit none
      include './mpif.h'

      integer :: m,cell,le,la
      integer :: nqte,nqti
      integer :: rds1,rds2,rds3

      character*(5) :: node

      real(kind=8) :: tb,dtb,tbenv
      real(kind=8) :: xmin,xmax,ymin,ymax,zmin,zmax
      real(kind=8) :: neb,xe,ye,ze,re,de,ae
      real(kind=8) :: nib,xi,yi,zi,ri,di,ai
      real(kind=8) :: phi,theta,ph,th,ang
      real(kind=8) :: pp,lb
      real(kind=8) :: nqe,nqi,nqet,nqit
      real(kind=8) :: nbe,nbi,snbe,snbi,nbet,nbit
      real(kind=8) :: Ibe,Ibi,sIbe,sIbi,Ibet,Ibit
      real(kind=8) :: Ixbe,Iybe,Izbe,Ixbi,Iybi,Izbi
      real(kind=8) :: sIxbe,sIybe,sIzbe,sIxbi,sIybi,sIzbi
      real(kind=8) :: ran1,ran2,ran3,ran4,proe,proi
      real(kind=8) :: x,y,z,px,py,pz
      real(kind=8) :: qni,tni,pni,mni,cni,lni,wni
      real(kind=8) :: qne,tne,pne,mne,cne,lne,wne
      real(kind=8) :: nl,nx,ny,nz,u,v,w,vx,vy,vz,rr

      real(kind=8) :: s_cpub
      real(kind=8) :: s_cpud
      real(kind=8) :: s_cpuh

      real(kind=8),allocatable,dimension(:) :: nte,nti,ntn
      real(kind=8),allocatable,dimension(:) :: beam
      real(kind=8),allocatable,dimension(:,:) :: beam_ary
      real(kind=8),allocatable,dimension(:) :: rndmv


      allocate(nte(i3mn:i3mx))
      allocate(nti(i3mn:i3mx))
      allocate(ntn(i3mn:i3mx))

      allocate(beam(1:6))
      allocate(beam_ary(1:6,0:npe-1))


      s_cpub=0.0d0
      s_cpud=0.0d0
      s_cpuh=0.0d0
      call SERV_systime(cpug)

 
c BOUNDING BOX


      xmin=(i1mn-0.5)*dx
      xmax=(i1mx+0.5)*dx
      ymin=(i2mn-0.5)*dy
      ymax=(i2mx+0.5)*dy
      zmin=(i3mn-0.5)*dz
      zmax=(i3mx+0.5)*dz


 
c  -rd <= i1 <= i1mx+rd-1
c  -rd <= i2 <= i2mx+rd-1
c  -rd <= i3 <= i3mx+rd-1


       rds1=rd1
       rds2=rd2
       rds3=rd3

       if (i1n==i1x) rds1=0
       if (i2n==i2x) rds2=0
       if (i3n==i3x) rds3=0


c TEMPORAL ENVELOPE OF PARTICLE BEAM


c dtb: temporal width of beam in sec
c tb: time when peak beam density is reached in sec


      dtb=2.0e-6
      tb=2.0e-6

      dtb=dtb*wl
      tb=tb*wl

      if (n*dt.lt.tb) then
         tbenv=1.0
      else if (n*dt.ge.tb) then
         tbenv=0.0
      endif


c SHAPE, LOCATION, AND TOTAL PARTICLE NUMBER OF SPHERICAL EMITTER 

c |                           
c |                           
c |  .               re
c |            .              
c |  de                  .    
c |               theta             .      |  
c -----------------------------------------x-------->
c |                                        |       z
c |                                        ze
c |
c |
c |

c neb: average beam density in m^3
c wne: weight of quasi-particle
c xe: location of shell center in x in m
c ye: location of shell center in y in m
c ze: location of shell center in z in m
c re: radius of emitter surface in m
c de: radius of emitter disk in m


      neb=2.0e10*tbenv
      wne=0.05
      xe=50.0/ld
      ye=50.0/ld
      ze=-50000.0/ld
      re=50001.0/ld
      de=0.0/ld


c |                           
c |                           
c |  .               ri
c |            .              
c |  di                  .    
c |               theta             .      |  
c -----------------------------------------x-------->
c |                                        |       z
c |                                        zi
c |
c |
c |

c nib: average beam density in m^3
c wni: weight of quasi-particle
c xi: location of shell center in x in m
c yi: location of shell center in y in m
c zi: location of shell center in z in m
c ri: radius of emitter surface in m
c di: radius of emitter disk in m


      nib=2.0e10*tbenv
      wni=0.05
      xi=50.0/ld
      yi=50.0/ld
      zi=-50000.0/ld
      ri=50001.0/ld
      di=0.0/ld


c: nqe: number of negative quasi-particles 
c: nqi: number of positive quasi-particles 


      nqe=0.0
      nqi=0.0


c: snbe: number of negative real particles 
c: sIxbe: negative current along x-axis in A 
c: sIybe: negative current along y-axis in A 
c: sIzbe: negative current along z-axis in A 
c: snbi: number of positive real particles 
c: sIxbi: positive current along x-axis in A 
c: sIybi: positive current along y-axis in A 
c: sIzbi: positive current along z-axis in A 


      snbe=0.0
      sIxbe=0.0
      sIybe=0.0
      sIzbe=0.0
      sIbe=0.0

      snbi=0.0
      sIxbi=0.0
      sIybi=0.0
      sIzbi=0.0
      sIbi=0.0

      beam=0.0
      beam_ary=0.0


c BEAM INJECTOR


      if (i1x-i1n.gt.0.and.i2x-i2n.gt.0) then
         ai=6.2831853072*ri*ri*ld*ld*(1.0-sqrt(1.0-(di/ri)*(di/ri)))      ! area of spherical shell
         nib=nib*ai*dz*ld                                                 ! volume of spherical shell
         nqti=nint(nib/(wni*n0*cori*dx*dy*dz*ld*ld*ld))                   ! number of quasi-particles
      elseif (i1x-i1n.gt.0) then
         ai=2.0*di*dy*ld*ld                                               ! area of a line along x
         nib=nib*ai*dz*ld                                                 ! volume of a line along x
         nqti=nint(nib/(wni*n0*cori*dx*dy*dz*ld*ld*ld))                   ! number of quasi-particles
      elseif (i2x-i2n.gt.0) then
         ai=2.0*di*dx*ld*ld
         nib=nib*ai*dz*ld
         nqti=nint(nib/(wni*n0*cori*dx*dy*dz*ld*ld*ld))
      else
         ai=dx*dy*ld*ld
         nib=nib*ai*dz*ld
         nqti=nint(nib/(wni*n0*cori*dx*dy*dz*ld*ld*ld))
      endif


      if (nqti.gt.0) then

         allocate(rndmv(1:6*nqti))
         call random_number(rndmv)

         do m=1,nqti

            ran1=max(1.0d-20,rndmv(6*m-5))
            ran2=rndmv(6*m-4)

            if (i1x-i1n.gt.0.and.i2x-i2n.gt.0) then
               lb=cos(asin(de/re))             ! homogenious distribution
               theta=acos(1.0+(lb-1.0)*ran1)   ! homogenious distribution
               phi=6.2831853*ran2              ! homogenious distribution
               x=xe+re*sin(theta)*cos(phi)
               y=ye+re*sin(theta)*sin(phi)
               z=ze+re*cos(theta)
            elseif (i1x-i1n.gt.0) then
               lb=cos(asin(de/re))             ! homogenious distribution
               theta=acos(1.0+(lb-1.0)*ran1)   ! homogenious distribution
               ye=i2n*dy
               x=xe+re*sin(theta)
               y=ye
               z=ze+re*cos(theta)
            elseif (i2x-i2n.gt.0) then
               lb=cos(asin(de/re))             ! homogenious distribution
               theta=acos(1.0+(lb-1.0)*ran1)   ! homogenious distribution
               xe=i1n*dx
               x=xe
               y=ye+re*sin(theta)
               z=ze+re*cos(theta)
            else
               xe=i1n*dx
               ye=i2n*dy
               x=xe
               y=ye
               z=ze+re
            endif

            if (i1x-i1n.gt.0.and.i2x-i2n.gt.0) then
               lb=cos(asin(di/ri))             ! homogenious distribution
               theta=acos(1.0+(lb-1.0)*ran1)   ! homogenious distribution
               phi=6.2831853*ran2              ! homogenious distribution
               x=xi+ri*sin(theta)*cos(phi)
               y=yi+ri*sin(theta)*sin(phi)
               z=zi+ri*cos(theta)
            elseif (i1x-i1n.gt.0) then
               lb=cos(asin(di/ri))             ! homogenious distribution
               theta=acos(1.0+(lb-1.0)*ran1)   ! homogenious distribution
               yi=i2n*dy
               x=xi+ri*sin(theta)
               y=yi
               z=zi+ri*cos(theta)
            elseif (i2x-i2n.gt.0) then
               lb=cos(asin(di/ri))             ! homogenious distribution
               theta=acos(1.0+(lb-1.0)*ran1)   ! homogenious distribution
               xi=i1n*dx
               x=xi
               y=yi+ri*sin(theta)
               z=zi+ri*cos(theta)
            else
               xi=i1n*dx
               yi=i2n*dy
               x=xi
               y=yi
               z=zi+ri
            endif

            if (zmin<=z.and.z<zmax) then
               if (ymin<=y.and.y<ymax) then
                  if (xmin<=x.and.x<xmax) then


                     nl=sqrt((x-xi)*(x-xi)+(y-yi)*(y-yi)+(z-zi)*(z-zi))
                     nx=(x-xi)/nl
                     ny=(y-yi)/nl
                     nz=(z-zi)/nl

                     ran1=max(1.0d-20,rndmv(6*m-3))
                     ran2=rndmv(6*m-2)
                     ran3=rndmv(6*m-1)
                     ran4=rndmv(6*m)

                     pni=+1.732                      ! particles i
                     ang=+0.0
                     mni=+1836.0
                     tni=+0.0


c Particles are emitted normal to a spherical surface such that the prescibed
c density at the emitting surface remains constant. The emitted current is
c determined by the prescribed particle spectrum.


                     pp=pni+sqrt((1.0-tni*beta**2*log(ran1)/mni)**2-1.0)
                     ang=cos(ang*3.1415926536/180.0)
                     th=acos(1.0+(ang-1.0)*ran2)
                     ph=6.2831853*ran3
                     nl=sqrt(nx*nx+ny*ny+1.0e-20)
                     px=pp*(cos(th)*nx
     &                      +sin(th)*sin(ph)*ny/nl
     &                      -sin(th)*cos(ph)*nx*nz/nl)
                     py=pp*(cos(th)*ny
     &                      -sin(th)*sin(ph)*nx/nl
     &                      -sin(th)*cos(ph)*ny*nz/nl)
                     pz=pp*(cos(th)*nz
     &                      +sin(th)*cos(ph)*nx*nx/nl
     &                      +sin(th)*cos(ph)*ny*ny/nl)

                     rr=sqrt(1.0+px*px+py*py+pz*pz)
                     vx=px/rr
                     vy=py/rr
                     vz=pz/rr

                     proi=1.0-(1.0-abs(vx*dt/dx))
     &                       *(1.0-abs(vy*dt/dy))
     &                       *(1.0-abs(vz*dt/dz))


c Depending on the velocity of the particle to be injected it is determined
c how big the average contribution of this particle to the density at the
c emitting surface is. If the contribution is low the particle is placed
c with high probability.


                     if (ran4<proi) then
                        niloc_n=niloc+1
                        call REALLOC_particles(niloc_n)
                        niloc=niloc_n


                        u=x/dx
                        v=y/dy
                        w=z/dz
                        i1=int(u)
                        i2=int(v)
                        i3=int(w)
                        if (u<0.0) i1=i1-1
                        if (v<0.0) i2=i2-1
                        if (w<0.0) i3=i3-1


                        cell=(i1-i1mn+rds1+1)
     &                        +(i1mx-i1mn+2*rds1+1)
     &                        *(i2-i2mn+rds2)
     &                        +(i1mx-i1mn+2*rds1+1)
     &                        *(i2mx-i2mn+2*rds2+1)
     &                        *(i3-i3mn+rds3)


                        qni=+1.0
                        cni=+cell
                        lni=+1


                        p_niloc(11*niloc+0)=x
                        p_niloc(11*niloc+1)=y
                        p_niloc(11*niloc+2)=z
                        p_niloc(11*niloc+3)=px
                        p_niloc(11*niloc+4)=py
                        p_niloc(11*niloc+5)=pz
                        p_niloc(11*niloc+6)=qni
                        p_niloc(11*niloc+7)=mni
                        p_niloc(11*niloc+8)=cni
                        p_niloc(11*niloc+9)=lni
                        p_niloc(11*niloc+10)=wni

                        nqi=nqi+1

                        nbi=n0*cori*dx*dy*dz*ld*ld*ld
                        Ixbi=qni*alpha*alpha*cori
     &                       *vx*dy*dz*ld*ld*j0/eta
                        Iybi=qni*alpha*alpha*cori
     &                       *vy*dx*dz*ld*ld*j0/eta
                        Izbi=qni*alpha*alpha*cori
     &                       *vz*dx*dy*ld*ld*j0/eta
                        Ibi=Ixbi*nx+Iybi*ny+Izbi*nz

                        snbi=snbi+nbi
                        sIxbi=sIxbi+Ixbi
                        sIybi=sIybi+Iybi
                        sIzbi=sIzbi+Izbi
                        sIbi=sIbi+Ibi

                     endif

                  endif
               endif
            endif

         enddo

         deallocate(rndmv)

      endif



      if (i1x-i1n.gt.0.and.i2x-i2n.gt.0) then
         ae=6.2831853072*re*re*ld*ld*(1.0-sqrt(1.0-(de/re)*(de/re)))      ! area of spherical shell
         neb=neb*ae*dz*ld                                                 ! volume of spherical shell
         nqte=nint(neb/(wne*n0*cori*dx*dy*dz*ld*ld*ld))                   ! number of quasi-particles
      elseif (i1x-i1n.gt.0) then
         ae=2.0*de*dy*ld*ld                                               ! area of a line along x
         neb=neb*ae*dz*ld                                                 ! volume of a line along x
         nqte=nint(neb/(wne*n0*cori*dx*dy*dz*ld*ld*ld))                   ! number of quasi-particles
      elseif (i2x-i2n.gt.0) then
         ae=2.0*de*dx*ld*ld
         neb=neb*ae*dz*ld
         nqte=nint(neb/(wne*n0*cori*dx*dy*dz*ld*ld*ld))
      else
         ae=dx*dy*ld*ld
         neb=neb*ae*dz*ld
         nqte=nint(neb/(wne*n0*cori*dx*dy*dz*ld*ld*ld))
      endif


      if (nqte.gt.0) then

         allocate(rndmv(1:6*nqte))
         call random_number(rndmv)

         do m=1,nqte

            ran1=rndmv(6*m-5)
            ran2=rndmv(6*m-4)

            if (i1x-i1n.gt.0.and.i2x-i2n.gt.0) then
               lb=cos(asin(de/re))             ! homogenious distribution
               theta=acos(1.0+(lb-1.0)*ran1)   ! homogenious distribution
               phi=6.2831853*ran2              ! homogenious distribution
               x=xe+re*sin(theta)*cos(phi)
               y=ye+re*sin(theta)*sin(phi)
               z=ze+re*cos(theta)
            elseif (i1x-i1n.gt.0) then
               lb=cos(asin(de/re))             ! homogenious distribution
               theta=acos(1.0+(lb-1.0)*ran1)   ! homogenious distribution
               ye=i2n*dy
               x=xe+re*sin(theta)
               y=ye
               z=ze+re*cos(theta)
            elseif (i2x-i2n.gt.0) then
               lb=cos(asin(de/re))             ! homogenious distribution
               theta=acos(1.0+(lb-1.0)*ran1)   ! homogenious distribution
               xe=i1n*dx
               x=xe
               y=ye+re*sin(theta)
               z=ze+re*cos(theta)
            else
               xe=i1n*dx
               ye=i2n*dy
               x=xe
               y=ye
               z=ze+re
            endif

            if (zmin<=z.and.z<zmax) then
               if (ymin<=y.and.y<ymax) then
                  if (xmin<=x.and.x<xmax) then

                     nl=sqrt((x-xe)*(x-xe)+(y-ye)*(y-ye)+(z-ze)*(z-ze))
                     nx=(x-xe)/nl
                     ny=(y-ye)/nl
                     nz=(z-ze)/nl

                     ran1=max(1.0d-20,rndmv(6*m-3))
                     ran2=rndmv(6*m-2)
                     ran3=rndmv(6*m-1)
                     ran4=rndmv(6*m)

                     pne=+1.732                      ! particles e
                     ang=+0.0
                     mne=+1.0
                     tne=+0.2


c Particles are emitted normal to a spherical surface such that the prescibed
c density at the emitting surface remains constant. The emitted current is
c determined by the prescribed particle spectrum.


                     pp=pne+sqrt((1.0-tne*beta**2*log(ran1)/mne)**2-1.0)
                     ang=cos(ang*3.1415926536/180.0)
                     th=acos(1.0+(ang-1.0)*ran2)
                     ph=6.2831853*ran3
                     nl=sqrt(nx*nx+ny*ny+1.0e-20)
                     px=pp*(cos(th)*nx
     &                      +sin(th)*sin(ph)*ny/nl
     &                      -sin(th)*cos(ph)*nx*nz/nl)
                     py=pp*(cos(th)*ny
     &                      -sin(th)*sin(ph)*nx/nl
     &                      -sin(th)*cos(ph)*ny*nz/nl)
                     pz=pp*(cos(th)*nz
     &                      +sin(th)*cos(ph)*nx*nx/nl
     &                      +sin(th)*cos(ph)*ny*ny/nl)

                     rr=sqrt(1.0+px*px+py*py+pz*pz)
                     vx=px/rr
                     vy=py/rr
                     vz=pz/rr

                     proe=1.0-(1.0-abs(vx*dt/dx))
     &                       *(1.0-abs(vy*dt/dy))
     &                       *(1.0-abs(vz*dt/dz))


c Depending on the velocity of the particle to be injected it is determined
c how big the average contribution of this particle to the density at the
c emitting surface is. If the contribution is low the particle is placed
c with high probability.


                     if (ran4<proe) then
                        niloc_n=niloc+1
                        call REALLOC_particles(niloc_n)
                        niloc=niloc_n


                        u=x/dx
                        v=y/dy
                        w=z/dz
                        i1=int(u)
                        i2=int(v)
                        i3=int(w)
                        if (u<0.0) i1=i1-1
                        if (v<0.0) i2=i2-1
                        if (w<0.0) i3=i3-1


                        cell=(i1-i1mn+rds1+1)
     &                        +(i1mx-i1mn+2*rds1+1)
     &                        *(i2-i2mn+rds2)
     &                        +(i1mx-i1mn+2*rds1+1)
     &                        *(i2mx-i2mn+2*rds2+1)
     &                        *(i3-i3mn+rds3)


                        qne=-1.0
                        cne=+cell
                        lne=+1


                        p_niloc(11*niloc+0)=x
                        p_niloc(11*niloc+1)=y
                        p_niloc(11*niloc+2)=z
                        p_niloc(11*niloc+3)=px
                        p_niloc(11*niloc+4)=py
                        p_niloc(11*niloc+5)=pz
                        p_niloc(11*niloc+6)=qne
                        p_niloc(11*niloc+7)=mne
                        p_niloc(11*niloc+8)=cne
                        p_niloc(11*niloc+9)=lne
                        p_niloc(11*niloc+10)=wne

                        nqe=nqe+1

                        nbe=wne*n0*cori*dx*dy*dz*ld*ld*ld
                        Ixbe=wne*qne*alpha*alpha*cori
     &                       *vx*dy*dz*ld*ld*j0/eta
                        Iybe=wne*qne*alpha*alpha*cori
     &                       *vy*dx*dz*ld*ld*j0/eta
                        Izbe=wne*qne*alpha*alpha*cori
     &                       *vz*dx*dy*ld*ld*j0/eta
                        Ibe=Ixbe*nx+Iybe*ny+Izbe*nz

                        snbe=snbe+nbe
                        sIxbe=sIxbe+Ixbe
                        sIybe=sIybe+Iybe
                        sIzbe=sIzbe+Izbe
                        sIbe=sIbe+Ibe

                     endif

                  endif
               endif
            endif

         enddo

         deallocate(rndmv)

      endif


c CALCULATION OF BEAM DENSITY VIA AREA DETERMINATION


      if (nqti.gt.0.or.nqte.gt.0) then

         snbe=snbe/(ae*dz*ld)
         snbi=snbi/(ai*dz*ld)

         beam_ary(1,mpe)=nqe
         beam_ary(2,mpe)=snbe
         beam_ary(3,mpe)=sIbe
         beam_ary(4,mpe)=nqi
         beam_ary(5,mpe)=snbi
         beam_ary(6,mpe)=sIbi


         call SERV_systime(cpua)
         do pec=0,npe-1
            beam(1)=beam_ary(1,pec)
            beam(2)=beam_ary(2,pec)
            beam(3)=beam_ary(3,pec)
            beam(4)=beam_ary(4,pec)
            beam(5)=beam_ary(5,pec)
            beam(6)=beam_ary(6,pec)
            call MPI_BCAST(beam,6,MPI_DOUBLE_PRECISION,pec,
     &                     MPI_COMM_WORLD,info)
            beam_ary(1,pec)=beam(1)
            beam_ary(2,pec)=beam(2)
            beam_ary(3,pec)=beam(3)
            beam_ary(4,pec)=beam(4)
            beam_ary(5,pec)=beam(5)
            beam_ary(6,pec)=beam(6)
         enddo
         call SERV_systime(cpub)
         s_cpub=s_cpub+cpub-cpua


         nqet=0.0
         nbet=0.0
         Ibet=0.0
         nqit=0.0
         nbit=0.0
         Ibit=0.0


         do pec=0,npe-1
            nqet=nqet+beam_ary(1,pec)
            nbet=nbet+beam_ary(2,pec)
            Ibet=Ibet+beam_ary(3,pec)
            nqit=nqit+beam_ary(4,pec)
            nbit=nbit+beam_ary(5,pec)
            Ibit=Ibit+beam_ary(6,pec)
         enddo


         call SERV_systime(cpuc)
         if (mpe.eq.0) then
            write(6,*) ' '
            write(6,*) '!!OUTPUT FROM "PIC_bpulse.f"!!'
            write(6,*) 'timestep:',n
            write(6,*) 'neb: ',neb
            write(6,*) 'nqte: ',nqte
            write(6,*) 'nqet: ',nqet
            write(6,*) 'nbet in 1/m^3: ',nbet
            write(6,*) 'Ibet in A: ',Ibet
            write(6,*) 'nib: ',nib
            write(6,*) 'nqti: ',nqti
            write(6,*) 'nqit: ',nqit
            write(6,*) 'nbit in 1/m^3: ',nbit
            write(6,*) 'Ibit in A: ',Ibit
            write(6,*) ' '
         endif
         call SERV_systime(cpud)
         s_cpud=s_cpud+cpud-cpuc

      endif


      deallocate(beam,beam_ary)
      deallocate(nte,nti,ntn)

      call SERV_systime(cpuh)
      s_cpuh=s_cpuh+cpuh-cpug

      cpuinou=cpuinou+s_cpud
      cpucomp=cpucomp+s_cpuh-s_cpud-s_cpub


      end subroutine PIC_bpulse
