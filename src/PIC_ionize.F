c THIS SUBROUTINE IS THE IONIZER!

c korrigierte Version


      subroutine PIC_ionize

      use PIC_variables
      use VLA_variables
      use fields

      implicit none
      integer :: m,le,la
      integer :: j1,j2,j3
      integer :: l1,l2,l3
      integer :: lph,lph_n

      character*(5) :: node

      real(kind=8) :: dxi,dyi,dzi
      real(kind=8) :: pxi,pyi,pzi
      real(kind=8) :: xi,yi,zi
      real(kind=8) :: vxi,vyi,vzi
      real(kind=8) :: qni,mni,cni,wni,lni
      real(kind=8) :: xl,yl,zl
      real(kind=8) :: hmx,h0x,h1x,hmy,h0y,h1y,hmz,h0z,h1z
      real(kind=8) :: gmx,g0x,g1x,gmy,g0y,g1y,gmz,g0z,g1z

      real(kind=8) :: eeq,root
      real(kind=8) :: wad,ran
      real(kind=8) :: ebn,nef,cnl,wh1,wh2,wh3,wh4,wh5
      real(kind=8) :: h1,h2,h3
      real(kind=8) :: exq,eyq,ezq
      real(kind=8) :: u,v,w

      real(kind=8) :: s_cpud
      real(kind=8) :: s_cpuh

      real(kind=8),allocatable,dimension(:) :: eb0
      real(kind=8),allocatable,dimension(:) :: rndmv

      allocate(eb0(1:100))
      allocate(rndmv(1:npe))


      s_cpud=0.0d0
      s_cpuh=0.0d0
      call SERV_systime(cpug)


c ADK IONIZATION


      xl=0.5*dt
      yl=0.5*dt
      zl=0.5*dt
      dxi=1.0/dx
      dyi=1.0/dy
      dzi=1.0/dz
      eb0=1.0e10
      lph=niloc


      if (niloc.gt.0) then
         do l=1,niloc

            xi=p_niloc(11*l)
            yi=p_niloc(11*l+1)
            zi=p_niloc(11*l+2)
            pxi=p_niloc(11*l+3)
            pyi=p_niloc(11*l+4)
            pzi=p_niloc(11*l+5)
            qni=p_niloc(11*l+6)
            mni=p_niloc(11*l+7)
            cni=p_niloc(11*l+8)
            lni=p_niloc(11*l+9)
            wni=p_niloc(11*l+10)


            root=1.0/dsqrt(1.0+pxi*pxi+pyi*pyi+pzi*pzi)
            vxi=pxi*root
            vyi=pyi*root
            vzi=pzi*root

            xi=xi+vxi*xl
            yi=yi+vyi*yl
            zi=zi+vzi*zl


c FIELD INTERPOLATION


            u=xi*dxi
            v=yi*dyi
            w=zi*dzi
            j1=nint(u)
            j2=nint(v)
            j3=nint(w)
            h1=j1-u
            h2=j2-v
            h3=j3-w
            gmx=0.5*(0.5+h1)*(0.5+h1)
            gmy=0.5*(0.5+h2)*(0.5+h2)
            gmz=0.5*(0.5+h3)*(0.5+h3)
            g0x=0.75-h1*h1
            g0y=0.75-h2*h2
            g0z=0.75-h3*h3
            g1x=0.5*(0.5-h1)*(0.5-h1)
            g1y=0.5*(0.5-h2)*(0.5-h2)
            g1z=0.5*(0.5-h3)*(0.5-h3)

            u=xi*dxi-0.5
            v=yi*dyi-0.5
            w=zi*dzi-0.5
            l1=nint(u)
            l2=nint(v)
            l3=nint(w)
            h1=l1-u
            h2=l2-v
            h3=l3-w
            hmx=0.5*(0.5+h1)*(0.5+h1)
            hmy=0.5*(0.5+h2)*(0.5+h2)
            hmz=0.5*(0.5+h3)*(0.5+h3)
            h0x=0.75-h1*h1
            h0y=0.75-h2*h2
            h0z=0.75-h3*h3
            h1x=0.5*(0.5-h1)*(0.5-h1)
            h1y=0.5*(0.5-h2)*(0.5-h2)
            h1z=0.5*(0.5-h3)*(0.5-h3)

            exq=gmz*(gmy*(hmx*ex(l1-1,j2-1,j3-1)
     &                   +h0x*ex(l1,j2-1,j3-1)
     &                   +h1x*ex(l1+1,j2-1,j3-1))
     &              +g0y*(hmx*ex(l1-1,j2,j3-1)
     &                   +h0x*ex(l1,j2,j3-1)
     &                   +h1x*ex(l1+1,j2,j3-1))
     &              +g1y*(hmx*ex(l1-1,j2+1,j3-1)
     &                   +h0x*ex(l1,j2+1,j3-1)
     &                   +h1x*ex(l1+1,j2+1,j3-1)))
     &         +g0z*(gmy*(hmx*ex(l1-1,j2-1,j3)
     &                   +h0x*ex(l1,j2-1,j3)
     &                   +h1x*ex(l1+1,j2-1,j3))
     &              +g0y*(hmx*ex(l1-1,j2,j3)
     &                    +h0x*ex(l1,j2,j3)
     &                    +h1x*ex(l1+1,j2,j3))
     &              +g1y*(hmx*ex(l1-1,j2+1,j3)
     &                   +h0x*ex(l1,j2+1,j3)
     &                   +h1x*ex(l1+1,j2+1,j3)))
     &         +g1z*(gmy*(hmx*ex(l1-1,j2-1,j3+1)
     &                   +h0x*ex(l1,j2-1,j3+1)
     &                   +h1x*ex(l1+1,j2-1,j3+1))
     &              +g0y*(hmx*ex(l1-1,j2,j3+1)
     &                   +h0x*ex(l1,j2,j3+1)
     &                   +h1x*ex(l1+1,j2,j3+1))
     &              +g1y*(hmx*ex(l1-1,j2+1,j3+1)
     &                   +h0x*ex(l1,j2+1,j3+1)
     &                   +h1x*ex(l1+1,j2+1,j3+1)))

            eyq=gmz*(hmy*(gmx*ey(j1-1,l2-1,j3-1)
     &                   +g0x*ey(j1,l2-1,j3-1)
     &                   +g1x*ey(j1+1,l2-1,j3-1))
     &              +h0y*(gmx*ey(j1-1,l2,j3-1)
     &                   +g0x*ey(j1,l2,j3-1)
     &                   +g1x*ey(j1+1,l2,j3-1))
     &              +h1y*(gmx*ey(j1-1,l2+1,j3-1)
     &                   +g0x*ey(j1,l2+1,j3-1)
     &                   +g1x*ey(j1+1,l2+1,j3-1)))
     &         +g0z*(hmy*(gmx*ey(j1-1,l2-1,j3)
     &                   +g0x*ey(j1,l2-1,j3)
     &                   +g1x*ey(j1+1,l2-1,j3))
     &              +h0y*(gmx*ey(j1-1,l2,j3)
     &                   +g0x*ey(j1,l2,j3)
     &                   +g1x*ey(j1+1,l2,j3))
     &              +h1y*(gmx*ey(j1-1,l2+1,j3)
     &                   +g0x*ey(j1,l2+1,j3)
     &                   +g1x*ey(j1+1,l2+1,j3)))
     &         +g1z*(hmy*(gmx*ey(j1-1,l2-1,j3+1)
     &                   +g0x*ey(j1,l2-1,j3+1)
     &                   +g1x*ey(j1+1,l2-1,j3+1))
     &              +h0y*(gmx*ey(j1-1,l2,j3+1)
     &                   +g0x*ey(j1,l2,j3+1)
     &                   +g1x*ey(j1+1,l2,j3+1))
     &              +h1y*(gmx*ey(j1-1,l2+1,j3+1)
     &                   +g0x*ey(j1,l2+1,j3+1)
     &                   +g1x*ey(j1+1,l2+1,j3+1)))

            ezq=hmz*(gmy*(gmx*ez(j1-1,j2-1,l3-1)
     &                   +g0x*ez(j1,j2-1,l3-1)
     &                   +g1x*ez(j1+1,j2-1,l3-1))
     &              +g0y*(gmx*ez(j1-1,j2,l3-1)
     &                   +g0x*ez(j1,j2,l3-1)
     &                   +g1x*ez(j1+1,j2,l3-1))
     &              +g1y*(gmx*ez(j1-1,j2+1,l3-1)
     &                   +g0x*ez(j1,j2+1,l3-1)
     &                   +g1x*ez(j1+1,j2+1,l3-1)))
     &         +h0z*(gmy*(gmx*ez(j1-1,j2-1,l3)
     &                   +g0x*ez(j1,j2-1,l3)
     &                   +g1x*ez(j1+1,j2-1,l3))
     &              +g0y*(gmx*ez(j1-1,j2,l3)
     &                   +g0x*ez(j1,j2,l3)
     &                   +g1x*ez(j1+1,j2,l3))
     &              +g1y*(gmx*ez(j1-1,j2+1,l3)
     &                   +g0x*ez(j1,j2+1,l3)
     &                   +g1x*ez(j1+1,j2+1,l3)))
     &         +h1z*(gmy*(gmx*ez(j1-1,j2-1,l3+1)
     &                   +g0x*ez(j1,j2-1,l3+1)
     &                   +g1x*ez(j1+1,j2-1,l3+1))
     &              +g0y*(gmx*ez(j1-1,j2,l3+1)
     &                   +g0x*ez(j1,j2,l3+1)
     &                   +g1x*ez(j1+1,j2,l3+1))
     &              +g1y*(gmx*ez(j1-1,j2+1,l3+1)
     &                   +g0x*ez(j1,j2+1,l3+1)
     &                   +g1x*ez(j1+1,j2+1,l3+1)))


c            u=xq*dxi-0.5
c            v=yq*dyi
c            w=zq*dzi
c            i1=int(u)
c            i2=int(v)
c            i3=int(w)
c            if (u<0.0) i1=i1-1
c            if (v<0.0) i2=i2-1
c            if (w<0.0) i3=i3-1
c            u=u-i1
c            v=v-i2
c            w=w-i3

c            h1=(1.0-u)*(1.0-v)*(1.0-w)
c            h2=u*(1.0-v)*(1.0-w)
c            h3=(1.0-u)*v*(1.0-w)
c            h4=(1.0-u)*(1.0-v)*w
c            h5=u*v*(1.0-w)
c            h6=u*(1.0-v)*w
c            h7=(1.0-u)*v*w
c            h8=u*v*w

c            exq=h1*ex(i1,i2,i3)+h2*ex(i1+1,i2,i3)
c     &          +h3*ex(i1,i2+1,i3)+h4*ex(i1,i2,i3+1)
c     &          +h5*ex(i1+1,i2+1,i3)+h6*ex(i1+1,i2,i3+1)
c     &          +h7*ex(i1,i2+1,i3+1)+h8*ex(i1+1,i2+1,i3+1)

c            u=xq*dxi
c            v=yq*dyi-0.5
c            w=zq*dzi
c            i1=int(u)
c            i2=int(v)
c            i3=int(w)
c            if (u<0.0) i1=i1-1
c            if (v<0.0) i2=i2-1
c            if (w<0.0) i3=i3-1
c            u=u-i1
c            v=v-i2
c            w=w-i3

c            h1=(1.0-u)*(1.0-v)*(1.0-w)
c            h2=u*(1.0-v)*(1.0-w)
c            h3=(1.0-u)*v*(1.0-w)
c            h4=(1.0-u)*(1.0-v)*w
c            h5=u*v*(1.0-w)
c            h6=u*(1.0-v)*w
c            h7=(1.0-u)*v*w
c            h8=u*v*w

c            eyq=h1*ey(i1,i2,i3)+h2*ey(i1+1,i2,i3)
c     &          +h3*ey(i1,i2+1,i3)+h4*ey(i1,i2,i3+1)
c     &          +h5*ey(i1+1,i2+1,i3)+h6*ey(i1+1,i2,i3+1)
c     &          +h7*ey(i1,i2+1,i3+1)+h8*ey(i1+1,i2+1,i3+1)

c            u=xq*dxi
c            v=yq*dyi
c            w=zq*dzi-0.5
c            i1=int(u)
c            i2=int(v)
c            i3=int(w)
c            if (u<0.0) i1=i1-1
c            if (v<0.0) i2=i2-1
c            if (w<0.0) i3=i3-1
c            u=u-i1
c            v=v-i2
c            w=w-i3

c            h1=(1.0-u)*(1.0-v)*(1.0-w)
c            h2=u*(1.0-v)*(1.0-w)
c            h3=(1.0-u)*v*(1.0-w)
c            h4=(1.0-u)*(1.0-v)*w
c            h5=u*v*(1.0-w)
c            h6=u*(1.0-v)*w
c            h7=(1.0-u)*v*w
c            h8=u*v*w

c            ezq=h1*ez(i1,i2,i3)+h2*ez(i1+1,i2,i3)
c     &          +h3*ez(i1,i2+1,i3)+h4*ez(i1,i2,i3+1)
c     &          +h5*ez(i1+1,i2+1,i3)+h6*ez(i1+1,i2,i3+1)
c     &          +h7*ez(i1,i2+1,i3+1)+h8*ez(i1+1,i2+1,i3+1)


c IONIZATION OF HYDROGEN


            if (1835.0.le.mni.and.mni.le.1837.0) then

               eb0(1)=13.598/27.21        ! binding energies of H in au

               do m=1,1
                  if (qni.eq.+1.0*(m-1)) then

                     eeq=e0*sqrt(exq*exq+eyq*eyq+ezq*ezq)
     &                   /5.1422e11+1.0e-10

                     ebn=eb0(m)
                     nef=1.0/sqrt(2.0*ebn)
                     cnl=(1.0/sqrt(6.2831*nef))*(5.43656/nef)**nef
                     wh1=(2.0*ebn)**1.5
                     wh2=sqrt(3.0/(3.141562*wh1))
                     wh3=(2.0*wh1)**(2.0*nef-1.0)
                     wh4=dexp(-2.0*wh1/3.0)
                     wh5=4.1341e16*(dt/wl)*cnl*cnl*ebn

                     wad=wh2*wh5*sqrt(eeq)*(wh3/eeq**(2.0*nef-1.0))
     &                   *wh4**(1.0/eeq)
                     wad=1.0-exp(-wad)

                     call random_number(rndmv)
                     ran=rndmv(mpe+1)
 
                     if (wad.ge.ran) then
                        lph_n=lph+1
                        call REALLOC_particles(lph_n)
                        lph=lph_n

                        p_niloc(11*l+6)=+1.0*m
                        p_niloc(11*l+7)=mni
                        p_niloc(11*lph)=xi
                        p_niloc(11*lph+1)=yi
                        p_niloc(11*lph+2)=zi
                        p_niloc(11*lph+3)=pxi
                        p_niloc(11*lph+4)=pyi
                        p_niloc(11*lph+5)=pzi
                        p_niloc(11*lph+6)=-1.0
                        p_niloc(11*lph+7)=+1.0
                        p_niloc(11*lph+8)=cni
                        p_niloc(11*lph+9)=lph
                        p_niloc(11*lph+10)=wni
                     endif

                  endif
               enddo

            endif


c IONIZATION OF HELIUM


            if (7343.0.le.mni.and.mni.le.7345.0) then

               eb0(1)=24.587/27.21        ! binding energies of He in au
               eb0(2)=54.416/27.21        ! binding energies of He in au

               do m=1,2
                  if (qni.eq.+1.0*(m-1)) then

                     eeq=e0*sqrt(exq*exq+eyq*eyq+ezq*ezq)
     &                   /5.1422e11+1.0e-10

                     ebn=eb0(m)
                     nef=1.0/sqrt(2.0*ebn)
                     cnl=(1.0/sqrt(6.2831*nef))*(5.43656/nef)**nef
                     wh1=(2.0*ebn)**1.5
                     wh2=sqrt(3.0/(3.141562*wh1))
                     wh3=(2.0*wh1)**(2.0*nef-1.0)
                     wh4=dexp(-2.0*wh1/3.0)
                     wh5=4.1341e16*(dt/wl)*cnl*cnl*ebn

                     wad=wh2*wh5*sqrt(eeq)*(wh3/eeq**(2.0*nef-1.0))
     &                   *wh4**(1.0/eeq)
                     wad=1.0-exp(-wad)

                     call random_number(rndmv)
                     ran=rndmv(mpe+1)
 
                     if (wad.ge.ran) then
                        lph_n=lph+1
                        call REALLOC_particles(lph_n)
                        lph=lph_n

                        p_niloc(11*l+6)=+1.0*m
                        p_niloc(11*l+7)=mni
                        p_niloc(11*lph)=xi
                        p_niloc(11*lph+1)=yi
                        p_niloc(11*lph+2)=zi
                        p_niloc(11*lph+3)=pxi
                        p_niloc(11*lph+4)=pyi
                        p_niloc(11*lph+5)=pzi
                        p_niloc(11*lph+6)=-1.0
                        p_niloc(11*lph+7)=+1.0
                        p_niloc(11*lph+8)=cni
                        p_niloc(11*lph+9)=lph
                        p_niloc(11*lph+10)=wni
                     endif

                  endif
               enddo

            endif


c IONIZATION OF CARBON


            if (22031.0.le.mni.and.mni.le.22033.0) then

               eb0(1)=11.260/27.21        ! binding energies of C in au
               eb0(2)=24.380/27.21        ! binding energies of C in au
               eb0(3)=47.890/27.21        ! binding energies of C in au
               eb0(4)=64.500/27.21        ! binding energies of C in au
               eb0(5)=392.090/27.21       ! binding energies of C in au
               eb0(6)=490.000/27.21       ! binding energies of C in au

               do m=1,6
                  if (qni.eq.+1.0*(m-1)) then

                     eeq=e0*sqrt(exq*exq+eyq*eyq+ezq*ezq)
     &                   /5.1422e11+1.0e-10

                     ebn=eb0(m)
                     nef=1.0/sqrt(2.0*ebn)
                     cnl=(1.0/sqrt(6.2831*nef))*(5.43656/nef)**nef
                     wh1=(2.0*ebn)**1.5
                     wh2=sqrt(3.0/(3.141562*wh1))
                     wh3=(2.0*wh1)**(2.0*nef-1.0)
                     wh4=dexp(-2.0*wh1/3.0)
                     wh5=4.1341e16*(dt/wl)*cnl*cnl*ebn

                     wad=wh2*wh5*sqrt(eeq)*(wh3/eeq**(2.0*nef-1.0))
     &                   *wh4**(1.0/eeq)
                     wad=1.0-exp(-wad)

                     call random_number(rndmv)
                     ran=rndmv(mpe+1)

                     if (wad.ge.ran) then
                        lph_n=lph+1
                        call REALLOC_particles(lph_n)
                        lph=lph_n

                        p_niloc(11*l+6)=+1.0*m
                        p_niloc(11*l+7)=mni
                        p_niloc(11*lph)=xi
                        p_niloc(11*lph+1)=yi
                        p_niloc(11*lph+2)=zi
                        p_niloc(11*lph+3)=pxi
                        p_niloc(11*lph+4)=pyi
                        p_niloc(11*lph+5)=pzi
                        p_niloc(11*lph+6)=-1.0
                        p_niloc(11*lph+7)=+1.0
                        p_niloc(11*lph+8)=cni
                        p_niloc(11*lph+9)=lph
                        p_niloc(11*lph+10)=wni

                     endif

                  endif
               enddo

            endif


c IONIZATION OF OXYGEN


            if (29375.0.le.mni.and.mni.le.29377.0) then

               eb0(1)=13.618/27.21        ! binding energies of O in au
               eb0(2)=35.118/27.21        ! binding energies of O in au
               eb0(3)=54.936/27.21        ! binding energies of O in au
               eb0(4)=77.414/27.21        ! binding energies of O in au
               eb0(5)=113.900/27.21       ! binding energies of O in au
               eb0(6)=138.120/27.21       ! binding energies of O in au
               eb0(7)=739.300/27.21       ! binding energies of O in au
               eb0(8)=871.420/27.21       ! binding energies of O in au

               do m=1,8
                  if (qni.eq.+1.0*(m-1)) then

                     eeq=e0*sqrt(exq*exq+eyq*eyq+ezq*ezq)
     &                   /5.1422e11+1.0e-10

                     ebn=eb0(m)
                     nef=1.0/sqrt(2.0*ebn)
                     cnl=(1.0/sqrt(6.2831*nef))*(5.43656/nef)**nef
                     wh1=(2.0*ebn)**1.5
                     wh2=sqrt(3.0/(3.141562*wh1))
                     wh3=(2.0*wh1)**(2.0*nef-1.0)
                     wh4=dexp(-2.0*wh1/3.0)
                     wh5=4.1341e16*(dt/wl)*cnl*cnl*ebn

                     wad=wh2*wh5*sqrt(eeq)*(wh3/eeq**(2.0*nef-1.0))
     &                   *wh4**(1.0/eeq)
                     wad=1.0-exp(-wad)

                     call random_number(rndmv)
                     ran=rndmv(mpe+1)

                     if (wad.ge.ran) then
                        lph_n=lph+1
                        call REALLOC_particles(lph_n)
                        lph=lph_n

                        p_niloc(11*l+6)=+1.0*m
                        p_niloc(11*l+7)=mni
                        p_niloc(11*lph)=xi
                        p_niloc(11*lph+1)=yi
                        p_niloc(11*lph+2)=zi
                        p_niloc(11*lph+3)=pxi
                        p_niloc(11*lph+4)=pyi
                        p_niloc(11*lph+5)=pzi
                        p_niloc(11*lph+6)=-1.0
                        p_niloc(11*lph+7)=+1.0
                        p_niloc(11*lph+8)=cni
                        p_niloc(11*lph+9)=lph
                        p_niloc(11*lph+10)=wni
                     endif

                  endif
               enddo

            endif


c IONIZATION OF NEON


            if (36719.0.le.mni.and.mni.le.36721.0) then

               eb0(1)=21.56453/27.21        ! binding energies of Ne in au
               eb0(2)=40.96328/27.21        ! binding energies of Ne in au
               eb0(3)=63.45510/27.21        ! binding energies of Ne in au
               eb0(4)=97.11690/27.21        ! binding energies of Ne in au
               eb0(5)=126.21600/27.21       ! binding energies of Ne in au
               eb0(6)=157.93110/27.21       ! binding energies of Ne in au
               eb0(7)=207.27590/27.21       ! binding energies of Ne in au
               eb0(8)=239.09890/27.21       ! binding energies of Ne in au
               eb0(9)=1195.82200/27.21      ! binding energies of Ne in au
               eb0(10)=1362.19900/27.21     ! binding energies of Ne in au

               do m=1,10
                  if (qni.eq.+1.0*(m-1)) then

                     eeq=e0*sqrt(exq*exq+eyq*eyq+ezq*ezq)
     &                   /5.1422e11+1.0e-10

                     ebn=eb0(m)
                     nef=1.0/sqrt(2.0*ebn)
                     cnl=(1.0/sqrt(6.2831*nef))*(5.43656/nef)**nef
                     wh1=(2.0*ebn)**1.5
                     wh2=sqrt(3.0/(3.141562*wh1))
                     wh3=(2.0*wh1)**(2.0*nef-1.0)
                     wh4=dexp(-2.0*wh1/3.0)
                     wh5=4.1341e16*(dt/wl)*cnl*cnl*ebn

                     wad=wh2*wh5*sqrt(eeq)*(wh3/eeq**(2.0*nef-1.0))
     &                   *wh4**(1.0/eeq)
                     wad=1.0-exp(-wad)

                     call random_number(rndmv)
                     ran=rndmv(mpe+1)

                     if (wad.ge.ran) then
                        lph_n=lph+1
                        call REALLOC_particles(lph_n)
                        lph=lph_n

                        p_niloc(11*l+6)=+1.0*m
                        p_niloc(11*l+7)=mni
                        p_niloc(11*lph)=xi
                        p_niloc(11*lph+1)=yi
                        p_niloc(11*lph+2)=zi
                        p_niloc(11*lph+3)=pxi
                        p_niloc(11*lph+4)=pyi
                        p_niloc(11*lph+5)=pzi
                        p_niloc(11*lph+6)=-1.0
                        p_niloc(11*lph+7)=+1.0
                        p_niloc(11*lph+8)=cni
                        p_niloc(11*lph+9)=lph
                        p_niloc(11*lph+10)=wni

                     endif

                  endif
               enddo

            endif

c IONIZATION OF ARGON


            if (73439.0.le.mni.and.mni.le.73441.0) then

               eb0(1)=15.75996/27.21        ! binding energies of Ar in au
               eb0(2)=27.62917/27.21        ! binding energies of Ar in au
               eb0(3)=40.74208/27.21        ! binding energies of Ar in au
               eb0(4)=59.81241/27.21        ! binding energies of Ar in au
               eb0(5)=75.01684/27.21        ! binding energies of Ar in au
               eb0(6)=91.00897/27.21        ! binding energies of Ar in au
               eb0(7)=124.31984/27.21        ! binding energies of Ar in au
               eb0(8)=143.46271/27.21        ! binding energies of Ar in au
               

               do m=1,8
                  if (qni.eq.+1.0*(m-1)) then

                     eeq=e0*sqrt(exq*exq+eyq*eyq+ezq*ezq)
     &                   /5.1422e11+1.0e-10

                     ebn=eb0(m)
                     nef=1.0/sqrt(2.0*ebn)
                     cnl=(1.0/sqrt(6.2831*nef))*(5.43656/nef)**nef
                     wh1=(2.0*ebn)**1.5
                     wh2=sqrt(3.0/(3.141562*wh1))
                     wh3=(2.0*wh1)**(2.0*nef-1.0)
                     wh4=dexp(-2.0*wh1/3.0)
                     wh5=4.1341e16*(dt/wl)*cnl*cnl*ebn

                     wad=wh2*wh5*sqrt(eeq)*(wh3/eeq**(2.0*nef-1.0))
     &                   *wh4**(1.0/eeq)
                     wad=1.0-exp(-wad)

                     call random_number(rndmv)
                     ran=rndmv(mpe+1)

                     if (wad.ge.ran) then
                        lph_n=lph+1
                        call REALLOC_particles(lph_n)
                        lph=lph_n

                        p_niloc(11*l+6)=+1.0*m
                        p_niloc(11*l+7)=mni
                        p_niloc(11*lph)=xi
                        p_niloc(11*lph+1)=yi
                        p_niloc(11*lph+2)=zi
                        p_niloc(11*lph+3)=pxi
                        p_niloc(11*lph+4)=pyi
                        p_niloc(11*lph+5)=pzi
                        p_niloc(11*lph+6)=-1.0
                        p_niloc(11*lph+7)=+1.0
                        p_niloc(11*lph+8)=cni
                        p_niloc(11*lph+9)=lph
                        p_niloc(11*lph+10)=wni
                     endif

                  endif
               enddo

            endif


c IONIZATION OF NITROGEN


            if (25703.0.le.mni.and.mni.le.25705.0) then

               eb0(1)=14.531/27.21        ! binding energies of Ar in au
               eb0(2)=29.600/27.21        ! binding energies of Ar in au
               eb0(3)=47.449/27.21        ! binding energies of Ar in au
               eb0(4)=77.4732/27.21        ! binding energies of Ar in au
               eb0(5)=97.890/27.21        ! binding energies of Ar in au
               eb0(6)=522.071/27.21        ! binding energies of Ar in au
               eb0(7)=667.047/27.21        ! binding energies of Ar in au
               

               do m=1,7
                  if (qni.eq.+1.0*(m-1)) then

                     eeq=e0*sqrt(exq*exq+eyq*eyq+ezq*ezq)
     &                   /5.1422e11+1.0e-10

                     ebn=eb0(m)
                     nef=1.0/sqrt(2.0*ebn)
                     cnl=(1.0/sqrt(6.2831*nef))*(5.43656/nef)**nef
                     wh1=(2.0*ebn)**1.5
                     wh2=sqrt(3.0/(3.141562*wh1))
                     wh3=(2.0*wh1)**(2.0*nef-1.0)
                     wh4=dexp(-2.0*wh1/3.0)
                     wh5=4.1341e16*(dt/wl)*cnl*cnl*ebn

                     wad=wh2*wh5*sqrt(eeq)*(wh3/eeq**(2.0*nef-1.0))
     &                   *wh4**(1.0/eeq)
                     wad=1.0-exp(-wad)

                     call random_number(rndmv)
                     ran=rndmv(mpe+1)

                     if (wad.ge.ran) then
                        lph_n=lph+1
                        call REALLOC_particles(lph_n)
                        lph=lph_n

                        p_niloc(11*l+6)=+1.0*m
                        p_niloc(11*l+7)=mni
                        p_niloc(11*lph)=xi
                        p_niloc(11*lph+1)=yi
                        p_niloc(11*lph+2)=zi
                        p_niloc(11*lph+3)=pxi
                        p_niloc(11*lph+4)=pyi
                        p_niloc(11*lph+5)=pzi
                        p_niloc(11*lph+6)=-1.0
                        p_niloc(11*lph+7)=+1.0
                        p_niloc(11*lph+8)=cni
                        p_niloc(11*lph+9)=lph
                        p_niloc(11*lph+10)=wni
                     endif

                  endif
               enddo

            endif


c IONIZATION OF COPPER


            if (116669.0.le.mni.and.mni.le.116671.0) then

               eb0(1)=7.726381/27.21        ! binding energies of Cu in au
               eb0(2)=20.29240/27.21        ! binding energies of Cu in au
               eb0(3)=36.84000/27.21        ! binding energies of Cu in au
               eb0(4)=57.37000/27.21        ! binding energies of Cu in au
               eb0(5)=79.80000/27.21        ! binding energies of Cu in au
               eb0(6)=103.0000/27.21        ! binding energies of Cu in au
               eb0(7)=138.0000/27.21        ! binding energies of Cu in au
               eb0(8)=170.0000/27.21        ! binding energies of Cu in au
               eb0(9)=200.0000/27.21        ! binding energies of Cu in au
               eb0(10)=231.0000/27.21       ! binding energies of Cu in au
               eb0(11)=265.0000/27.21       ! binding energies of Cu in au
               eb0(12)=368.0000/27.21       ! binding energies of Cu in au
               eb0(13)=400.0000/27.21       ! binding energies of Cu in au
               eb0(14)=434.0000/27.21       ! binding energies of Cu in au
               eb0(15)=483.0000/27.21       ! binding energies of Cu in au
               eb0(16)=519.0000/27.21       ! binding energies of Cu in au
               eb0(17)=557.0000/27.21       ! binding energies of Cu in au
               eb0(18)=632.9000/27.21       ! binding energies of Cu in au
               eb0(19)=670.5880/27.21       ! binding energies of Cu in au
               eb0(20)=1689.0000/27.21      ! binding energies of Cu in au
               eb0(21)=1803.0000/27.21      ! binding energies of Cu in au
               eb0(22)=1915.0000/27.21      ! binding energies of Cu in au
               eb0(23)=2060.0000/27.21      ! binding energies of Cu in au
               eb0(24)=2180.0000/27.21      ! binding energies of Cu in au
               eb0(25)=2308.0000/27.21      ! binding energies of Cu in au
               eb0(26)=2477.9000/27.21      ! binding energies of Cu in au
               eb0(27)=2587.0000/27.21      ! binding energies of Cu in au
               eb0(28)=11062.0000/27.21     ! binding energies of Cu in au
               eb0(29)=11567.6100/27.21     ! binding energies of Cu in au

               do m=1,29
                  if (qni.eq.+1.0*(m-1)) then

                     eeq=e0*sqrt(exq*exq+eyq*eyq+ezq*ezq)
     &                   /5.1422e11+1.0e-10

                     ebn=eb0(m)
                     nef=1.0/sqrt(2.0*ebn)
                     cnl=(1.0/sqrt(6.2831*nef))*(5.43656/nef)**nef
                     wh1=(2.0*ebn)**1.5
                     wh2=sqrt(3.0/(3.141562*wh1))
                     wh3=(2.0*wh1)**(2.0*nef-1.0)
                     wh4=dexp(-2.0*wh1/3.0)
                     wh5=4.1341e16*(dt/wl)*cnl*cnl*ebn

                     wad=wh2*wh5*sqrt(eeq)*(wh3/eeq**(2.0*nef-1.0))
     &                   *wh4**(1.0/eeq)
                     wad=1.0-exp(-wad)

                     call random_number(rndmv)
                     ran=rndmv(mpe+1)

                     if (wad.ge.ran) then
                        lph_n=lph+1
                        call REALLOC_particles(lph_n)
                        lph=lph_n

                        p_niloc(11*l+6)=+1.0*m
                        p_niloc(11*l+7)=mni
                        p_niloc(11*lph)=xi
                        p_niloc(11*lph+1)=yi
                        p_niloc(11*lph+2)=zi
                        p_niloc(11*lph+3)=pxi
                        p_niloc(11*lph+4)=pyi
                        p_niloc(11*lph+5)=pzi
                        p_niloc(11*lph+6)=-1.0
                        p_niloc(11*lph+7)=+1.0
                        p_niloc(11*lph+8)=cni
                        p_niloc(11*lph+9)=lph
                        p_niloc(11*lph+10)=wni
                     endif

                  endif
               enddo

            endif

         enddo
      endif

      niloc=lph

      deallocate(eb0)
      deallocate(rndmv)


      call SERV_systime(cpuh)
      s_cpuh=s_cpuh+cpuh-cpug


      cpuinou=cpuinou+s_cpud
      cpucomp=cpucomp+s_cpuh-s_cpud


      end subroutine PIC_ionize
