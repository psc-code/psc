
c     ----------------------------------------------------------------------
      subroutine PIC_push_part_xz(niloc, p_niloc, p2A, p2B,
     &     jxi, jyi, jzi, ex, ey, ez, bx, by, bz)
c     ----------------------------------------------------------------------
c     particle pusher yz, w/o the following boundary exchange

      use PSC_globals
      use PSC_patch

      implicit none

      integer      :: niloc
      real(kind=8) :: p_niloc(0:*)
      real(kind=8) :: p2A, p2B
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: jxi, jyi, jzi
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: ex, ey, ez
      real(kind=8),dimension(i1mn-rd1:i1mx+rd1,
     &                       i2mn-rd2:i2mx+rd2,
     &                       i3mn-rd3:i3mx+rd3) :: bx, by, bz

      integer :: j1,j2,j3,k1,k3,l1,l2,l3
      integer :: l1min,l1max,l3min,l3max
      integer :: l

      real(kind=8) :: dxi,dyi,dzi
      real(kind=8) :: pxi,pyi,pzi
      real(kind=8) :: pxm,pym,pzm,pxp,pyp,pzp
      real(kind=8) :: qni,mni,cni,lni,wni
      real(kind=8) :: xi,yi,zi,vxi,vyi,vzi,root
      real(kind=8) :: xl,yl,zl
      
      real(kind=8) :: dqs,fnqs,fnqxs,fnqys,fnqzs
      real(kind=8) :: dq,fnq,fnqx,fnqy,fnqz

      real(kind=8) :: h1,h2,h3
      real(kind=8) :: hmx,h0x,h1x,hmz,h0z,h1z
      real(kind=8) :: gmx,g0x,g1x,gmz,g0z,g1z
      real(kind=8) :: wx,wy,wz
      real(kind=8) :: exq,eyq,ezq
      real(kind=8) :: bxq,byq,bzq
      real(kind=8) :: taux,tauy,tauz,tau
      real(kind=8) :: u,v,w

      real(kind=8),allocatable,dimension(:) :: s0x,s0z


      real(kind=8),allocatable,dimension(:) :: s1x,s1z
      real(kind=8),allocatable,dimension(:,:,:) :: jxh,jyh,jzh

      allocate(s0x(-2:2))
      allocate(s0z(-2:2))
      allocate(s1x(-2:2))
      allocate(s1z(-2:2))

      allocate(jxh(-3:2,-2:2,-2:2))
      allocate(jyh(-2:2,-3:2,-2:2))
      allocate(jzh(-2:2,-2:2,-3:2))


c INITIALIZATION


      xl=0.5*dt
      yl=0.5*dt
      zl=0.5*dt
      dqs=0.5*eta*dt
      fnqs=alpha*alpha*cori/eta
      fnqxs=dx*fnqs/dt
      fnqys=dy*fnqs/dt
      fnqzs=dz*fnqs/dt
      dxi=1.0/dx
      dyi=1.0/dy
      dzi=1.0/dz


      jxi=0.0d0
      jyi=0.0d0
      jzi=0.0d0


      p2A=0.0d0
      p2B=0.0d0


c PARTICLE LOOP


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

c CHARGE DENSITY FORM FACTOR AT (n+0.5)*dt 
c x^n, p^n -> x^(n+0.5), p^n

            root=1.0/dsqrt(1.0+pxi*pxi+pyi*pyi+pzi*pzi)
            vxi=pxi*root
            vyi=pyi*root
            vzi=pzi*root

            p2A=p2A+mni*fnqs*(1.0d0/root-1.0d0)/eta

            xi=xi+vxi*xl
            zi=zi+vzi*zl

            s0x=0.0
            s0z=0.0
            s1x=0.0
            s1z=0.0

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
            gmz=0.5*(0.5+h3)*(0.5+h3)
            g0x=0.75-h1*h1
            g0z=0.75-h3*h3
            g1x=0.5*(0.5-h1)*(0.5-h1)
            g1z=0.5*(0.5-h3)*(0.5-h3)

            s0x(-1)=0.5*(1.5-abs(h1-1.0))*(1.5-abs(h1-1.0))
            s0x(+0)=0.75-abs(h1)*abs(h1)
            s0x(+1)=0.5*(1.5-abs(h1+1.0))*(1.5-abs(h1+1.0))
            s0z(-1)=0.5*(1.5-abs(h3-1.0))*(1.5-abs(h3-1.0))
            s0z(+0)=0.75-abs(h3)*abs(h3)
            s0z(+1)=0.5*(1.5-abs(h3+1.0))*(1.5-abs(h3+1.0))

            u=xi*dxi-0.5
            v=yi*dyi
            w=zi*dzi-0.5
            l1=nint(u)
            l2=nint(v)
            l3=nint(w)
            h1=l1-u
            h2=l2-v
            h3=l3-w
            hmx=0.5*(0.5+h1)*(0.5+h1)
            hmz=0.5*(0.5+h3)*(0.5+h3)
            h0x=0.75-h1*h1
            h0z=0.75-h3*h3
            h1x=0.5*(0.5-h1)*(0.5-h1)
            h1z=0.5*(0.5-h3)*(0.5-h3)

c     FIELD INTERPOLATION

            exq=gmz*(hmx*ex(l1-1,j2,j3-1)
     &              +h0x*ex(l1,j2,j3-1)
     &              +h1x*ex(l1+1,j2,j3-1))
     &         +g0z*(hmx*ex(l1-1,j2,j3)
     &              +h0x*ex(l1,j2,j3)
     &              +h1x*ex(l1+1,j2,j3))
     &         +g1z*(hmx*ex(l1-1,j2,j3+1)
     &              +h0x*ex(l1,j2,j3+1)
     &              +h1x*ex(l1+1,j2,j3+1))

            eyq=gmz*(gmx*ey(j1-1,l2,j3-1)
     &              +g0x*ey(j1,l2,j3-1)
     &              +g1x*ey(j1+1,l2,j3-1))
     &         +g0z*(gmx*ey(j1-1,l2,j3)
     &              +g0x*ey(j1,l2,j3)
     &              +g1x*ey(j1+1,l2,j3))
     &         +g1z*(gmx*ey(j1-1,l2,j3+1)
     &              +g0x*ey(j1,l2,j3+1)
     &              +g1x*ey(j1+1,l2,j3+1))

            ezq=hmz*(gmx*ez(j1-1,j2,l3-1)
     &              +g0x*ez(j1,j2,l3-1)
     &              +g1x*ez(j1+1,j2,l3-1))
     &         +h0z*(gmx*ez(j1-1,j2,l3)
     &              +g0x*ez(j1,j2,l3)
     &              +g1x*ez(j1+1,j2,l3))
     &         +h1z*(gmx*ez(j1-1,j2,l3+1)
     &              +g0x*ez(j1,j2,l3+1)
     &              +g1x*ez(j1+1,j2,l3+1))

            bxq=hmz*(gmx*bx(j1-1,l2,l3-1)
     &              +g0x*bx(j1,l2,l3-1)
     &              +g1x*bx(j1+1,l2,l3-1))
     &         +h0z*(gmx*bx(j1-1,l2,l3)
     &              +g0x*bx(j1,l2,l3)
     &              +g1x*bx(j1+1,l2,l3))
     &         +h1z*(gmx*bx(j1-1,l2,l3+1)
     &              +g0x*bx(j1,l2,l3+1)
     &              +g1x*bx(j1+1,l2,l3+1))

            byq=hmz*(hmx*by(l1-1,j2,l3-1)
     &              +h0x*by(l1,j2,l3-1)
     &              +h1x*by(l1+1,j2,l3-1))
     &         +h0z*(hmx*by(l1-1,j2,l3)
     &              +h0x*by(l1,j2,l3)
     &              +h1x*by(l1+1,j2,l3))
     &         +h1z*(hmx*by(l1-1,j2,l3+1)
     &              +h0x*by(l1,j2,l3+1)
     &              +h1x*by(l1+1,j2,l3+1))

            bzq=gmz*(hmx*bz(l1-1,l2,j3-1)
     &              +h0x*bz(l1,l2,j3-1)
     &              +h1x*bz(l1+1,l2,j3-1))
     &         +g0z*(hmx*bz(l1-1,l2,j3)
     &              +h0x*bz(l1,l2,j3)
     &              +h1x*bz(l1+1,l2,j3))
     &         +g1z*(hmx*bz(l1-1,l2,j3+1)
     &              +h0x*bz(l1,l2,j3+1)
     &              +h1x*bz(l1+1,l2,j3+1))

c x^(n+0.5), p^n -> x^(n+1.0), p^(n+1.0) 

            dq=qni*dqs/mni
            pxm=pxi+dq*exq
            pym=pyi+dq*eyq
            pzm=pzi+dq*ezq

            root=dq/dsqrt(1.0+pxm*pxm+pym*pym+pzm*pzm)
            taux=bxq*root
            tauy=byq*root
            tauz=bzq*root

            tau=1.0/(1.0+taux*taux+tauy*tauy+tauz*tauz)
            pxp=((1.0+taux*taux-tauy*tauy-tauz*tauz)*pxm
     &          +(2.0*taux*tauy+2.0*tauz)*pym
     &          +(2.0*taux*tauz-2.0*tauy)*pzm)*tau
            pyp=((2.0*taux*tauy-2.0*tauz)*pxm
     &          +(1.0-taux*taux+tauy*tauy-tauz*tauz)*pym
     &          +(2.0*tauy*tauz+2.0*taux)*pzm)*tau
            pzp=((2.0*taux*tauz+2.0*tauy)*pxm
     &          +(2.0*tauy*tauz-2.0*taux)*pym
     &          +(1.0-taux*taux-tauy*tauy+tauz*tauz)*pzm)*tau

            pxi=pxp+dq*exq
            pyi=pyp+dq*eyq
            pzi=pzp+dq*ezq

            root=1.0/dsqrt(1.0+pxi*pxi+pyi*pyi+pzi*pzi)
            vxi=pxi*root
            vyi=pyi*root
            vzi=pzi*root

            xi=xi+vxi*xl
            zi=zi+vzi*zl

            p_niloc(11*l)=xi
            p_niloc(11*l+1)=yi
            p_niloc(11*l+2)=zi
            p_niloc(11*l+3)=pxi
            p_niloc(11*l+4)=pyi
            p_niloc(11*l+5)=pzi

            p2B=p2B+mni*fnqs*(1.0d0/root-1.0d0)/eta

c CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
c x^(n+1), p^(n+1) -> x^(n+1.5), p^(n+1)

            xi=xi+vxi*xl
            zi=zi+vzi*zl

            u=xi*dxi
            w=zi*dzi
            k1=nint(u)
            k3=nint(w)
            h1=k1-u
            h3=k3-w

            s1x(k1-j1-1)=0.5*(1.5-abs(h1-1.0))*(1.5-abs(h1-1.0))
            s1x(k1-j1+0)=0.75-abs(h1)*abs(h1)
            s1x(k1-j1+1)=0.5*(1.5-abs(h1+1.0))*(1.5-abs(h1+1.0))
            s1z(k3-j3-1)=0.5*(1.5-abs(h3-1.0))*(1.5-abs(h3-1.0))
            s1z(k3-j3+0)=0.75-abs(h3)*abs(h3)
            s1z(k3-j3+1)=0.5*(1.5-abs(h3+1.0))*(1.5-abs(h3+1.0))

c CURRENT DENSITY AT (n+1.0)*dt

            s1x=s1x-s0x
            s1z=s1z-s0z

            if (k1==j1) then
               l1min=-1
               l1max=+1
            else if (k1==j1-1) then
               l1min=-2
               l1max=+1
            else if (k1==j1+1) then
               l1min=-1
               l1max=+2
            endif
            if (k3==j3) then
               l3min=-1
               l3max=+1
            else if (k3==j3-1) then
               l3min=-2
               l3max=+1
            else if (k3==j3+1) then
               l3min=-1
               l3max=+2
            endif

            jxh=0.0
            jyh=0.0
            jzh=0.0

            fnqx=qni*wni*fnqxs
            fnqy=vyi*qni*wni*fnqs
            fnqz=qni*wni*fnqzs
            do l3=l3min,l3max
               do l1=l1min,l1max
                  wx=s1x(l1)*(s0z(l3)+0.5*s1z(l3))
                  wy=s0x(l1)*s0z(l3)
     &               +0.5*s1x(l1)*s0z(l3)
     &               +0.5*s0x(l1)*s1z(l3)
     &               +0.3333333333*s1x(l1)*s1z(l3)
                  wz=s1z(l3)*(s0x(l1)+0.5*s1x(l1))

                  jxh(l1,0,l3)=jxh(l1-1,0,l3)-fnqx*wx
                  jyh(l1,0,l3)=fnqy*wy
                  jzh(l1,0,l3)=jzh(l1,0,l3-1)-fnqz*wz

                  jxi(j1+l1,j2,j3+l3)=jxi(j1+l1,j2,j3+l3)
     &                                   +jxh(l1,0,l3)
                  jyi(j1+l1,j2,j3+l3)=jyi(j1+l1,j2,j3+l3)
     &                                   +jyh(l1,0,l3)
                  jzi(j1+l1,j2,j3+l3)=jzi(j1+l1,j2,j3+l3)
     &                                   +jzh(l1,0,l3)
               enddo
            enddo

         enddo
      endif



      deallocate(s0x,s0z)
      deallocate(s1x,s1z)
      deallocate(jxh,jyh,jzh)


      end subroutine PIC_push_part_xz

