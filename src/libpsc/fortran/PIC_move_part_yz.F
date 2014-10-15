
c     ----------------------------------------------------------------------
      subroutine PIC_push_part_yz(niloc, p_niloc, p2A, p2B,
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

      integer :: j1,j2,j3,k2,k3,l1,l2,l3
      integer :: l2min,l2max,l3min,l3max
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
      real(kind=8) :: hmy,h0y,h1y,hmz,h0z,h1z
      real(kind=8) :: gmy,g0y,g1y,gmz,g0z,g1z
      real(kind=8) :: wx,wy,wz
      real(kind=8) :: exq,eyq,ezq
      real(kind=8) :: bxq,byq,bzq
      real(kind=8) :: taux,tauy,tauz,tau
      real(kind=8) :: u,v,w

      real(kind=8),dimension(-2:2) :: s0y,s0z
      real(kind=8),dimension(-2:2) :: s1y,s1z
      real(kind=8),dimension(-3:2,-2:2,-2:2) :: jxh
      real(kind=8),dimension(-2:2,-3:2,-2:2) :: jyh
      real(kind=8),dimension(-2:2,-2:2,-3:2) :: jzh


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

            yi=yi+vyi*yl
            zi=zi+vzi*zl

            s0y=0.0
            s0z=0.0
            s1y=0.0
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
            gmy=0.5*(0.5+h2)*(0.5+h2)
            gmz=0.5*(0.5+h3)*(0.5+h3)
            g0y=0.75-h2*h2
            g0z=0.75-h3*h3
            g1y=0.5*(0.5-h2)*(0.5-h2)
            g1z=0.5*(0.5-h3)*(0.5-h3)

            s0y(-1)=0.5*(1.5-abs(h2-1.0))*(1.5-abs(h2-1.0))
            s0y(+0)=0.75-abs(h2)*abs(h2)
            s0y(+1)=0.5*(1.5-abs(h2+1.0))*(1.5-abs(h2+1.0))
            s0z(-1)=0.5*(1.5-abs(h3-1.0))*(1.5-abs(h3-1.0))
            s0z(+0)=0.75-abs(h3)*abs(h3)
            s0z(+1)=0.5*(1.5-abs(h3+1.0))*(1.5-abs(h3+1.0))

            u=xi*dxi
            v=yi*dyi-0.5
            w=zi*dzi-0.5
            l1=nint(u)
            l2=nint(v)
            l3=nint(w)
            h1=l1-u
            h2=l2-v
            h3=l3-w
            hmy=0.5*(0.5+h2)*(0.5+h2)
            hmz=0.5*(0.5+h3)*(0.5+h3)
            h0y=0.75-h2*h2
            h0z=0.75-h3*h3
            h1y=0.5*(0.5-h2)*(0.5-h2)
            h1z=0.5*(0.5-h3)*(0.5-h3)

c     FIELD INTERPOLATION

            exq=gmz*(gmy*ex(l1,j2-1,j3-1)
     &              +g0y*ex(l1,j2,j3-1)
     &              +g1y*ex(l1,j2+1,j3-1))
     &         +g0z*(gmy*ex(l1,j2-1,j3)
     &              +g0y*ex(l1,j2,j3)
     &              +g1y*ex(l1,j2+1,j3))
     &         +g1z*(gmy*ex(l1,j2-1,j3+1)
     &              +g0y*ex(l1,j2,j3+1)
     &              +g1y*ex(l1,j2+1,j3+1))

            eyq=gmz*(hmy*ey(j1,l2-1,j3-1)
     &              +h0y*ey(j1,l2,j3-1)
     &              +h1y*ey(j1,l2+1,j3-1))
     &         +g0z*(hmy*ey(j1,l2-1,j3)
     &              +h0y*ey(j1,l2,j3)
     &              +h1y*ey(j1,l2+1,j3))
     &         +g1z*(hmy*ey(j1,l2-1,j3+1)
     &              +h0y*ey(j1,l2,j3+1)
     &              +h1y*ey(j1,l2+1,j3+1))

            ezq=hmz*(gmy*ez(j1,j2-1,l3-1)
     &              +g0y*ez(j1,j2,l3-1)
     &              +g1y*ez(j1,j2+1,l3-1))
     &         +h0z*(gmy*ez(j1,j2-1,l3)
     &              +g0y*ez(j1,j2,l3)
     &              +g1y*ez(j1,j2+1,l3))
     &         +h1z*(gmy*ez(j1,j2-1,l3+1)
     &              +g0y*ez(j1,j2,l3+1)
     &              +g1y*ez(j1,j2+1,l3+1))

            bxq=hmz*(hmy*bx(j1,l2-1,l3-1)
     &              +h0y*bx(j1,l2,l3-1)
     &              +h1y*bx(j1,l2+1,l3-1))
     &         +h0z*(hmy*bx(j1,l2-1,l3)
     &              +h0y*bx(j1,l2,l3)
     &              +h1y*bx(j1,l2+1,l3))
     &         +h1z*(hmy*bx(j1,l2-1,l3+1)
     &              +h0y*bx(j1,l2,l3+1)
     &              +h1y*bx(j1,l2+1,l3+1))

            byq=hmz*(gmy*by(l1,j2-1,l3-1)
     &              +g0y*by(l1,j2,l3-1)
     &              +g1y*by(l1,j2+1,l3-1))
     &         +h0z*(gmy*by(l1,j2-1,l3)
     &              +g0y*by(l1,j2,l3)
     &              +g1y*by(l1,j2+1,l3))
     &         +h1z*(gmy*by(l1,j2-1,l3+1)
     &              +g0y*by(l1,j2,l3+1)
     &              +g1y*by(l1,j2+1,l3+1))

            bzq=gmz*(hmy*bz(l1,l2-1,j3-1)
     &              +h0y*bz(l1,l2,j3-1)
     &              +h1y*bz(l1,l2+1,j3-1))
     &         +g0z*(hmy*bz(l1,l2-1,j3)
     &              +h0y*bz(l1,l2,j3)
     &              +h1y*bz(l1,l2+1,j3))
     &         +g1z*(hmy*bz(l1,l2-1,j3+1)
     &              +h0y*bz(l1,l2,j3+1)
     &              +h1y*bz(l1,l2+1,j3+1))

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

            yi=yi+vyi*yl
            zi=zi+vzi*zl
            
            p_niloc(11*l)=xi
            p_niloc(11*l+1)=yi
            p_niloc(11*l+2)=zi
            p_niloc(11*l+3)=pxi
            p_niloc(11*l+4)=pyi
            p_niloc(11*l+5)=pzi


            p2B=p2B+mni*fnqs*(1.0d0/root-1.0d0)/eta

c DETERMINE THE DENSITIES AT t=(n+1.0)*dt

            u=xi*dxi
            v=yi*dyi
            w=zi*dzi
            l1=nint(u)
            l2=nint(v)
            l3=nint(w)
            h2=l2-v
            h3=l3-w

            gmy=0.5*(1.5-abs(h2-1.0))*(1.5-abs(h2-1.0))
            g0y=0.75-abs(h2)*abs(h2)
            g1y=0.5*(1.5-abs(h2+1.0))*(1.5-abs(h2+1.0))
            gmz=0.5*(1.5-abs(h3-1.0))*(1.5-abs(h3-1.0))
            g0z=0.75-abs(h3)*abs(h3)
            g1z=0.5*(1.5-abs(h3+1.0))*(1.5-abs(h3+1.0))

c CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
c x^(n+1), p^(n+1) -> x^(n+1.5), p^(n+1)

            yi=yi+vyi*yl
            zi=zi+vzi*zl

            v=yi*dyi
            w=zi*dzi
            k2=nint(v)
            k3=nint(w)
            h2=k2-v
            h3=k3-w

            s1y(k2-j2-1)=0.5*(1.5-abs(h2-1.0))*(1.5-abs(h2-1.0))
            s1y(k2-j2+0)=0.75-abs(h2)*abs(h2)
            s1y(k2-j2+1)=0.5*(1.5-abs(h2+1.0))*(1.5-abs(h2+1.0))
            s1z(k3-j3-1)=0.5*(1.5-abs(h3-1.0))*(1.5-abs(h3-1.0))
            s1z(k3-j3+0)=0.75-abs(h3)*abs(h3)
            s1z(k3-j3+1)=0.5*(1.5-abs(h3+1.0))*(1.5-abs(h3+1.0))

c CURRENT DENSITY AT (n+1.0)*dt

            s1y=s1y-s0y
            s1z=s1z-s0z

            if (k2==j2) then
               l2min=-1
               l2max=+1
            else if (k2==j2-1) then
               l2min=-2
               l2max=+1
            else if (k2==j2+1) then
               l2min=-1
               l2max=+2
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

            fnqx=vxi*qni*wni*fnqs
            fnqy=qni*wni*fnqys
            fnqz=qni*wni*fnqzs
            do l3=l3min,l3max
               do l2=l2min,l2max
                  wx=s0y(l2)*s0z(l3)
     &               +0.5*s1y(l2)*s0z(l3)
     &               +0.5*s0y(l2)*s1z(l3)
     &               +0.3333333333*s1y(l2)*s1z(l3)
                  wy=s1y(l2)*(s0z(l3)+0.5*s1z(l3))
                  wz=s1z(l3)*(s0y(l2)+0.5*s1y(l2))

                  jxh(0,l2,l3)=fnqx*wx
                  jyh(0,l2,l3)=jyh(0,l2-1,l3)-fnqy*wy
                  jzh(0,l2,l3)=jzh(0,l2,l3-1)-fnqz*wz

                  jxi(j1,j2+l2,j3+l3)=jxi(j1,j2+l2,j3+l3)
     &                                   +jxh(0,l2,l3)
                  jyi(j1,j2+l2,j3+l3)=jyi(j1,j2+l2,j3+l3)
     &                                   +jyh(0,l2,l3)
                  jzi(j1,j2+l2,j3+l3)=jzi(j1,j2+l2,j3+l3)
     &                                   +jzh(0,l2,l3)
               enddo
            enddo

         enddo
      endif

      end subroutine



