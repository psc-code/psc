
c     ----------------------------------------------------------------------
      subroutine PIC_push_part_xyz(niloc, p_niloc, p2A, p2B,
     &     jxi, jyi, jzi, ex, ey, ez, hx, hy, hz)
c     ----------------------------------------------------------------------
c     particle pusher xyz, w/o the following boundary exchange

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
     &                       i3mn-rd3:i3mx+rd3) :: hx, hy, hz

      integer :: j1,j2,j3,k1,k2,k3,l1,l2,l3
      integer :: l1min,l1max,l2min,l2max,l3min,l3max
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
      real(kind=8) :: hmx,h0x,h1x,hmy,h0y,h1y,hmz,h0z,h1z
      real(kind=8) :: gmx,g0x,g1x,gmy,g0y,g1y,gmz,g0z,g1z
      real(kind=8) :: wx,wy,wz
      real(kind=8) :: exq,eyq,ezq
!      real(kind=8) :: bxq,byq,bzq
      real(kind=8) :: hxq,hyq,hzq
      real(kind=8) :: taux,tauy,tauz,tau
      real(kind=8) :: u,v,w

      real(kind=8) :: s_cpub
      real(kind=8) :: s_cpuh

      real(kind=8),allocatable,dimension(:) :: s0x,s0y,s0z


      real(kind=8),allocatable,dimension(:) :: s1x,s1y,s1z
      real(kind=8),allocatable,dimension(:,:,:) :: jxh,jyh,jzh


      allocate(s0x(-2:2))
      allocate(s0y(-2:2))
      allocate(s0z(-2:2))
      allocate(s1x(-2:2))
      allocate(s1y(-2:2))
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
            yi=yi+vyi*yl
            zi=zi+vzi*zl

            s0x=0.0
            s0y=0.0
            s0z=0.0
            s1x=0.0
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
            gmx=0.5*(0.5+h1)*(0.5+h1)
            gmy=0.5*(0.5+h2)*(0.5+h2)
            gmz=0.5*(0.5+h3)*(0.5+h3)
            g0x=0.75-h1*h1
            g0y=0.75-h2*h2
            g0z=0.75-h3*h3
            g1x=0.5*(0.5-h1)*(0.5-h1)
            g1y=0.5*(0.5-h2)*(0.5-h2)
            g1z=0.5*(0.5-h3)*(0.5-h3)

            s0x(-1)=0.5*(1.5-abs(h1-1.0))*(1.5-abs(h1-1.0))
            s0x(+0)=0.75-abs(h1)*abs(h1)
            s0x(+1)=0.5*(1.5-abs(h1+1.0))*(1.5-abs(h1+1.0))
            s0y(-1)=0.5*(1.5-abs(h2-1.0))*(1.5-abs(h2-1.0))
            s0y(+0)=0.75-abs(h2)*abs(h2)
            s0y(+1)=0.5*(1.5-abs(h2+1.0))*(1.5-abs(h2+1.0))
            s0z(-1)=0.5*(1.5-abs(h3-1.0))*(1.5-abs(h3-1.0))
            s0z(+0)=0.75-abs(h3)*abs(h3)
            s0z(+1)=0.5*(1.5-abs(h3+1.0))*(1.5-abs(h3+1.0))

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

c     FIELD INTERPOLATION

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

!            bxq=hmz*(hmy*(gmx*bx(j1-1,l2-1,l3-1)
!     &                   +g0x*bx(j1,l2-1,l3-1)
!     &                   +g1x*bx(j1+1,l2-1,l3-1))
!     &              +h0y*(gmx*bx(j1-1,l2,l3-1)
!     &                   +g0x*bx(j1,l2,l3-1)
!     &                   +g1x*bx(j1+1,l2,l3-1))
!     &              +h1y*(gmx*bx(j1-1,l2+1,l3-1)
!     &                   +g0x*bx(j1,l2+1,l3-1)
!     &                   +g1x*bx(j1+1,l2+1,l3-1)))
!     &         +h0z*(hmy*(gmx*bx(j1-1,l2-1,l3)
!     &                   +g0x*bx(j1,l2-1,l3)
!     &                   +g1x*bx(j1+1,l2-1,l3))
!     &              +h0y*(gmx*bx(j1-1,l2,l3)
!     &                   +g0x*bx(j1,l2,l3)
!     &                   +g1x*bx(j1+1,l2,l3))
!     &              +h1y*(gmx*bx(j1-1,l2+1,l3)
!     &                   +g0x*bx(j1,l2+1,l3)
!     &                   +g1x*bx(j1+1,l2+1,l3)))
!     &         +h1z*(hmy*(gmx*bx(j1-1,l2-1,l3+1)
!     &                   +g0x*bx(j1,l2-1,l3+1)
!     &                   +g1x*bx(j1+1,l2-1,l3+1))
!     &              +h0y*(gmx*bx(j1-1,l2,l3+1)
!     &                   +g0x*bx(j1,l2,l3+1)
!     &                   +g1x*bx(j1+1,l2,l3+1))
!     &              +h1y*(gmx*bx(j1-1,l2+1,l3+1)
!     &                   +g0x*bx(j1,l2+1,l3+1)
!     &                   +g1x*bx(j1+1,l2+1,l3+1)))

!            byq=hmz*(gmy*(hmx*by(l1-1,j2-1,l3-1)
!     &                   +h0x*by(l1,j2-1,l3-1)
!     &                   +h1x*by(l1+1,j2-1,l3-1))
!     &              +g0y*(hmx*by(l1-1,j2,l3-1)
!     &                   +h0x*by(l1,j2,l3-1)
!     &                   +h1x*by(l1+1,j2,l3-1))
!     &              +g1y*(hmx*by(l1-1,j2+1,l3-1)
!     &                   +h0x*by(l1,j2+1,l3-1)
!     &                   +h1x*by(l1+1,j2+1,l3-1)))
!     &         +h0z*(gmy*(hmx*by(l1-1,j2-1,l3)
!     &                   +h0x*by(l1,j2-1,l3)
!     &                   +h1x*by(l1+1,j2-1,l3))
!     &              +g0y*(hmx*by(l1-1,j2,l3)
!     &                   +h0x*by(l1,j2,l3)
!     &                   +h1x*by(l1+1,j2,l3))
!     &              +g1y*(hmx*by(l1-1,j2+1,l3)
!     &                   +h0x*by(l1,j2+1,l3)
!     &                   +h1x*by(l1+1,j2+1,l3)))
!     &         +h1z*(gmy*(hmx*by(l1-1,j2-1,l3+1)
!     &                   +h0x*by(l1,j2-1,l3+1)
!     &                   +h1x*by(l1+1,j2-1,l3+1))
!     &              +g0y*(hmx*by(l1-1,j2,l3+1)
!     &                   +h0x*by(l1,j2,l3+1)
!     &                   +h1x*by(l1+1,j2,l3+1))
!     &              +g1y*(hmx*by(l1-1,j2+1,l3+1)
!     &                   +h0x*by(l1,j2+1,l3+1)
!     &                   +h1x*by(l1+1,j2+1,l3+1)))

!            bzq=gmz*(hmy*(hmx*bz(l1-1,l2-1,j3-1)
!     &                   +h0x*bz(l1,l2-1,j3-1)
!     &                   +h1x*bz(l1+1,l2-1,j3-1))
!     &              +h0y*(hmx*bz(l1-1,l2,j3-1)
!     &                   +h0x*bz(l1,l2,j3-1)
!     &                   +h1x*bz(l1+1,l2,j3-1))
!     &              +h1y*(hmx*bz(l1-1,l2+1,j3-1)
!     &                   +h0x*bz(l1,l2+1,j3-1)
!     &                   +h1x*bz(l1+1,l2+1,j3-1)))
!     &         +g0z*(hmy*(hmx*bz(l1-1,l2-1,j3)
!     &                   +h0x*bz(l1,l2-1,j3)
!     &                   +h1x*bz(l1+1,l2-1,j3))
!     &              +h0y*(hmx*bz(l1-1,l2,j3)
!     &                   +h0x*bz(l1,l2,j3)
!     &                   +h1x*bz(l1+1,l2,j3))
!     &              +h1y*(hmx*bz(l1-1,l2+1,j3)
!     &                   +h0x*bz(l1,l2+1,j3)
!     &                   +h1x*bz(l1+1,l2+1,j3)))
!     &         +g1z*(hmy*(hmx*bz(l1-1,l2-1,j3+1)
!     &                   +h0x*bz(l1,l2-1,j3+1)
!     &                   +h1x*bz(l1+1,l2-1,j3+1))
!     &              +h0y*(hmx*bz(l1-1,l2,j3+1)
!     &                   +h0x*bz(l1,l2,j3+1)
!     &                   +h1x*bz(l1+1,l2,j3+1))
!     &              +h1y*(hmx*bz(l1-1,l2+1,j3+1)
!     &                   +h0x*bz(l1,l2+1,j3+1)
!     &                   +h1x*bz(l1+1,l2+1,j3+1)))

            hxq=hmz*(hmy*(gmx*hx(j1-1,l2-1,l3-1)
     &                   +g0x*hx(j1,l2-1,l3-1)
     &                   +g1x*hx(j1+1,l2-1,l3-1))
     &              +h0y*(gmx*hx(j1-1,l2,l3-1)
     &                   +g0x*hx(j1,l2,l3-1)
     &                   +g1x*hx(j1+1,l2,l3-1))
     &              +h1y*(gmx*hx(j1-1,l2+1,l3-1)
     &                   +g0x*hx(j1,l2+1,l3-1)
     &                   +g1x*hx(j1+1,l2+1,l3-1)))
     &         +h0z*(hmy*(gmx*hx(j1-1,l2-1,l3)
     &                   +g0x*hx(j1,l2-1,l3)
     &                   +g1x*hx(j1+1,l2-1,l3))
     &              +h0y*(gmx*hx(j1-1,l2,l3)
     &                   +g0x*hx(j1,l2,l3)
     &                   +g1x*hx(j1+1,l2,l3))
     &              +h1y*(gmx*hx(j1-1,l2+1,l3)
     &                   +g0x*hx(j1,l2+1,l3)
     &                   +g1x*hx(j1+1,l2+1,l3)))
     &         +h1z*(hmy*(gmx*hx(j1-1,l2-1,l3+1)
     &                   +g0x*hx(j1,l2-1,l3+1)
     &                   +g1x*hx(j1+1,l2-1,l3+1))
     &              +h0y*(gmx*hx(j1-1,l2,l3+1)
     &                   +g0x*hx(j1,l2,l3+1)
     &                   +g1x*hx(j1+1,l2,l3+1))
     &              +h1y*(gmx*hx(j1-1,l2+1,l3+1)
     &                   +g0x*hx(j1,l2+1,l3+1)
     &                   +g1x*hx(j1+1,l2+1,l3+1)))

            hyq=hmz*(gmy*(hmx*hy(l1-1,j2-1,l3-1)
     &                   +h0x*hy(l1,j2-1,l3-1)
     &                   +h1x*hy(l1+1,j2-1,l3-1))
     &              +g0y*(hmx*hy(l1-1,j2,l3-1)
     &                   +h0x*hy(l1,j2,l3-1)
     &                   +h1x*hy(l1+1,j2,l3-1))
     &              +g1y*(hmx*hy(l1-1,j2+1,l3-1)
     &                   +h0x*hy(l1,j2+1,l3-1)
     &                   +h1x*hy(l1+1,j2+1,l3-1)))
     &         +h0z*(gmy*(hmx*hy(l1-1,j2-1,l3)
     &                   +h0x*hy(l1,j2-1,l3)
     &                   +h1x*hy(l1+1,j2-1,l3))
     &              +g0y*(hmx*hy(l1-1,j2,l3)
     &                   +h0x*hy(l1,j2,l3)
     &                   +h1x*hy(l1+1,j2,l3))
     &              +g1y*(hmx*hy(l1-1,j2+1,l3)
     &                   +h0x*hy(l1,j2+1,l3)
     &                   +h1x*hy(l1+1,j2+1,l3)))
     &         +h1z*(gmy*(hmx*hy(l1-1,j2-1,l3+1)
     &                   +h0x*hy(l1,j2-1,l3+1)
     &                   +h1x*hy(l1+1,j2-1,l3+1))
     &              +g0y*(hmx*hy(l1-1,j2,l3+1)
     &                   +h0x*hy(l1,j2,l3+1)
     &                   +h1x*hy(l1+1,j2,l3+1))
     &              +g1y*(hmx*hy(l1-1,j2+1,l3+1)
     &                   +h0x*hy(l1,j2+1,l3+1)
     &                   +h1x*hy(l1+1,j2+1,l3+1)))

            hzq=gmz*(hmy*(hmx*hz(l1-1,l2-1,j3-1)
     &                   +h0x*hz(l1,l2-1,j3-1)
     &                   +h1x*hz(l1+1,l2-1,j3-1))
     &              +h0y*(hmx*hz(l1-1,l2,j3-1)
     &                   +h0x*hz(l1,l2,j3-1)
     &                   +h1x*hz(l1+1,l2,j3-1))
     &              +h1y*(hmx*hz(l1-1,l2+1,j3-1)
     &                   +h0x*hz(l1,l2+1,j3-1)
     &                   +h1x*hz(l1+1,l2+1,j3-1)))
     &         +g0z*(hmy*(hmx*hz(l1-1,l2-1,j3)
     &                   +h0x*hz(l1,l2-1,j3)
     &                   +h1x*hz(l1+1,l2-1,j3))
     &              +h0y*(hmx*hz(l1-1,l2,j3)
     &                   +h0x*hz(l1,l2,j3)
     &                   +h1x*hz(l1+1,l2,j3))
     &              +h1y*(hmx*hz(l1-1,l2+1,j3)
     &                   +h0x*hz(l1,l2+1,j3)
     &                   +h1x*hz(l1+1,l2+1,j3)))
     &         +g1z*(hmy*(hmx*hz(l1-1,l2-1,j3+1)
     &                   +h0x*hz(l1,l2-1,j3+1)
     &                   +h1x*hz(l1+1,l2-1,j3+1))
     &              +h0y*(hmx*hz(l1-1,l2,j3+1)
     &                   +h0x*hz(l1,l2,j3+1)
     &                   +h1x*hz(l1+1,l2,j3+1))
     &              +h1y*(hmx*hz(l1-1,l2+1,j3+1)
     &                   +h0x*hz(l1,l2+1,j3+1)
     &                   +h1x*hz(l1+1,l2+1,j3+1)))

c x^(n+0.5), p^n -> x^(n+1.0), p^(n+1.0) 

            dq=qni*dqs/mni
            pxm=pxi+dq*exq
            pym=pyi+dq*eyq
            pzm=pzi+dq*ezq

            root=dq/dsqrt(1.0+pxm*pxm+pym*pym+pzm*pzm)
!            taux=bxq*root
!            tauy=byq*root
!            tauz=bzq*root
            taux=hxq*root
            tauy=hyq*root
            tauz=hzq*root

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
            yi=yi+vyi*yl
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
            yi=yi+vyi*yl
            zi=zi+vzi*zl

            u=xi*dxi
            v=yi*dyi
            w=zi*dzi
            k1=nint(u)
            k2=nint(v)
            k3=nint(w)
            h1=k1-u
            h2=k2-v
            h3=k3-w

            s1x(k1-j1-1)=0.5*(1.5-abs(h1-1.0))*(1.5-abs(h1-1.0))
            s1x(k1-j1+0)=0.75-abs(h1)*abs(h1)
            s1x(k1-j1+1)=0.5*(1.5-abs(h1+1.0))*(1.5-abs(h1+1.0))
            s1y(k2-j2-1)=0.5*(1.5-abs(h2-1.0))*(1.5-abs(h2-1.0))
            s1y(k2-j2+0)=0.75-abs(h2)*abs(h2)
            s1y(k2-j2+1)=0.5*(1.5-abs(h2+1.0))*(1.5-abs(h2+1.0))
            s1z(k3-j3-1)=0.5*(1.5-abs(h3-1.0))*(1.5-abs(h3-1.0))
            s1z(k3-j3+0)=0.75-abs(h3)*abs(h3)
            s1z(k3-j3+1)=0.5*(1.5-abs(h3+1.0))*(1.5-abs(h3+1.0))

c CURRENT DENSITY AT (n+1.0)*dt

            s1x=s1x-s0x
            s1y=s1y-s0y
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

            fnqx=qni*wni*fnqxs
            fnqy=qni*wni*fnqys
            fnqz=qni*wni*fnqzs
            do l3=l3min,l3max
               do l2=l2min,l2max
                  do l1=l1min,l1max
                     wx=s1x(l1)*(s0y(l2)*s0z(l3)
     &                  +0.5*s1y(l2)*s0z(l3)
     &                  +0.5*s0y(l2)*s1z(l3)
     &                  +0.3333333333*s1y(l2)*s1z(l3))
                     wy=s1y(l2)*(s0x(l1)*s0z(l3)
     &                  +0.5*s1x(l1)*s0z(l3)
     &                  +0.5*s0x(l1)*s1z(l3)
     &                  +0.3333333333*s1x(l1)*s1z(l3))
                     wz=s1z(l3)*(s0x(l1)*s0y(l2)
     &                  +0.5*s1x(l1)*s0y(l2)
     &                  +0.5*s0x(l1)*s1y(l2)
     &                  +0.3333333333*s1x(l1)*s1y(l2))

                     jxh(l1,l2,l3)=jxh(l1-1,l2,l3)-fnqx*wx
                     jyh(l1,l2,l3)=jyh(l1,l2-1,l3)-fnqy*wy
                     jzh(l1,l2,l3)=jzh(l1,l2,l3-1)-fnqz*wz

                     jxi(j1+l1,j2+l2,j3+l3)=jxi(j1+l1,j2+l2,j3+l3)
     &                                      +jxh(l1,l2,l3)
                     jyi(j1+l1,j2+l2,j3+l3)=jyi(j1+l1,j2+l2,j3+l3)
     &                                      +jyh(l1,l2,l3)
                     jzi(j1+l1,j2+l2,j3+l3)=jzi(j1+l1,j2+l2,j3+l3)
     &                                      +jzh(l1,l2,l3)
                  enddo
               enddo
            enddo

         enddo
      endif

      deallocate(s0x,s0y,s0z)
      deallocate(s1x,s1y,s1z)
      deallocate(jxh,jyh,jzh)

      end subroutine PIC_push_part_xyz

