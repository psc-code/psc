
c     ----------------------------------------------------------------------
      subroutine PIC_push_part_z_vay(niloc, p_niloc, p2A, p2B,
     &     jxi, jyi, jzi, ex, ey, ez, hx, hy, hz)
c     ----------------------------------------------------------------------
c     particle pusher z, w/o the following boundary exchange

      use VLA_variables, only: n
      use PSC_globals
      use PSC_params
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

      real(kind=8) :: pi

      integer :: j1,j2,j3,k2,k3,l1,l2,l3
      integer :: l2min,l2max,l3min,l3max
      integer :: l

      real(kind=8) :: dxi,dyi,dzi
      real(kind=8) :: pxi,pyi,pzi
      real(kind=8) :: qni,mni,cni,lni,wni
      real(kind=8) :: xi,yi,zi,vxi,vyi,vzi,root
      real(kind=8) :: xl,yl,zl
			real(kind=8) :: usx, usy, usz
			real(kind=8) :: gammas2, sigma, sigma2
			real(kind=8) :: ustau, ustau2
			real(kind=8) :: uscrosstaux, uscrosstauy, uscrosstauz
      
      real(kind=8) :: dqs,fnqs,fnqxs,fnqys,fnqzs
      real(kind=8) :: dq,fnq,fnqx,fnqy,fnqz

      real(kind=8) :: h1,h2,h3
      real(kind=8) :: hmy,h0y,h1y,hmz,h0z,h1z
      real(kind=8) :: gmy,g0y,g1y,gmz,g0z,g1z
      real(kind=8) :: wx,wy,wz
      real(kind=8) :: exq,eyq,ezq
      real(kind=8) :: hxq,hyq,hzq
      real(kind=8) :: taux,tauy,tauz,tau
      real(kind=8) :: u,v,w

      real(kind=8),dimension(-2:2) :: s0y,s0z
      real(kind=8),dimension(-2:2) :: s1y,s1z
      real(kind=8),dimension(-3:2,-2:2,-2:2) :: jxh
      real(kind=8),dimension(-2:2,-3:2,-2:2) :: jyh
      real(kind=8),dimension(-2:2,-2:2,-3:2) :: jzh

c INITIALIZATION

      pi=3.1415927

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

            zi=zi+vzi*zl

            s0z=0.0
            s1z=0.0

            u=xi*dxi
            v=yi*dyi
            w=zi*dzi
            j1=nint(u)
            j2=nint(v)
            j3=nint(w)
            h3=j3-w
            gmz=0.5*(0.5+h3)*(0.5+h3)
            g0z=0.75-h3*h3
            g1z=0.5*(0.5-h3)*(0.5-h3)

            s0z(-1)=0.5*(1.5-abs(h3-1.0))*(1.5-abs(h3-1.0))
            s0z(+0)=0.75-abs(h3)*abs(h3)
            s0z(+1)=0.5*(1.5-abs(h3+1.0))*(1.5-abs(h3+1.0))

            u=xi*dxi
            v=yi*dyi
            w=zi*dzi-0.5
            l1=nint(u)
            l2=nint(v)
            l3=nint(w)
            h3=l3-w
            hmz=0.5*(0.5+h3)*(0.5+h3)
            h0z=0.75-h3*h3
            h1z=0.5*(0.5-h3)*(0.5-h3)

c     FIELD INTERPOLATION

            exq=gmz*ex(l1,j2,j3-1)
     &          +g0z*ex(l1,j2,j3)
     &          +g1z*ex(l1,j2,j3+1)

            eyq=gmz*ey(j1,l2,j3-1)
     &          +g0z*ey(j1,l2,j3)
     &          +g1z*ey(j1,l2,j3+1)

            ezq=hmz*ez(j1,j2,l3-1)
     &         +h0z*ez(j1,j2,l3)
     &         +h1z*ez(j1,j2,l3+1)

!            bxq=hmz*bx(j1,l2,l3-1)
!     &          +h0z*bx(j1,l2,l3)
!     &          +h1z*bx(j1,l2,l3+1)

!            byq=hmz*by(l1,j2,l3-1)
!     &         +h0z*by(l1,j2,l3)
!     &         +h1z*by(l1,j2,l3+1)

!            bzq=gmz*bz(l1,l2,j3-1)
!     &          +g0z*bz(l1,l2,j3)
!     &          +g1z*bz(l1,l2,j3+1)

            hxq=hmz*hx(j1,l2,l3-1)
     &          +h0z*hx(j1,l2,l3)
     &          +h1z*hx(j1,l2,l3+1)

            hyq=hmz*hy(l1,j2,l3-1)
     &         +h0z*hy(l1,j2,l3)
     &         +h1z*hy(l1,j2,l3+1)

            hzq=gmz*hz(l1,l2,j3-1)
     &          +g0z*hz(l1,l2,j3)
     &          +g1z*hz(l1,l2,j3+1)

c x^(n+0.5), p^n -> x^(n+1.0), p^(n+1.0) 

c            if (n*dt < 2.0*pi) then
c               ezq=ezq+0.01*n*dt*dsin(1.0*wp*n*dt/wl)/(2.0*pi)
c            else
c               ezq=ezq+0.01*dsin(1.0*wp*n*dt/wl)
c            endif


c this is the method from Vay, Physics of Plasmas 2008


						dq=qni/mni*dt*eta                 ! eta is essentially a0, so it is here, not in fields

						usx=pxi+dq*exq+0.5*dq*(vyi*hzq-vzi*hyq)  ! u' x component
						usy=pyi+dq*eyq+0.5*dq*(vzi*hxq-vxi*hzq)  ! u' y component
						usz=pzi+dq*ezq+0.5*dq*(vxi*hyq-vyi*hxq)  ! u' z component
						
						taux=0.5*dq*hxq   ! tau x component
						tauy=0.5*dq*hyq   ! tau y component
						tauz=0.5*dq*hzq   ! tau z component
						
						tau=taux*taux+tauy*tauy+tauz*tauz    ! tau squared
						
						gammas2=(1.0+usx*usx+usy*usy+usz*usz)     ! gamma' squared
						sigma=gammas2-tau        
						sigma2=sigma*sigma
						
						ustau=usx*taux+usy*tauy+usz*tauz           ! u' dot tau
						ustau2=ustau*ustau
						
						root=1/dsqrt(0.5*(sigma+dsqrt(sigma2+4.0*(tau+ustau2))))   ! see Vay for details, reciprocal gamma
						
						taux=taux*root															
						tauy=tauy*root
						tauz=tauz*root
						sigma=1./(1+taux*taux+tauy*tauy+tauz*tauz)
						
						ustau=usx*taux+usy*tauy+usz*tauz
						
						uscrosstaux=usy*tauz-usz*tauy
						uscrosstauy=usz*taux-usx*tauz
						uscrosstauz=usx*tauy-usy*taux
						
						pxi=sigma*(usx+taux*ustau+uscrosstaux)
						pyi=sigma*(usy+tauy*ustau+uscrosstauy)
						pzi=sigma*(usz+tauz*ustau+uscrosstauz)
						
						root=1.0/dsqrt(1.0+pxi*pxi+pyi*pyi+pzi*pzi)
            vxi=pxi*root
            vyi=pyi*root
            vzi=pzi*root

            zi=zi+vzi*zl
						
c end Vay method

            

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
            h3=l3-w

            gmz=0.5*(1.5-abs(h3-1.0))*(1.5-abs(h3-1.0))
            g0z=0.75-abs(h3)*abs(h3)
            g1z=0.5*(1.5-abs(h3+1.0))*(1.5-abs(h3+1.0))

c CHARGE DENSITY FORM FACTOR AT (n+1.5)*dt 
c x^(n+1), p^(n+1) -> x^(n+1.5), p^(n+1)

            zi=zi+vzi*zl

            w=zi*dzi
            k3=nint(w)
            h3=k3-w

            s1z(k3-j3-1)=0.5*(1.5-abs(h3-1.0))*(1.5-abs(h3-1.0))
            s1z(k3-j3+0)=0.75-abs(h3)*abs(h3)
            s1z(k3-j3+1)=0.5*(1.5-abs(h3+1.0))*(1.5-abs(h3+1.0))

c CURRENT DENSITY AT (n+1.0)*dt

            s1z=s1z-s0z

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
            fnqy=vyi*qni*wni*fnqs
            fnqz=qni*wni*fnqzs
            do l3=l3min,l3max
               wx=s0z(l3)+0.5*s1z(l3)
               wy=s0z(l3)+0.5*s1z(l3)
               wz=s1z(l3)

               jxh(0,0,l3)=fnqx*wx
               jyh(0,0,l3)=fnqy*wy
               jzh(0,0,l3)=jzh(0,0,l3-1)-fnqz*wz

               jxi(j1,j2,j3+l3)=jxi(j1,j2,j3+l3)
     &                                +jxh(0,0,l3)
               jyi(j1,j2,j3+l3)=jyi(j1,j2,j3+l3)
     &                                +jyh(0,0,l3)
               jzi(j1,j2,j3+l3)=jzi(j1,j2,j3+l3)
     &                                +jzh(0,0,l3)
            enddo

         enddo
      endif

      end subroutine

