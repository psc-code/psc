c CALCULATE IONIZATION X-SECTION FROM LOTZ FORMULA TO BE SIMPLE AND GENERAL
c interpolate x-section and sv == sigma v
c for material mat and initial charge state xs
c expect ek in eV
c return sigma and sv in SI units


      subroutine MCC_ixsection(mat,cs,ek,sigma,sv)

      use PIC_variables
      use VLA_variables
      use MCC_variables

      implicit none

      integer :: mat,cs                           ! expect real material and charge state index -- NOT mcc_list - index
      integer :: ux,lx,mx

      real(kind=8) :: sigma,sv,ek
      real(kind=8) :: va


      va=cc*sqrt(ek*ek+2.0*me*ek)/(ek+me)         ! va given in SI units
      ux=xstable_n(mat,cs)-1                      ! assume table is ordered for increasing energies
      lx=0


      if (ek < xstable_e(mat,cs,lx)) then         ! ekin < table boundary
         sigma=0.0d0
         sv=0.0d0
         return
      endif 

      if (ek > xstable_e(mat,cs,ux)) then         ! ekin > table boundary
         sigma=xstable_s(mat,cs,ux)
         sv=sigma*va
         return
      else 
         do while (ux > lx+1)                     ! perform linear bisection to find sigma
            mx=int((lx+ux)/2.0)
            if (xstable_e(mat,cs,mx) > ek) then
               ux=mx
            else
               lx=mx
            endif
         enddo
         sigma=((ek-xstable_e(mat,cs,lx))*xstable_s(mat,cs,ux)+
     &        (xstable_e(mat,cs,ux)-ek)*xstable_s(mat,cs,lx))/
     &        (xstable_e(mat,cs,ux)-xstable_e(mat,cs,lx))
         sv=sigma*va
         return
      endif

      end subroutine MCC_ixsection
