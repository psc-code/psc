
c THIS SUBROUTINE INITIALIZES THE ELECTROMAGNETIC FIELDS.


c     ----------------------------------------------------------------------
      subroutine SETUP_field
c     ----------------------------------------------------------------------
c     set up various things after the main fields have been initialized

      use PSC_domain

      implicit none


      if (use_pml) then
         call PML_coeff_alloc
         call PML_coeff_setup
      endif

      end subroutine SETUP_field

