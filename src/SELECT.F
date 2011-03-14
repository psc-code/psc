c SELECT is used for producing the graphics data files
c THE PROGRAM SELECTS DISTRIBUTED DATA FILES PRODUCED BY THE
c OUT-ROUTINES AND MERGES THEM INTO GLOBAL DATA FILES!


      program SELECT

      implicit none

      integer :: ierr
      character*(20) char

      call MPI_INIT(ierr)
      call INIT_basic
      call INIT_param
      call INIT_param_fortran
      

c SELECT DATA FILES


 110  write(6,*) ' '
      write(6,*) 'time evolution of fields -> a1'
      write(6,*) 'time evolution of time averaged fields -> b1'
      write(6,*) 'time evolution of conservation laws -> c1'
      write(6,*) 'time evolution of electrons -> d1'
      write(6,*) 'time evolution of ions -> e1'
      write(6,*) 'time evolution of atoms -> f1'
      write(6,*) ' '
      write(6,*) 'Quit -> q'
      write(6,*) ' '

      read(5,*) char

      if (trim(char).eq.'a1') then
         call SELECT_pfield_evol
         goto 110   ! Go back to selection menue
      endif

      if (trim(char).eq.'b1') then
         call SELECT_tfield_evol
         goto 110   ! Go back to selection menue
      endif

      if (trim(char).eq.'c1') then
         call SELECT_cl_evol
         goto 110   ! Go back to selection menue
      endif

      if (trim(char).eq.'d1') then
         call SELECT_electron_evol
         goto 110   ! Go back to selection menue
      endif


      if (trim(char).eq.'e1') then
         call SELECT_ion_evol
         goto 110   ! Go back to selection menue
      endif


      if (trim(char).eq.'f1') then
         call SELECT_atom_evol
         goto 110   ! Go back to selection menue
      endif


      stop
      end
