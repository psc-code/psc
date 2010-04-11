c ====================================================================
c MONTE-CARLO COLLISIONAL MODULE BY A.KEMP AND H.RUHL 03/2005
c CONTAINS SPECIES (FOR BINARY COLLISIONS) AND MATERIAL (FOR ATOMIC PHYSICS) 
c LIST INFORMATION USED ONLY IN THIS MODULE. THIS MODULE OPERATES 
c INDEPENDENT OF THE BINARY COLLISION MODULE PIC_bin_coll.  
c INITIALIZATION OF ROUTINES IS DONE IN MCC_IMPACT VIA INIT_MCC.  
c ====================================================================

      module MCC_variables

      implicit none

      character*8, allocatable,dimension(:) :: matname            ! material names

      integer :: NCS,NPMAX,NMAT
      integer :: n_unsorted                                       ! number of unsorted particles

      integer,allocatable,dimension(:)     :: mcc_elist           ! separate list of electrons
      integer,allocatable,dimension(:,:,:) :: mcc_ilist           ! material particle index list
      integer,allocatable,dimension(:)     :: mcc_np              ! number of particles of material 
      integer,allocatable,dimension(:,:)   :: mcc_nc              ! number of particles of mat/charge

      integer,allocatable,dimension(:)     :: mcc_matlist         ! list of material numbers in random order
      integer,allocatable,dimension(:)     :: mcc_cslist          ! list of charge states in random order
      integer :: max_imat,max_ics

      integer, allocatable, dimension(:,:)        :: xstable_n      ! number of elements in table 
      integer, allocatable, dimension(:)          :: n_xstable      ! number of ii-xs-tables for mat
 
      integer :: p_0_toolarge
      integer :: n_xsect_tables
      integer :: p_0_count


      real(kind=8), allocatable, dimension(:,:)   :: xstable_t      ! ii energy threshold in eV 
      real(kind=8), allocatable, dimension(:,:,:) :: xstable_e      ! energy list for charge state in eV
      real(kind=8), allocatable, dimension(:,:,:) :: xstable_s      ! sigma  list for charge state in SI
      real(kind=8), allocatable, dimension(:,:)   :: max_sigmav     ! for charge state
      real(kind=8), allocatable, dimension(:)     :: mpart          ! particle mass list for matter id
      
      real(kind=8) :: me,dtsi
      real(kind=8) :: nel_pro,p_0_sum

      end module MCC_variables
