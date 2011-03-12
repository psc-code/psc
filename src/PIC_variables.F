      module PIC_variables
      implicit none


      real(kind=8), dimension(:), allocatable :: p_niloc
      real(kind=8), dimension(:), allocatable :: p_e1i
      real(kind=8), dimension(:), allocatable :: p_e2i
      real(kind=8), dimension(:), allocatable :: p_e3i
      real(kind=8), dimension(:), allocatable :: p_e4i
      real(kind=8), dimension(:), allocatable :: p_e5i
      real(kind=8), dimension(:), allocatable :: p_e6i
      real(kind=8), dimension(:), allocatable :: p_c1i
      real(kind=8), dimension(:), allocatable :: p_c2i
      real(kind=8), dimension(:), allocatable :: p_c3i
      real(kind=8), dimension(:), allocatable :: p_c4i
      real(kind=8), dimension(:), allocatable :: p_re1i
      real(kind=8), dimension(:), allocatable :: p_re2i
      real(kind=8), dimension(:), allocatable :: p_re3i
      real(kind=8), dimension(:), allocatable :: p_re4i
      real(kind=8), dimension(:), allocatable :: p_re5i
      real(kind=8), dimension(:), allocatable :: p_re6i

      real(kind=8) :: nudt,coll_add,coll_tot
      real(kind=8), dimension(:,:), allocatable :: nudt_min
      real(kind=8), dimension(:,:), allocatable :: nudt_max
      real(kind=8), dimension(:,:), allocatable :: nudt_mean
      real(kind=8), dimension(:,:), allocatable :: coll_count

      integer :: e1i,e2i,e3i,e4i,e5i,e6i
      integer :: re1i,re2i,re3i,re4i,re5i,re6i

      integer :: niloc
      integer :: niloc_n
      integer :: nicell
      integer :: nitot
      integer :: nialloc

      integer :: nprpartio,nprparti
      integer :: dnprparti
      integer :: enprparti
      integer :: nistep

      real(kind=8) :: plin,plix,pli

      end module PIC_variables

      
