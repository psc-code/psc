      module PIC_variables
      implicit none


      real(kind=8), dimension(:), allocatable :: p_c1i
      real(kind=8), dimension(:), allocatable :: p_c2i
      real(kind=8), dimension(:), allocatable :: p_c3i
      real(kind=8), dimension(:), allocatable :: p_c4i

      real(kind=8) :: nudt,coll_add,coll_tot
      real(kind=8), dimension(:,:), allocatable :: nudt_min
      real(kind=8), dimension(:,:), allocatable :: nudt_max
      real(kind=8), dimension(:,:), allocatable :: nudt_mean
      real(kind=8), dimension(:,:), allocatable :: coll_count

      integer :: nicell

      integer :: nprpartio,nprparti
      integer :: dnprparti
      integer :: nistep

      real(kind=8) :: plin,plix

      end module PIC_variables

      
