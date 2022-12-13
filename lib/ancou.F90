#include "assert.h"
! module dealing with angular coupling matrices storage
module mod_ancou

   use mod_grovec, only: dgrovec_type, igrovec_type
   implicit none
   integer(8), public, protected :: g_num_ancou_set_calls = 0
   integer(8), public, protected :: g_num_ancou_get_calls = 0
   integer(8), public, protected :: g_num_ancouma_set_calls = 0
   integer(8), public, protected :: g_num_ancouma_get_calls = 0
   integer, public, protected :: g_num_ancou_instances = 0
   integer, public, protected :: g_num_ancouma_instances = 0
   ! ancouma_type stores the angular coupling matrix related to a singe lambda
   ! lower triangle of nonzero angular coupling matrix elements
   ! stored in packed column form that is :
   !      (1,1), (2,1), (3,1) ... (n,1),
   !             (2,2), (3,2) ... (n,2), etc.
   ! only nonzero elements are stored
   type                         :: ancouma_type
     logical                    :: is_allocated = .false.
     integer                    :: num_channels = 0
     type(dgrovec_type), allocatable :: v2d  ! the non zero element values
     type(igrovec_type), allocatable :: v2i  ! the non-zero elements indices 
     contains
     
     procedure                  :: get_num_nonzero_elements => ancouma_type_get_num_nonzero_elements
     procedure                  :: get_element => ancouma_type_get_element
     procedure                  :: set_element => ancouma_type_set_element
     procedure                  :: get_allocated_mem => ancouma_type_get_allocated_mem
     procedure                  :: get_used_mem => ancouma_type_get_used_mem
   end type ancouma_type

   ! ancou_type stores the angular coupling matrix elements for all lamdbas
   type, public                 :: ancou_type
     integer                    :: nlam = 0
     integer                    :: num_channels = 0
     type(ancouma_type), allocatable  :: ancouma(:)
     contains
     final                      :: ancou_type_destroy
     procedure                  :: set_element => ancou_type_set_element
     procedure                  :: get_element => ancou_type_get_element
     procedure                  :: get_num_nonzero_elements => ancou_type_get_num_nonzero_elements
     procedure                  :: get_angular_coupling_matrix => ancou_type_get_angular_coupling_matrix
     procedure                  :: ensure_ancouma_is_allocated => ancou_type_ensure_ancouma_is_allocated
     procedure                  :: empty => ancou_type_empty
     procedure                  :: print => ancou_type_print
     procedure                  :: print_summary => ancou_type_print_summary

   end type ancou_type

   ! interface to create an instance of ancou_type using a construct familiar to other languages:
   ! g1 = ancou_type()
   interface ancou_type
     module procedure create_ancou_type
   end interface ancou_type

   contains

   !
   ! ancouma_type implementation
   !

   function ancouma_type_get_num_nonzero_elements(this) result(n)
      class(ancouma_type)        :: this
      integer           :: n
      n = this%v2d%num_elements
   end function 

   subroutine ancouma_type_get_element(this, ielement, ij, vee)
      class(ancouma_type)        :: this
      integer, intent(in) :: ielement
      integer, intent(out) :: ij
      real(8), intent(out) :: vee
      ij = this%v2i%get_element(ielement)
      vee = this%v2d%get_element(ielement)
      g_num_ancouma_get_calls = g_num_ancouma_get_calls + 1
   end subroutine 

   subroutine ancouma_type_set_element(this, irow, icol, vee)
      class(ancouma_type)       :: this
      integer, intent(in) :: irow
      integer, intent(in) :: icol
      real(8), intent(in) :: vee
      integer :: ij
      ij = this%num_channels * (icol - 1) +irow
      call this%v2d%append(vee)
      call this%v2i%append(ij)
      g_num_ancouma_set_calls = g_num_ancouma_set_calls + 1
   end subroutine

   ! returns the memory allocated by this instance of ancouma_type
   function ancouma_type_get_allocated_mem(this) result(num_bytes)
      class(ancouma_type) :: this
      integer             :: num_bytes

      integer, parameter :: size_of_int = 4
      integer, parameter :: size_of_double = 8
      num_bytes = this%v2d%num_allocated_blocks * this%v2d%block_size * size_of_double & 
                + this%v2i%num_allocated_blocks * this%v2i%block_size * size_of_int
   end function

   ! returns the memory used by this instance of ancouma_type
   function ancouma_type_get_used_mem(this) result(num_bytes)
      class(ancouma_type) :: this
      integer             :: num_bytes

      integer, parameter :: size_of_int = 4
      integer, parameter :: size_of_double = 8
      num_bytes = this%get_num_nonzero_elements() * (size_of_int + size_of_double)
   end function

   !
   ! ancou_type implementation
   !

   function create_ancou_type(nlam, num_channels) result(v2)
     integer, intent(in) :: nlam
     integer, intent(in) :: num_channels

     type(ancou_type)  :: v2
     integer :: ilam

     v2%nlam = nlam
     v2%num_channels = num_channels
     allocate(v2%ancouma(nlam))
     do ilam = 1, nlam
       v2%ancouma(ilam)%is_allocated = .false.
     !  v2%ancouma(ilam)%num_channels = num_channels
     !  v2%ancouma(ilam)%v2d = dgrovec_type(block_size=1024*4, num_blocks=1024*8)
     !  v2%ancouma(ilam)%v2i = igrovec_type(block_size=1024*4, num_blocks=1024*8)
     end do
     g_num_ancou_instances = g_num_ancou_instances + 1
   end function

   subroutine ancou_type_destroy(this)
      type(ancou_type)        :: this
      if (allocated(this%ancouma)) then
         deallocate(this%ancouma)
      end if
   end subroutine

   function ancou_type_get_num_nonzero_elements(this) result(num_nz_elements)
      class(ancou_type), target  :: this
      integer :: ilam
      integer :: lam_nz_els, num_nz_elements
      type(ancouma_type), pointer :: ancouma 
      do ilam = 1, this%nlam
         ancouma => this%get_angular_coupling_matrix(ilam)
         lam_nz_els = ancouma%get_num_nonzero_elements()
         num_nz_elements = num_nz_elements + lam_nz_els
      end do
   end function

   function ancou_type_get_angular_coupling_matrix(this, ilam) result(ancouma)
      class(ancou_type), target  :: this
      integer, intent(in) :: ilam
      type(ancouma_type), pointer :: ancouma
      call this%ensure_ancouma_is_allocated(ilam)
      ancouma => this%ancouma(ilam)
   end function

   subroutine ancou_type_ensure_ancouma_is_allocated(this, ilam)
      class(ancou_type)        :: this
      integer, intent(in) :: ilam
      integer :: block_size
      integer :: num_blocks
      integer :: max_num_elements
      real(8) :: expected_sparsity = 0.1
      integer :: expected_num_elements
      ASSERT(ilam > 0)
      ASSERT(ilam-1 <= this%nlam)
      ! v2 matrix is triangular
      max_num_elements = this%num_channels * (this%num_channels + 1)
      expected_num_elements = expected_sparsity * max_num_elements
      ! optimize block size so that one block should be enough in most cases
      block_size = expected_num_elements
      num_blocks = (max_num_elements / block_size) + 1
      if (.not. this%ancouma(ilam)%is_allocated) then
         allocate(this%ancouma(ilam)%v2d)
         allocate(this%ancouma(ilam)%v2i)
         this%ancouma(ilam)%v2d = dgrovec_type(block_size=block_size, num_blocks=num_blocks)
         this%ancouma(ilam)%v2i = igrovec_type(block_size=block_size, num_blocks=num_blocks)
         this%ancouma(ilam)%num_channels = this%num_channels
         this%ancouma(ilam)%is_allocated = .true.
         g_num_ancouma_instances = g_num_ancouma_instances + 1
      end if
    end subroutine

   subroutine ancou_type_set_element(this, ilam, irow, icol, vee)
      implicit none
      class(ancou_type), target :: this
      integer, intent(in) :: ilam
      integer, intent(in) :: irow
      integer, intent(in) :: icol
      real(8), intent(in) :: vee
      type(ancouma_type), pointer :: ancouma
      ancouma => this%get_angular_coupling_matrix(ilam)
      call ancouma%set_element(irow, icol, vee)
      g_num_ancou_set_calls = g_num_ancou_set_calls + 1
   end subroutine 

   subroutine ancou_type_get_element(this, ilam, ielement, ij, vee)
      implicit none
      class(ancou_type)        :: this
      integer, intent(in) :: ilam
      integer, intent(in) :: ielement
      integer, intent(out) :: ij
      real(8), intent(out) :: vee

      call this%ancouma(ilam)%get_element(ielement, ij, vee)
      g_num_ancou_get_calls = g_num_ancou_get_calls + 1
   end subroutine 

   subroutine ancou_type_empty(this)
      class(ancou_type)        :: this
      integer :: ilam
      do ilam = 1, this%nlam
         call this%ancouma(ilam)%v2d%empty()
         call this%ancouma(ilam)%v2i%empty()
      end do
   end subroutine 

   ! prints the contents of this ancou_type
   subroutine ancou_type_print(this, unit)
      class(ancou_type), intent(in)    :: this
      integer, intent(in) :: unit
      integer :: ilam, num_lam_nz_elements, num_channels, lam_num_elements
      integer :: num_nz_elements, num_elements, nz_el_index, irow, icol
      integer :: ij
      real(8) :: vee
      num_nz_elements = int(0, 8)
      num_elements = int(0, 8)
      do ilam = 1, this%nlam
         if (this%ancouma(ilam)%is_allocated) then
            num_channels = this%ancouma(ilam)%num_channels
            lam_num_elements = num_channels * num_channels
            num_lam_nz_elements = this%ancouma(ilam)%get_num_nonzero_elements()
            if (num_lam_nz_elements > 0) then
               write (unit,'(A, I3, A, I10, A, I10, A, F6.2, A)') 'ilam = ', ilam, ' : ', num_lam_nz_elements, '/', lam_num_elements ,' non zero elements (', real(num_lam_nz_elements, 8) / lam_num_elements * 100.d0, '%)'
               num_nz_elements = num_nz_elements + num_lam_nz_elements
            end if
            do nz_el_index = 1, num_lam_nz_elements
               call this%ancouma(ilam)%get_element(nz_el_index, ij, vee)
               icol = ij / num_channels
               irow = modulo(ij, num_channels)
               write(unit, fmt='(A, I4, A, I4, A, F6.2, A)') '(', irow, ', ', icol, ')=', vee, ' '
            end do
            num_elements = num_elements + lam_num_elements
         end if
      end do
      write (unit,'(A, I10, A, I10, A, F6.2, A)') 'total : ', num_nz_elements, '/', num_elements ,' non zero elements (', real(num_nz_elements, 8) / num_elements * 100.d0, '%)'
   end subroutine 

   ! prints a summary of this instance of ancou_type
   subroutine ancou_type_print_summary(this, unit)
      class(ancou_type), intent(in)    :: this
      integer, intent(in) :: unit
      integer :: ilam, num_lam_nz_elements, num_channels, lam_num_elements
      integer :: num_nz_elements, num_elements
      integer(8) :: allocated_memory  ! in bytes
      integer(8) :: used_memory  ! in bytes
      allocated_memory = 0
      used_memory = 0
      num_nz_elements = int(0, 8)
      num_elements = int(0, 8)
      do ilam = 1, this%nlam
         if (this%ancouma(ilam)%is_allocated) then
            allocated_memory = allocated_memory + this%ancouma(ilam)%get_allocated_mem()
            used_memory = used_memory + this%ancouma(ilam)%get_used_mem()
            num_channels = this%ancouma(ilam)%num_channels
            lam_num_elements = num_channels * num_channels
            num_lam_nz_elements = this%ancouma(ilam)%get_num_nonzero_elements()
            if (num_lam_nz_elements > 0) then
               write (unit,'(A, I3, A, I10, A, I10, A, F6.2, A)') 'ilam = ', ilam, ' : ', num_lam_nz_elements, '/', lam_num_elements ,' non zero elements (sparsity: ', real(num_lam_nz_elements, 8) / lam_num_elements * 100.d0, '%)'
               num_nz_elements = num_nz_elements + num_lam_nz_elements
            end if
            num_elements = num_elements + lam_num_elements
         end if
      end do
      write (unit,'(A, I10, A, I10, A, F6.2, A)') 'total : ', num_nz_elements, '/', num_elements ,' non zero elements (sparsity :', real(num_nz_elements, 8) / num_elements * 100.d0, '%)'
      write (unit,'(A, F6.2, A, I15, A, I15, A)') 'storage efficiency : ', real(used_memory, 8) / real(allocated_memory, 8) * 100.d0,'% (', used_memory,' used bytes / ', allocated_memory,' allocated bytes)'
      
   end subroutine 

   ! end of ancou_type implmentation

   subroutine print_ancou_stats(unit)
      implicit none
      integer, intent(in) :: unit
      write(unit,*) "number of calls to ancou_type%set_element : ", g_num_ancou_set_calls
      write(unit,*) "number of calls to ancou_type%get_element : ", g_num_ancou_get_calls
      write(unit,*) "number of calls to ancouma_type%set_element : ", g_num_ancouma_set_calls
      write(unit,*) "number of calls to ancouma_type%get_element : ", g_num_ancouma_get_calls
      write(unit,*) "number of ancou_type instances created : ", g_num_ancou_instances
      write(unit,*) "number of ancouma_type instances created : ", g_num_ancouma_instances
   end subroutine

end module mod_ancou