#include "assert.h"
!#define TEST_V2MAT_USE_ASSOCIATE

#if defined(HIB_UNIX_IFORT) || defined(HIB_UNIX_X86)
#define SYSTEM_MEM_USAGE_WORKS
#endif

#if defined(SYSTEM_MEM_USAGE_WORKS)
subroutine system_mem_usage(valueRSS)
#if defined(HIB_UNIX_IFORT)
use ifport ! needed to use getpid on ifort compiler
#endif 
implicit none
integer, intent(out) :: valueRSS

character(len=200):: filename=' '
character(len=80) :: line
character(len=8)  :: pid_char=' '
integer :: pid
logical :: ifxst

valueRSS=-1    ! return negative number if not found

!--- get process ID

pid=getpid()
write(pid_char,'(I8)') pid
filename='/proc/'//trim(adjustl(pid_char))//'/status'

!--- read system file

inquire (file=filename,exist=ifxst)
if (.not.ifxst) then
  write (*,*) 'system file does not exist'
  return
endif

open(unit=100, file=filename, action='read')
do
  read (100,'(a)',end=120) line
  if (line(1:6).eq.'VmRSS:') then
     read (line(7:),*) valueRSS
     exit
  endif
enddo
120 continue
close(100)

return
end subroutine system_mem_usage
#endif


#define ANCOUMA_READ_METHOD_NORMAL 0
#define ANCOUMA_READ_METHOD_INLINE_LEVEL1 1
#define ANCOUMA_READ_METHOD_INLINE_LEVEL2 2
#define ANCOUMA_READ_METHOD_INLINE_LEVEL3 3

#define ANCOUMA_READ_METHOD ANCOUMA_READ_METHOD_NORMAL


#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL1)

#define ANCOUMA_GET_ELEMENT(ancouma, element_index, ij, vee) \
ij = ancouma%v2i%get_element(iv2_element) ; \
vee = ancouma%v2d%get_element(iv2_element)

#endif


subroutine test_alloc_simpler
   implicit none
   integer :: i, block_size
   integer :: nloops = 3
   integer :: iloop
#if defined(SYSTEM_MEM_USAGE_WORKS)
   integer :: mem_used_by_this_process
#endif
   integer, allocatable :: b3(:)

   block_size = 10000000
   do iloop = 1, nloops
      if (allocated(b3)) then
       write(6, *) 'deallocating block'
       deallocate(b3)
      end if
#if defined(SYSTEM_MEM_USAGE_WORKS)
      call system_mem_usage(mem_used_by_this_process)
      write (6,*) 'mem used by this process : ', mem_used_by_this_process, ' kbytes'
#endif
      write(6, *) 'allocating block'
      allocate(b3(block_size))
      do i = 1, block_size
         b3(i) = 42
      end do
      block_size = block_size - 1000000
#if defined(SYSTEM_MEM_USAGE_WORKS)
      call system_mem_usage(mem_used_by_this_process)
      write (6,*) 'mem used by this process : ', mem_used_by_this_process, ' kbytes'
#endif
   end do
end subroutine test_alloc_simpler

! ancou statistics from test nh3h2_qma_long
!  total :   24452296/ 834406375 non zero elements (sparsity :  2.93%)
!  storage efficiency :  29.30% (      293427552 used bytes /      1001544720 allocated bytes)
!   number of calls to ancou_type%set_element :                     0
!   number of calls to ancou_type%get_element :                     0
!   number of calls to ancouma_type%set_element :              24452296
!   number of calls to ancouma_type%get_element :                     0
!   number of ancou_type instances created :            1
!   number of ancouma_type instances created :           55

program test_ancou_type
    use mod_ancou, only: ancou_type, ancouma_type, print_ancou_stats
    use mod_grovec, only: print_grovec_stats, igrovec_type_block, dgrovec_type_block

    implicit none
    integer :: nlam = 55  ! same as nh3h2_qma_long test
    integer :: num_channels = 3895  ! same as nh3h2_qma_long test
    integer :: ilam
    integer :: irow, icol
    integer :: nloops = 1
    integer :: nread_loops = 120  ! generates roughly the same number of calls to  ancouma_type%get_element (3003864600) than nh3h2_qma_long test (3032084704)
    integer :: iloop, iread_loop
    integer :: iv2_element, num_nz_elements
    integer :: num_elements
    integer :: ij
#if defined(SYSTEM_MEM_USAGE_WORKS)
    integer :: mem_used_by_this_process
#endif
    real(8) :: vee
    real(8) :: sum
    real(8) :: fill_ratio = 0.03  ! roughly the same as nh3h2_qma_long test (2.93%)
#ifdef TEST_V2MAT_USE_ASSOCIATE
    type(ancou_type) :: v2
#else
    type(ancou_type) :: v2, target
    type(ancouma_type), pointer :: ancouma
#endif
#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL2) || (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL3)
integer :: block_size
integer :: block_index
integer :: el_index_in_block
#endif

#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL3)
integer :: num_full_blocks
integer :: num_remaining_nz_elements
type(igrovec_type_block), pointer :: blocki
type(dgrovec_type_block), pointer :: blockd
#endif
   !call test_alloc_simpler()
   !stop

    do iloop = 1, nloops
      v2 = ancou_type(nlam = nlam, num_channels=num_channels)
       do ilam = 1, v2%nlam
          ancouma => v2%get_angular_coupling_matrix(ilam)
          num_elements = num_channels * num_channels
          num_nz_elements = int(fill_ratio * num_elements, 8)
          ij = 0
          do irow = 1, num_channels
             do icol = 1, num_channels
                call ancouma%set_element(irow=irow, icol=icol, vee=42.d0)
                !call v2%set_element(ilam=ilam, irow=irow, icol=icol, vee=42.d0)
                ij = ij + 1
                if (ij > num_nz_elements) then
                  exit
                end if
             end do
             if (ij > num_nz_elements) then
               exit
             end if
          end do
       end do

       do iread_loop =  1, nread_loops
          sum = 0.d0
          do ilam = 1, v2%nlam
             ! write(6,*) 'ilam=', ilam, 'v2%get_angular_coupling_matrix(ilam)%get_num_nonzero_elements()=', v2%get_angular_coupling_matrix(ilam)%get_num_nonzero_elements()
             !ancouma => v2%ancouma(ilam)
#ifdef TEST_V2MAT_USE_ASSOCIATE
             associate( ancouma => v2%get_angular_coupling_matrix(ilam) )
#else
             ancouma => v2%get_angular_coupling_matrix(ilam)
#endif
             num_nz_elements = ancouma%get_num_nonzero_elements()
             ASSERT(num_nz_elements >= 0)

#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_NORMAL )
             do iv2_element = 1, num_nz_elements
               call ancouma%get_element(iv2_element, ij, vee)  
               sum = sum + vee
             end do
#endif

#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL1)
             do iv2_element = 1, num_nz_elements
               ANCOUMA_GET_ELEMENT(ancouma, iv2_element, ij, vee)
               sum = sum + vee
             end do
#endif

#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL2)
             block_size = ancouma%v2i%block_size
             do iv2_element = 1, num_nz_elements
               block_index =  (iv2_element-1) / block_size
               el_index_in_block = iv2_element - 1 - (block_index * block_size)
               ij = ancouma%v2i%blocks(block_index)%p(el_index_in_block)
               vee = ancouma%v2d%blocks(block_index)%p(el_index_in_block)
               sum = sum + vee
             end do
#endif

#if (ANCOUMA_READ_METHOD == ANCOUMA_READ_METHOD_INLINE_LEVEL3)
             block_size = ancouma%v2i%block_size
             num_full_blocks = ((num_nz_elements - 1) / block_size)
             iv2_element = 0
             do block_index = 0, num_full_blocks-1
               blocki => ancouma%v2i%blocks(block_index)
               blockd => ancouma%v2d%blocks(block_index)
               do el_index_in_block = 0, block_size-1
                 ij = blocki%p(el_index_in_block)
                 vee = blockd%p(el_index_in_block)
                 iv2_element = iv2_element + 1
                 sum = sum + vee
               end do
             end do
             num_remaining_nz_elements = num_nz_elements - iv2_element
             if (num_remaining_nz_elements > 0) then
               blocki => ancouma%v2i%blocks(num_full_blocks)
               blockd => ancouma%v2d%blocks(num_full_blocks)
               do el_index_in_block = 0, num_remaining_nz_elements-1
                 ij = blocki%p(el_index_in_block)
                 vee = blockd%p(el_index_in_block)
                 iv2_element = iv2_element + 1
                 sum = sum + vee
               end do
             end if
#endif

#ifdef TEST_V2MAT_USE_ASSOCIATE
             end associate
#endif
         end do
      end do
      call v2%print_summary(unit=6)
#if defined(SYSTEM_MEM_USAGE_WORKS)
      call system_mem_usage(mem_used_by_this_process)
      write (6,*) 'mem used by this process : ', mem_used_by_this_process, ' kbytes'
#endif      
   end do
   call print_grovec_stats(unit=6)
   call print_ancou_stats(unit=6)
end program test_ancou_type
