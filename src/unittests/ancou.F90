#include "assert.h"
!#define TEST_V2MAT_USE_ASSOCIATE

#if defined(HIB_UNIX_X86)
#define SYSTEM_MEM_USAGE_WORKS
#endif

#if defined(SYSTEM_MEM_USAGE_WORKS)
subroutine system_mem_usage(valueRSS)
implicit none
!use ifport !if on intel compiler
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

program test_ancou_type
    use mod_cov2, only: ancou_type, ancouma_type, print_ancou_stats
    use mod_grovec, only: print_grovec_stats
    implicit none
    integer :: nlam = 80
    integer :: num_channels = 70
    integer :: ilam
    integer :: irow, icol
    integer :: nloops = 6
    integer :: nread_loops = 200
    integer :: iloop, iread_loop
    integer :: iv2_element, num_nz_elements
    integer :: num_elements
    integer :: ij
#if defined(SYSTEM_MEM_USAGE_WORKS)
    integer :: mem_used_by_this_process
#endif
    real(8) :: vee
    real(8) :: sum
    real(8) :: fill_ratio = 0.18
#ifdef TEST_V2MAT_USE_ASSOCIATE
    type(ancou_type) :: v2
#else
    type(ancou_type) :: v2, target
    type(ancouma_type), pointer :: ancouma
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
             do iv2_element = 1, num_nz_elements
               call ancouma%get_element(iv2_element, ij, vee)
               sum = sum + vee
             end do
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
   call print_grovec_stats()
   call print_ancou_stats()
end program test_ancou_type
