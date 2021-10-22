#include "assert.h"
!#define TEST_V2MAT_USE_ASSOCIATE

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
    real(8) :: vee
    real(8) :: sum
    real(8) :: fill_ratio = 0.18
#ifdef TEST_V2MAT_USE_ASSOCIATE
    type(ancou_type) :: v2
#else
    type(ancou_type) :: v2, target
    type(ancouma_type), pointer :: ancouma
#endif
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
   end do
   call print_grovec_stats()
   call print_ancou_stats()
end program test_ancou_type
