#include "assert.h"
!#define TEST_V2MAT_USE_ASSOCIATE

program test_ancou_type
    use mod_cov2, only: ancou_type, ancouma_type
    implicit none
    integer :: nlam = 1000
    integer :: num_channels = 200
    integer :: ilam
    integer :: irow, icol
    integer :: nloops = 1
    integer :: iloop
    integer(8) :: iv2_element, num_nz_elements
    integer :: ij
    real(8) :: vee
    real(8) :: sum
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
       do irow = 1, num_channels
          do icol = 1, num_channels
             call ancouma%set_element(irow=irow, icol=icol, vee=42.d0)
             !call v2%set_element(ilam=ilam, irow=irow, icol=icol, vee=42.d0)
          end do
       end do
    end do

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
end program test_ancou_type
