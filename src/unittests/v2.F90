!#define TEST_V2MAT_USE_ASSOCIATE

program test_v2mat
    use mod_cov2, only: v2mat, lamv2t
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
    type(v2mat) :: v2
#else
    type(v2mat) :: v2, target
    type(lamv2t), pointer :: lamv2
#endif
    do iloop = 1, nloops
    v2 = v2mat(nlam = nlam, num_channels=num_channels)
    do ilam = 1, v2%nlam
       lamv2 => v2%get_lamv2(ilam)
       do irow = 1, num_channels
          do icol = 1, num_channels
             call lamv2%set_element(irow=irow, icol=icol, vee=42.d0)
             !call v2%set_element(ilam=ilam, irow=irow, icol=icol, vee=42.d0)
          end do
       end do
    end do

    sum = 0.d0
    do ilam = 1, v2%nlam
       ! write(6,*) 'ilam=', ilam, 'v2%get_lamv2(ilam)%get_num_elements()=', v2%get_lamv2(ilam)%get_num_elements()
       !lamv2 => v2%lamv2(ilam)
#ifdef TEST_V2MAT_USE_ASSOCIATE
       associate( lamv2 => v2%get_lamv2(ilam) )
#else
       lamv2 => v2%get_lamv2(ilam)
#endif
       num_nz_elements = lamv2%get_num_elements()
       ASSERT(num_nz_elements >= 0)
       do iv2_element = 1, num_nz_elements
         call lamv2%get_element(iv2_element, ij, vee)
         sum = sum + vee
       end do
#ifdef TEST_V2MAT_USE_ASSOCIATE
       end associate
#endif       
   end do
   end do
end program test_v2mat
