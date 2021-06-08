! DESCRIPTION
!   This short program compares outputs from hibtest (in hibriddir/testnew)
!   with reference outputs (in hibriddir/tests)
  
program comp_tests
    implicit none
    character(len=200) :: ref, test
    integer :: nlines, ncols, it, i, j
    real(8), dimension(:,:,:), allocatable :: xs

    call get_command_argument(1,ref)
    call get_command_argument(2,test)

    open(unit=1, status='old', file=trim(ref))
    open(unit=2, status='old', file=trim(test))


! READ REF AND TEST ICS
do it = 1, 2
    do i = 1, 6 ; read(it,*) ; enddo ! Skip the first 6 lines
    read(it,*) ncols, nlines ! Read the number of lines and columns to read
    ! Allocate and initialize xs
    if(.not.allocated(xs)) then ; allocate( xs(2, nlines+1, ncols) ) ; xs = -1d0 ; endif
    read(it,*) ! Skip the next line
    ! Read the ICS
    do i = 1, nlines+1
        read(it,*) xs(it,i,:)
    enddo
enddo
close(1)
close(2)

! COMPARE XS

    do i=1,nlines+1
    do j=1,ncols
        ! stop with code 2 if difference is more than 1% for one of the XSs.
        if( 100d0 * abs(xs(2,i,j) - xs(1,i,j))/max(xs(1,i,j),1d-300) > 1d0 ) then
         stop 1
        endif
    enddo
    enddo
    return
end program comp_tests