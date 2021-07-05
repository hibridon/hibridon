! DESCRIPTION
!   This short program compares outputs from hibtest (in hibriddir/testnew)
!   with reference outputs (in hibriddir/tests)

module iso_fortran_env

    ! Nonintrinsic version for Lahey/Fujitsu Fortran for Linux. 
    ! See Subclause 13.8.2 of the Fortran 2003 standard. 
  
    implicit NONE 
    public 
  
    integer, parameter :: Character_Storage_Size = 8 
    integer, parameter :: Error_Unit = 0 
    integer, parameter :: File_Storage_Size = 8 
    integer, parameter :: Input_Unit = 5 
    integer, parameter :: IOSTAT_END = -1 
    integer, parameter :: IOSTAT_EOR = -2 
    integer, parameter :: Numeric_Storage_Size = 32 
    integer, parameter :: Output_Unit = 6 
  
end module iso_fortran_env

program comp_tests
    use iso_fortran_env, only: Error_Unit
    implicit none
    character(len=200) :: ref, test, line
    character(len=10) :: ext
    integer :: nstateso, nstates, it, i, j, ios
    real(8), dimension(:,:,:), allocatable :: xs
    real(8) :: dumr

    call get_command_argument(1,ref)
    call get_command_argument(2,test)

    ! Extract file extension in ext string variable
    j = scan(trim(ref),".", BACK= .true.)
    ext = ref(j+1:len(ref))

    if(ext=="xsc") then
        write(Error_Unit, *) "warning: unhandled file extension : ", ext
        stop 
    endif

    open(unit=1, status='old', file=trim(ref))
    open(unit=2, status='old', file=trim(test))


    if(ext=="evl") then 
        j=0
        ! READ REF AND TEST BOUND OUTPUT FILE
        ! First get the number of evls
        do 
            read(1,'(A)',iostat=ios) line 
            if(ios/=0) exit
            if(index(line,'** EIGENVALUES (CM-1)') /=0 ) then
                do; read(1,*,iostat=ios) dumr ; if(ios/=0) exit ; j=j+1 ; enddo
            endif
        enddo

        ! Allocate and initialize xs
        allocate(xs(2,j,1))
        rewind(1)
        ! Now read evls values from ref and test files
        do it = 1, 2
            do 
                read(it,'(A)',iostat=ios) line 
                if(ios/=0) exit
                if(index(line,'** EIGENVALUES (CM-1)') /=0 ) then
                    do i = 1, j ; read(it,*) xs(it,i,1) ; write(*,*) xs(it,i,1) ; enddo
                endif
            enddo
        enddo

    elseif(ext=="ics") then
        ! READ REF AND TEST ICSOUTPUT FILE
        do it = 1, 2
            do i = 1, 6 ; read(it,*) ; enddo ! Skip the first 6 lines
            read(it,*) nstates, nstateso ! Read the total number of states and open channels
            ! Allocate and initialize xs
            if(.not.allocated(xs)) then ; allocate( xs(2, nstateso, nstateso) ) ; xs = -1d0 ; endif
            j=ceiling(nstates/12.0d0)
            do i=1,j ; read(it,*) ; enddo ! Skip the next lines containing states labels
            j=ceiling(nstates/8.0d0)
            do i=1,j ; read(it,*) ; enddo ! Skip the next lines containing states energies
            ! Read the ICS
            do i = 1, nstateso
                read(it,*) xs(it,i,:)
            enddo
        enddo
    
    endif
    
    close(1)
    close(2)


    ! COMPARE XS

    do i=1,size(xs,2)
        do j=1,size(xs,3)
            ! stop with code 2 if difference is more than 1% for one of the XSs.
            if( 100d0 * abs(xs(2,i,j) - xs(1,i,j))/max(xs(1,i,j),1d-300) > 1d0 ) then
             stop 1
            endif
        enddo
    enddo
    if(allocated(xs)) then ; deallocate( xs ) ; endif
    return
end program comp_tests