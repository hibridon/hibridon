! DESCRIPTION
!   This short program compares outputs from hibtest (in hibriddir/testnew)
!   with reference outputs (in hibriddir/tests)

program comp_tests
    use iso_fortran_env, only: Error_Unit, Output_Unit
    use m_diff, only: result_files_differ
    implicit none
    character(len=200)  :: ref, test
    character(len=10)   :: ext
    real(8)             :: tolerance
    real(8)             :: min_significant_value ! values below this are ignored in the comparison
    integer             :: num_header_lines(2)
    ! Get test and reference files paths from command arguments
    call get_command_argument(1,ref)
    call get_command_argument(2,test)

    ! Extract file extension in ext string variable
    ext = ref(scan(trim(ref),".", BACK= .true.)+1:len(ref))

    write(Output_Unit, "(4a)") 'comparing test file ', trim(test), ' to reference file ', trim(ref)

    ! Tolerance is set to 1%
    tolerance=0.01d0

    min_significant_value = 1e-30

    ! Set the number of header lines depending on the type of output file
    num_header_lines = 0
    select case (ext)
    case("ics") ; num_header_lines = 3
    case("dcs") ; num_header_lines = 15
    case("hfx") ; num_header_lines = 6
    case("ppb") ; num_header_lines = 6
    case("pcs") ; num_header_lines = 6
        min_significant_value = 1e-20
    case("trn") ; num_header_lines = 7
    case("xxsc"); num_header_lines = 3
    case("tcb"); num_header_lines = 2
    case ("flx") ! Header ends at first occurence of "R (BOHR) AND OUTGOING FLUXES"
        num_header_lines(1) = get_first_occ_of("R (BOHR) AND OUTGOING FLUXES",ref)
        num_header_lines(2) = get_first_occ_of("R (BOHR) AND OUTGOING FLUXES",test)
        ! the values of the flx files can be highly sensitive to compilers (see doc/hib_html/tests.html and issue #37)
        ! so we ignore more values
        min_significant_value = 1e-10
    case("evl") ! Header ends at first occurence of "** EIGENVALUES"
        num_header_lines(1) = get_first_occ_of("** EIGENVALUES",ref)
        num_header_lines(2) = get_first_occ_of("** EIGENVALUES",test)
    case("xsc") ! Header ends at first occurence of "INTEGRAL CROSS SECTIONS"
        num_header_lines(1) = get_first_occ_of("INTEGRAL CROSS SECTIONS",ref)
        num_header_lines(2) = get_first_occ_of("INTEGRAL CROSS SECTIONS",test)
    case("stdout") ! Header ends at first occurence of "Hibridon>"
        num_header_lines(1) = get_first_occ_of("Hibridon>",ref)
        num_header_lines(2) = get_first_occ_of("Hibridon>",test)
    case("out")
        num_header_lines(1) = get_first_occ_of("** LABEL:",ref)
        num_header_lines(2) = get_first_occ_of("** LABEL:",test)
    end select

    ! Compare numeric values between reference and test files
    if(result_files_differ(ref, test, num_header_lines, tolerance, min_significant_value)) stop 1


contains 

! This function returns the line where the first occurence of a substring is find.
! If the substring is not detected, STOP 1.
function get_first_occ_of(sub,file) result(iline)
    implicit none
    character(len=*), intent(in) :: sub, file
    integer :: iline, ierr
    character(len=200) :: line
    iline = 0
    open(unit=1, file=file, status='old')
    do
        read(1,'(a)', iostat=ierr) line ; if(ierr.ne.0) exit ! Try to read the current line. If it fails, exit the do loop.
        iline = iline+1 ! Increment the line number
        if (index(line, trim(sub)).ne.0) then ; close(1) ; return ; endif ! If the substring is found, close the file and return.
    enddo
    close(1)
    STOP 1 ! If the function did not return within the do loop, then the substring is absent -> STOP 1.
end function get_first_occ_of





end program comp_tests
