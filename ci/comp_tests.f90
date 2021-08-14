! DESCRIPTION
!   This short program compares outputs from hibtest (in hibriddir/testnew)
!   with reference outputs (in hibriddir/tests)

program comp_tests
    use iso_fortran_env, only: Error_Unit, Output_Unit
    use m_diff, only: result_files_differ
    implicit none
    character(len=200) :: ref, test
    character(len=10) :: ext
    real(8) :: tolerance
    integer :: num_header_lines(2)
    ! Get test and reference files paths from command arguments
    call get_command_argument(1,ref)
    call get_command_argument(2,test)

    ! Extract file extension in ext string variable
    ext = ref(scan(trim(ref),".", BACK= .true.)+1:len(ref))

    write(Output_Unit, *) 'comparing test file ', trim(test), ' to reference file ', trim(ref)

    ! Tolerance is set to 1%
    tolerance=0.01d0 

    select case (ext)
    case("ics") ; num_header_lines = 3
    case("dcs") ; num_header_lines = 15
    case("hfx") ; num_header_lines = 6
    case("ppb") ; num_header_lines = 6
    case("pcs") ; num_header_lines = 6
    case("trn") ; num_header_lines = 7
    case("xxsc"); num_header_lines = 3
    case ("flx") ! Header ends at first occurence of "R (BOHR) AND OUTGOING FLUXES"
        num_header_lines(1) = get_first_occ_of("R (BOHR) AND OUTGOING FLUXES",ref)
        num_header_lines(2) = get_first_occ_of("R (BOHR) AND OUTGOING FLUXES",test)
    case("evl") ! Header ends at first occurence of "** EIGENVALUES"
        num_header_lines(1) = get_first_occ_of("** EIGENVALUES",ref)
        num_header_lines(2) = get_first_occ_of("** EIGENVALUES",test)
    case("xsc") ! Header ends at first occurence of "INTEGRAL CROSS SECTIONS"
        num_header_lines(1) = get_first_occ_of("INTEGRAL CROSS SECTIONS",ref)
        num_header_lines(2) = get_first_occ_of("INTEGRAL CROSS SECTIONS",test)
    case("stdout") ! Header ends at first occurence of "Hibridon>"
        num_header_lines(1) = get_first_occ_of("Hibridon>",ref)
        num_header_lines(2) = get_first_occ_of("Hibridon>",test)
    end select


    if(result_files_differ(ref, test, num_header_lines, tolerance)) stop 1


contains 

function get_first_occ_of(sub,file) result(iline)
    implicit none
    character(len=*), intent(in) :: sub, file
    integer :: iline
    character(len=200) :: line
    integer :: ierr, i

    iline=-1
    i=0
    open(unit=1, file=file, status='old')
    do
        read(1,'(a)', iostat=ierr) line
        if(ierr.ne.0) exit
        i = i+1
        if (index(line, trim(sub)).ne.0) then
            iline=i ; exit
        endif
    enddo
    if(iline==-1) STOP 1
    close(1)
end function get_first_occ_of





end program comp_tests
