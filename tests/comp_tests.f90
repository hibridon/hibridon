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
    character(len=200)  :: pattern
    integer :: num_patterns
    integer :: pattern_index
    logical :: pattern_is_found
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
        num_patterns = 2
        pattern_is_found = .false.
        do pattern_index=1, num_patterns
            if (pattern_index == 1) then
                pattern = "R (BOHR) AND OUTGOING FLUXES"
            else if (pattern_index == 2) then
                pattern = "TRANSFORMATION MATRIX AT R = "
            else
                write(Error_Unit, "(a)") "Error: Invalid pattern index."
                stop 1
            end if
            write(Output_Unit, "(a)") "testing against pattern: " // trim(pattern)

            call get_header_lines_counts(ref, test, pattern, num_header_lines(1), num_header_lines(2), pattern_is_found)
            if (pattern_is_found) then
                exit ! Exit the loop if the pattern is found
            end if
        end do

        if (.not. pattern_is_found) then
            write(Error_Unit, "(a)") "Error: Pattern not found in either reference or test file."
            stop 1
        end if

        ! the values of the flx files can be highly sensitive to compilers (see doc/hib_html/tests.html and issue #37)
        ! so we ignore more values
        min_significant_value = 1e-10
    case("evl") ! Header ends at first occurence of "** EIGENVALUES"
        pattern = "** EIGENVALUES"
        call get_header_lines_counts(ref, test, pattern, num_header_lines(1), num_header_lines(2), pattern_is_found)
    case("xsc") ! Header ends at first occurence of "INTEGRAL CROSS SECTIONS"
        pattern = "INTEGRAL CROSS SECTIONS"
        call get_header_lines_counts(ref, test, pattern, num_header_lines(1), num_header_lines(2), pattern_is_found)
    case("stdout") ! Header ends at first occurence of "Hibridon>"
        pattern = "Hibridon>"
        call get_header_lines_counts(ref, test, pattern, num_header_lines(1), num_header_lines(2), pattern_is_found)
    case("out")
        pattern = "** LABEL:"
        call get_header_lines_counts(ref, test, pattern, num_header_lines(1), num_header_lines(2), pattern_is_found)
    end select

    ! Compare numeric values between reference and test files
    if(result_files_differ(ref, test, num_header_lines, tolerance, min_significant_value)) stop 1


contains 

! This function returns the line where the first occurence of a pattern is find.
! If the pattern is not detected, pattern_is_found is set to .false
function get_first_occ_of(pattern, file_path, pattern_is_found) result(iline)
    implicit none
    character(len=*), intent(in) :: pattern
    character(len=*), intent(in) :: file_path
    logical, intent(out) :: pattern_is_found
    integer :: iline, ierr
    character(len=200) :: line
    iline = 0
    open(unit=1, file=file_path, status='old')
    do
        read(1,'(a)', iostat=ierr) line ; if(ierr.ne.0) exit ! Try to read the current line. If it fails, exit the do loop.
        iline = iline+1 ! Increment the line number
        if (index(line, trim(pattern)).ne.0) then
            ! If the substring is found, close the file and return.
            close(1)
            pattern_is_found = .true.
            return
        endif
    enddo
    close(1)
    pattern_is_found = .false. ! If the substring was not found, set the flag to false.
end function get_first_occ_of

subroutine get_header_lines_counts(ref_data_file_path, test_data_file_path, data_line_pattern, num_ref_header_lines, num_test_header_lines, pattern_is_found)
implicit none
character(len=*), intent(in) :: ref_data_file_path  ! eg "/work/hibridon/issues/issue249/hibridon.git/tests/cahe/Job.flx"
character(len=*), intent(in) :: test_data_file_path  ! eg "/work/hibridon/issues/issue249/build/tests/cahe/Job.flx"
character(len=*), intent(in) :: data_line_pattern  ! Pattern for the first data line, eg "R (BOHR) AND OUTGOING FLUXES"
integer, intent(out) :: num_ref_header_lines
integer, intent(out) :: num_test_header_lines
logical, intent(out) :: pattern_is_found
  num_ref_header_lines = get_first_occ_of(data_line_pattern, ref_data_file_path, pattern_is_found)
  if (.not. pattern_is_found) then
    write(Error_Unit, "(a)") "Pattern not found in reference file " // trim(ref_data_file_path)
    return
  end if
  num_test_header_lines = get_first_occ_of(data_line_pattern, test_data_file_path, pattern_is_found)
  if (.not. pattern_is_found) then
    write(Error_Unit, "(a)") "Pattern not found in test file " // trim(test_data_file_path)
    return
  end if
end subroutine






end program comp_tests
