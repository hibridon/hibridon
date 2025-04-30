
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

module m_vector
    implicit none
    private

    type, public :: t_vector
        integer :: vec_size = 0
        real(8), dimension(:), allocatable :: array
    end type

end module m_vector

module m_stringsplitter
contains
! from https://fpm.fortran-lang.org/proc/split.html
subroutine split(input_line,array,delimiters,order,nulls)
    !! given a line of structure " par1 par2 par3 ... parn " store each par(n) into a separate variable in array.
    !!
    !! * by default adjacent delimiters in the input string do not create an empty string in the output array
    !! * no quoting of delimiters is supported
    character(len=*),intent(in)              :: input_line  !! input string to tokenize
    character(len=*),optional,intent(in)     :: delimiters  !! list of delimiter characters
    character(len=*),optional,intent(in)     :: order       !! order of output array sequential|[reverse|right]
    character(len=*),optional,intent(in)     :: nulls       !! return strings composed of delimiters or not ignore|return|ignoreend
    character(len=:),allocatable,intent(out) :: array(:)    !! output array of tokens

    integer                       :: n                      ! max number of strings INPUT_LINE could split into if all delimiter
    integer,allocatable           :: ibegin(:)              ! positions in input string where tokens start
    integer,allocatable           :: iterm(:)               ! positions in input string where tokens end
    character(len=:),allocatable  :: dlim                   ! string containing delimiter characters
    character(len=:),allocatable  :: ordr                   ! string containing order keyword
    character(len=:),allocatable  :: nlls                   ! string containing nulls keyword
    integer                       :: ii,iiii                ! loop parameters used to control print order
    integer                       :: icount                 ! number of tokens found
    integer                       :: ilen                   ! length of input string with trailing spaces trimmed
    integer                       :: i10,i20,i30            ! loop counters
    integer                       :: icol                   ! pointer into input string as it is being parsed
    integer                       :: idlim                  ! number of delimiter characters
    integer                       :: ifound                 ! where next delimiter character is found in remaining input string data
    integer                       :: inotnull               ! count strings not composed of delimiters
    integer                       :: ireturn                ! number of tokens returned
    integer                       :: imax                   ! length of longest token

    ! decide on value for optional DELIMITERS parameter
    if (present(delimiters)) then                                     ! optional delimiter list was present
        if(delimiters.ne.'')then                                       ! if DELIMITERS was specified and not null use it
            dlim=delimiters
        else                                                           ! DELIMITERS was specified on call as empty string
            dlim=' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0) ! use default delimiter when not specified
        endif
    else                                                              ! no delimiter value was specified
        dlim=' '//char(9)//char(10)//char(11)//char(12)//char(13)//char(0)    ! use default delimiter when not specified
    endif
    idlim=len(dlim)                                                   ! dlim a lot of blanks on some machines if dlim is a big string

    if(present(order))then; ordr=(adjustl(order)); else; ordr='sequential'; endif ! decide on value for optional ORDER parameter
    if(present(nulls))then; nlls=(adjustl(nulls)); else; nlls='ignore'    ; endif ! optional parameter

    n=len(input_line)+1                        ! max number of strings INPUT_LINE could split into if all delimiter
    allocate(ibegin(n))                        ! allocate enough space to hold starting location of tokens if string all tokens
    allocate(iterm(n))                         ! allocate enough space to hold ending location of tokens if string all tokens
    ibegin(:)=1
    iterm(:)=1

    ilen=len(input_line)                                           ! ILEN is the column position of the last non-blank character
    icount=0                                                       ! how many tokens found
    inotnull=0                                                     ! how many tokens found not composed of delimiters
    imax=0                                                         ! length of longest token found

    select case (ilen)

    case (0)                                                      ! command was totally blank

    case default                                                   ! there is at least one non-delimiter in INPUT_LINE if get here
        icol=1                                                      ! initialize pointer into input line
        INFINITE: do i30=1,ilen,1                                   ! store into each array element
            ibegin(i30)=icol                                         ! assume start new token on the character
            if(index(dlim(1:idlim),input_line(icol:icol)).eq.0)then  ! if current character is not a delimiter
            iterm(i30)=ilen                                       ! initially assume no more tokens
            do i10=1,idlim                                        ! search for next delimiter
                ifound=index(input_line(ibegin(i30):ilen),dlim(i10:i10))
                IF(ifound.gt.0)then
                    iterm(i30)=min(iterm(i30),ifound+ibegin(i30)-2)
                endif
            enddo
            icol=iterm(i30)+2                                     ! next place to look as found end of this token
            inotnull=inotnull+1                                   ! increment count of number of tokens not composed of delimiters
            else                                                     ! character is a delimiter for a null string
            iterm(i30)=icol-1                                     ! record assumed end of string. Will be less than beginning
            icol=icol+1                                           ! advance pointer into input string
            endif
            imax=max(imax,iterm(i30)-ibegin(i30)+1)
            icount=i30                                               ! increment count of number of tokens found
            if(icol.gt.ilen)then                                     ! no text left
            exit INFINITE
            endif
        enddo INFINITE

    end select

    select case (trim(adjustl(nlls)))
    case ('ignore','','ignoreend')
        ireturn=inotnull
    case default
        ireturn=icount
    end select
    allocate(character(len=imax) :: array(ireturn))                ! allocate the array to return
    !allocate(array(ireturn))                                       ! allocate the array to turn

    select case (trim(adjustl(ordr)))                              ! decide which order to store tokens
    case ('reverse','right') ; ii=ireturn ; iiii=-1                ! last to first
    case default             ; ii=1       ; iiii=1                 ! first to last
    end select

    do i20=1,icount                                                ! fill the array with the tokens that were found
        if(iterm(i20).lt.ibegin(i20))then
            select case (trim(adjustl(nlls)))
            case ('ignore','','ignoreend')
            case default
            array(ii)=' '
            ii=ii+iiii
            end select
        else
            array(ii)=input_line(ibegin(i20):iterm(i20))
            ii=ii+iiii
        endif
    enddo
    deallocate(iterm)
    deallocate(ibegin)
end subroutine split
end module m_stringsplitter

function contains_time_information(line)
    use iso_fortran_env, only: Output_Unit
    use m_stringsplitter, only: split
    ! returns true if the line contains time information (eg parts of a date, elapsed time, etc.)
    character(len=*), intent(in) :: line
    logical :: contains_time_information

    character(len=:), allocatable :: tokens(:)
    integer :: token_index
    call split(line, tokens)
    token_index = 1
    do while(token_index <= size(tokens))
        ! write(Output_Unit,*) 'token : ', tokens(token_index)
        if( tokens(token_index) == 'DATE:' ) then
            ! write(Output_Unit,*) 'DATE : ', line
            contains_time_information = .true.
            return
        endif
        if( tokens(token_index) == 'WRITTEN:' ) then
            ! write(Output_Unit,*) 'WRITTEN : ', line
            contains_time_information = .true.
            return
        endif
        if( tokens(token_index) == 'TIMING:' ) then
            ! write(Output_Unit,*) 'TIMING : ', line
            contains_time_information = .true.
            return
        endif
        token_index = token_index + 1
    enddo

    contains_time_information = .false.
end function contains_time_information

FUNCTION is_numeric(string)
    USE ieee_arithmetic
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: string
    LOGICAL :: is_numeric
    REAL :: x
    INTEGER :: e
    if (string == 'INF') then ! INF happens to be the name of a column
      is_numeric = .false.
    else
      x = 0
      READ(string,*,IOSTAT=e) x
      is_numeric = ((e == 0) .and. (.NOT. ISNAN(X)))
    end if
    END FUNCTION is_numeric

subroutine get_file_numbers(file_name, numbers, num_header_lines)
    use iso_fortran_env, only: Output_Unit
    use m_vector, only: t_vector
    use m_stringsplitter, only: split
    implicit none
    character(len=200) :: file_name
    type(t_vector) :: numbers
    integer :: num_header_lines
    integer ios
    integer num_tokens
    integer num_numbers
    real(8) :: number
    character(len=200) :: line
    character(len=:), allocatable :: tokens(:)
    integer :: file_id = 1
    integer :: line_index = 1
    integer :: number_index = 1
    integer :: token_index = 1
    logical :: is_numeric
    logical :: contains_time_information
    integer, parameter :: PASS_FIRST = 0
    integer, parameter :: PASS_COUNT_NUMBERS = 0
    integer, parameter :: PASS_FILL_VECTOR = 1
    integer, parameter :: PASS_LAST = 1
    
    integer :: pass = PASS_COUNT_NUMBERS

    do pass = PASS_FIRST, PASS_LAST
        if (pass == PASS_COUNT_NUMBERS) then
            num_numbers = 0
            num_tokens = 0
            open(unit=file_id, status='old', file=trim(file_name))
        endif

        if (pass == PASS_FILL_VECTOR) then
            number_index = 1
            rewind(file_id)
        endif

        line_index = 1
        do
            ! write(Output_Unit, *) 'num_numbers = ', num_numbers
            read(file_id, '(A)', iostat=ios) line
            if(ios /= 0) then
                exit
            endif
            ! line = 'rowan atkinson'
            ! write(Output_Unit, '(A)') 'line = ', line        
            ! write(Output_Unit, *) 'num_tokens = ', size(tokens)
            if (( line_index > num_header_lines ) .and. (.not. contains_time_information(line))) then
                call split(line, tokens)
                token_index = 1
                do while(token_index <= size(tokens))
                    if( is_numeric(tokens(token_index)) ) then
                        if (pass == PASS_COUNT_NUMBERS) then
                            num_numbers = num_numbers + 1
                        endif
                        if (pass == PASS_FILL_VECTOR) then
                            read(tokens(token_index),*) number
                            numbers%array(number_index) = number
                            number_index = number_index + 1
                        endif
                    endif
                    token_index = token_index + 1
                    if (pass == PASS_COUNT_NUMBERS) then
                        num_tokens = num_tokens + 1
                    endif
                enddo
            endif
            ! num_numbers = num_numbers + size(tokens)
            line_index = line_index + 1
            ! num_numbers = num_numbers + 1
        enddo
        if (pass == PASS_COUNT_NUMBERS) then
            write(Output_Unit, "(a,a,i0,a,i0)") trim(file_name), ' : total num_tokens = ', num_tokens, &
            ', total num_numbers = ', num_numbers
            allocate(numbers%array(num_numbers))
        endif
        if (pass == PASS_FILL_VECTOR) then
            close(file_id)
        endif
    enddo

end

module m_diff
    implicit none
    contains
    function vectors_differ(va, vb, tolerance, min_significant_value)
        use iso_fortran_env, only: Error_Unit, Output_Unit
        implicit none
        logical :: vectors_differ
        real(8), intent(in) :: va(:)
        real(8), intent(in) :: vb(:)
        real(8), intent(in) :: tolerance
        real(8), intent(in) :: min_significant_value
        integer :: number_index

        vectors_differ = .FALSE.

        if ( size(va) /= size(vb) ) then
            write(Error_Unit, "(a,i0,a,i0)") "different number of numbers found in results files: ", &
                & size(va), ' <> ', size(vb)
            vectors_differ = .TRUE.
            return
        else
            number_index = 1
            do while(number_index <= size(va))
                if(abs(va(number_index)) > min_significant_value ) then
                    if( abs(va(number_index) - vb(number_index))/max(abs(va(number_index)),1d-300) > tolerance ) then
                        write(Error_Unit, "(a,i0,a,e12.5,a,e12.5,a,f3.1,a)") "at number_index ", number_index, ": ",&
                        va(number_index), " and ", vb(number_index), " differ by more than ", tolerance*100.0d0,&
                        " percent"
                        vectors_differ = .TRUE.
                        !return
                    endif
                endif
                number_index = number_index + 1
            enddo
        endif
        return
        end function vectors_differ

    function result_files_differ(results1_file_name, results2_file_name, num_header_lines, tolerance, min_significant_value)
        use m_vector, only: t_vector
        ! use m_diff, only: vectors_differ
        implicit none
        character(len=200), intent(in) :: results1_file_name
        character(len=200), intent(in) :: results2_file_name
        integer, intent(in) :: num_header_lines(2)
        real(8), intent(in) :: tolerance
        real(8), intent(in) :: min_significant_value
        logical :: result_files_differ
    
        type(t_vector) :: ref_numbers
        type(t_vector) :: test_numbers
        ! integer :: number_index
        ! real(8) :: number
    
        call get_file_numbers(results1_file_name, ref_numbers, num_header_lines(1))
        ! number_index = 1
        ! do while(number_index <= size(ref_numbers%array))
        !     number = ref_numbers%array(number_index)
        !     write(Output_Unit, *) "number(", number_index, ") = ", number
        !     number_index = number_index + 1
        ! enddo
    
        call get_file_numbers(results2_file_name, test_numbers, num_header_lines(2))
    
        result_files_differ = vectors_differ(ref_numbers%array, test_numbers%array, tolerance, min_significant_value)
    
        end function result_files_differ
        
endmodule m_diff

