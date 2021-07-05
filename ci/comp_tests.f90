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
end subroutine split
end module m_stringsplitter

FUNCTION is_numeric(string)
    USE ieee_arithmetic
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: string
    LOGICAL :: is_numeric
    REAL :: x
    INTEGER :: e
    x = 0
    READ(string,*,IOSTAT=e) x
    is_numeric = ((e == 0) .and. (.NOT. ISNAN(X)))
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

    num_numbers = 0
    num_tokens = 0
    line_index = 1
    open(unit=1, status='old', file=trim(file_name))
    do
        ! write(Output_Unit, *) 'num_numbers = ', num_numbers
        read(file_id, '(A)', iostat=ios) line
        if(ios /= 0) then
            exit
        endif
        call split(line, tokens)
        ! line = 'rowan atkinson'
        ! write(Output_Unit, '(A)') 'line = ', line        
        ! write(Output_Unit, *) 'num_tokens = ', size(tokens)
        if ( line_index > num_header_lines ) then
            token_index = 1
            do while(token_index <= size(tokens))
                if( is_numeric(tokens(token_index)) ) then
                    num_numbers = num_numbers + 1
                endif
                token_index = token_index + 1
                num_tokens = num_tokens + 1
            enddo
        endif
        ! num_numbers = num_numbers + size(tokens)
        line_index = line_index + 1
        ! num_numbers = num_numbers + 1
    enddo
    write(Output_Unit, *) 'total num_tokens = ', num_tokens
    write(Output_Unit, *) 'total num_numbers = ', num_numbers

    allocate(numbers%array(num_numbers))

    number_index = 1
    rewind(file_id)
    line_index = 1
    do
        read(file_id, '(A)', iostat=ios) line
        if(ios /= 0) then
            exit
        endif
        call split(line, tokens)
        ! line = 'rowan atkinson'
        ! write(Output_Unit, '(A)') 'line = ', line
        ! write(Output_Unit, *) 'num_tokens = ', size(tokens)
        if ( line_index > num_header_lines ) then
            token_index = 1
            do while(token_index <= size(tokens))
                if( is_numeric(tokens(token_index)) ) then
                    read(tokens(token_index),*) number
                    numbers%array(number_index) = number
                    number_index = number_index + 1
                endif
                token_index = token_index + 1
            enddo
        endif
        line_index = line_index + 1
        ! num_numbers = num_numbers + 1
    enddo

    close(file_id)
    write(Output_Unit, *) 'at end, num_numbers = ', num_numbers
end

module m_diff
    implicit none
    contains
    function vectors_differ(va, vb, tolerance)
        use iso_fortran_env, only: Error_Unit, Output_Unit
        implicit none
        logical :: vectors_differ
        real(8), intent(in) :: va(:)
        real(8), intent(in) :: vb(:)
        real(8) :: tolerance
        integer :: number_index

        if ( size(va) /= size(vb) ) then
            write(Error_Unit, *) "different number of numbers found in results files : ", &
                & size(va), ' <> ', size(vb)
            vectors_differ = .TRUE.
            return
        else
            number_index = 1
            do while(number_index <= size(va))
                if( abs(va(number_index) - vb(number_index))/max(va(number_index),1d-300) > tolerance ) then
                    write(Error_Unit, *) "at number_index ", number_index, " : ", va(number_index), &
                      & " and ", vb(number_index), " differ for more than ", tolerance*100.0d0, " percent"
                    vectors_differ = .TRUE.
                    return
                endif
                number_index = number_index + 1
            enddo
            vectors_differ = .TRUE.
            return
        endif
        end function vectors_differ

    function result_files_differ(results1_file_name, results2_file_name, num_header_lines, tolerance)
        use m_vector, only: t_vector
        ! use m_diff, only: vectors_differ
        implicit none
        character(len=200), intent(in) :: results1_file_name
        character(len=200), intent(in) :: results2_file_name
        integer, intent(in) :: num_header_lines
        real(8) :: tolerance
        logical :: result_files_differ
    
        type(t_vector) :: ref_numbers
        type(t_vector) :: test_numbers
        ! integer :: number_index
        ! real(8) :: number
    
        call get_file_numbers(results1_file_name, ref_numbers, num_header_lines)
        ! number_index = 1
        ! do while(number_index <= size(ref_numbers%array))
        !     number = ref_numbers%array(number_index)
        !     write(Output_Unit, *) "number(", number_index, ") = ", number
        !     number_index = number_index + 1
        ! enddo
    
        call get_file_numbers(results2_file_name, test_numbers, num_header_lines)
    
        result_files_differ = vectors_differ(ref_numbers%array, test_numbers%array, tolerance=tolerance)
    
        end function result_files_differ
        
endmodule m_diff


program comp_tests
    use iso_fortran_env, only: Error_Unit, Output_Unit
    use m_diff, only: result_files_differ

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
        write(Error_Unit, *) "error: unhandled file extension : ", ext
        stop 1
    endif

    if(ext=="pcs") then
        if(result_files_differ(ref, test, num_header_lines=3, tolerance=0.01d0)) then
            stop 1
        else
            stop 0
        endif
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