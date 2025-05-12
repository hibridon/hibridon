#include "assert.h"
#include "unused.h"
#include "command.inc.F90"

module mod_hicommands
use mod_command, only: command_type, command_mgr_type
implicit none
  class(command_mgr_type), pointer :: command_mgr  ! singleton

  ! check
  type, extends(command_type) :: check_command_type
  contains
    procedure :: execute => check_execute
  end type check_command_type

  ! eadiab
  type, extends(command_type) :: eadiab_command_type
  contains
    procedure :: execute => eadiab_execute
  end type eadiab_command_type

  ! flux
  type, extends(command_type) :: flux_command_type
  contains
    procedure :: execute => flux_execute
  end type flux_command_type

  ! hypxsc
  type, extends(command_type) :: hypxsc_command_type
  contains
    procedure :: execute => hypxsc_execute
  end type hypxsc_command_type

  ! indout
  type, extends(command_type) :: indout_command_type
  contains
    procedure :: execute => indout_execute
  end type indout_command_type

  ! input
  type, extends(command_type) :: input_command_type
  contains
    procedure :: execute => input_execute
  end type input_command_type

  ! intcrs
  type, extends(command_type) :: intcrs_command_type
  contains
    procedure :: execute => intcrs_execute
  end type intcrs_command_type

  ! j1j2
  type, extends(command_type) :: j1j2_command_type
  contains
    procedure :: execute => j1j2_execute
  end type j1j2_command_type

  ! job
  type, extends(command_type) :: job_command_type
  contains
    procedure :: execute => job_execute
  end type job_command_type

  ! label
  type, extends(command_type) :: label_command_type
  contains
    procedure :: execute => label_execute
  end type label_command_type

  ! output
  type, extends(command_type) :: output_command_type
  contains
    procedure :: execute => output_execute
  end type output_command_type

  ! partc
  type, extends(command_type) :: partc_command_type
  contains
    procedure :: execute => partc_execute
  end type partc_command_type

  ! prsbr
  type, extends(command_type) :: prsbr_command_type
  contains
    procedure :: execute => prsbr_execute
  end type prsbr_command_type

  ! read
  type, extends(command_type) :: read_command_type
  contains
    procedure :: execute => read_execute
  end type read_command_type

  ! run
  type, extends(command_type) :: run_command_type
  contains
    procedure :: execute => run_execute
  end type run_command_type

  ! show
  type, extends(command_type) :: show_command_type
  contains
    procedure :: execute => show_execute
  end type show_command_type

  ! stmix
  type, extends(command_type) :: stmix_command_type
  contains
    procedure :: execute => stmix_execute
  end type stmix_command_type

  ! sysconf
  type, extends(command_type) :: sysconf_command_type
  contains
    procedure :: execute => sysconf_execute
  end type sysconf_command_type

  ! trnprt
  type, extends(command_type) :: trnprt_command_type
  contains
    procedure :: execute => trnprt_execute
  end type trnprt_command_type

  ! turn
  type, extends(command_type) :: turn_command_type
  contains
    procedure :: execute => turn_execute
  end type turn_command_type

contains

  subroutine update_nu_params()
    !  numin and numax should be 0 if cc calculation, if not, then set them
    !  equal to zero
    use mod_par, only: lpar, ipar
    use ipar_enum
    use lpar_enum
    use mod_selb, only: ibasty

    if (.not. lpar(LPAR_CSFLAG)) then
      lpar(LPAR_NUCROS)=.false.
      if (ipar(IPAR_NUMAX) .ne. 0) then
        write (6, 105)
        105     format ('  CC calculation, numax set to zero')
        ipar(IPAR_NUMAX) = 0
      end if
      ! NB this is disabled currently for 2P atom + homonuclear
      if (ipar(IPAR_NUMIN) .ne. 0.and.ibasty.ne.12) then
        write (6, 106)
        106     format ('  CC calculation, numin set to zero')
        ipar(IPAR_NUMIN) = 0
      end if
    end if
  end subroutine update_nu_params

  subroutine print_main_params(out_unit)
    use mod_hibasis, only: basknd
    use mod_cosys, only: scod
    use mod_cosysl, only: islcod, lspar
    use mod_cosysi, only: isicod, ispar
    use mod_cosysr, only: isrcod, rspar
    use mod_hinput_state, only: lindx
    use mod_si_params, only: lcode
    use mod_par, only: lpar
    use mod_parbas, only: lammin, lammax, mproj
    use lpar_enum
    use mod_selb, only: ibasty
    use mod_par, only: fcod
    use mod_two, only: nj1j2
    implicit none
    integer, intent(in) :: out_unit

    integer :: j
    integer :: length
    integer :: numj

    if (ibasty .lt. 99) then
      length = index(basknd(ibasty),' ') - 1
      if (length .eq. -1) length=9
#if defined(HIB_UNIX) || defined(HIB_CRAY) || defined(HIB_MAC)
      write(out_unit,710) basknd(ibasty)(1:length)//' system parameters:', &
                 (scod(j),ispar(j),j = 1,isicod)
#endif
    else
      write(out_unit,710) 'user defined system parameters:', &
                   (scod(j),ispar(j),j = 1,isicod)
    endif
    if(isrcod.gt.0) &
      write(out_unit,720) (scod(isicod+j),rspar(j),j = 1,isrcod)
    if(islcod.gt.0) &
      write(out_unit,735) (scod(isicod+isrcod+j),lspar(j),j = 1,islcod)
    if (.not. lpar(LPAR_TWOMOL) ) then
      write(out_unit,736) 'LAMMIN: ',(lammin(j),j=1,ispar(1))
      write(out_unit,736) 'LAMMAX: ',(lammax(j),j=1,ispar(1))
      write(out_unit,736) 'MPROJ:  ',(mproj(j),j=1,ispar(1))
    else if (lpar(LPAR_TWOMOL)) then
      if (allocated(nj1j2)) then
        numj = size(nj1j2)
      else
        numj = 0
      end if
      write (out_unit, 738)'J1/J2: ',(nj1j2(j)/10,mod(nj1j2(j),10), &
                        j=1,numj)
    end if
    write(out_unit,730) 'Flags:',(fcod(j),lpar(lindx(j)),j = 1,lcode)
    710 format(5x,'*** ',(a)/(4(1x,a7,'=',i4,7x)))
    720 format(4(1x,a7,'=',1pg11.4))
    730 format(5x,'*** ',(a)/(6(1x,a6,'=',l2,3x)))
    735 format(3(1x,a7,'=',l2,9x))
    736   format(1x,(a),10i4,/,9x,10i4)
    738   format (1x,(a),1x,20(2i1,'  ') )

  end subroutine print_main_params

  function get_max_energy(success) result(max_energy)
    use mod_par, only: ipar
    use mod_coener, only: energ
    use ipar_enum, only: IPAR_NERG

    logical, intent(out) :: success
    real(8) :: max_energy
    real(8) :: e
    integer :: i

    success = .true.
    e = 0
    do i = 1, ipar(IPAR_NERG)
      e = max(e,energ(i))
    end do
    max_energy = e
    if(e .eq. 0) then
      write(6,1610)
      1610   format(' Total energy has not been given a value !')
      success = .false.
    end if
  end function get_max_energy

  subroutine parse_int_array(statement_parser, array_name, int_array, success)
    use mod_statement_parser, only: statement_parser_type
    use mod_command, only: k_post_action_interpret_next_statement
    use mod_hiutil, only: get_token, assignment_parse

    class(statement_parser_type), intent(inout) :: statement_parser
    character(len=*), intent(in) :: array_name
    integer, allocatable, intent(out) :: int_array(:)
    logical, intent(out) :: success

    character(len=:), allocatable :: argument
    integer :: num_parsed_elements, j
    real(8) :: arg_val
    character*8 empty_var_list(0)
    integer :: num_elements

    success = .true.
    num_parsed_elements = 0
    argument = statement_parser%get_token(equal_is_delimiter=.false.)  ! disable-warnings:maybe-uninitialized (argument)
    call assignment_parse(argument, empty_var_list, j, arg_val)
    num_elements = arg_val
    if(num_elements .lt. 0) then
      write (6, 485), array_name
    485   format(' ** YOU MUST SPECIFY THE SIZE OF ARRAY ', (a))
    end if

    if(num_elements .ne. 0) then
      if ( allocated(int_array) ) deallocate(int_array)
      allocate(int_array(num_elements))
      do while(.not. statement_parser%statement_end_reached())
        if(statement_parser%prev_char_is(';')) exit
        argument = statement_parser%get_token(equal_is_delimiter=.false.)
        call assignment_parse(argument, empty_var_list, j, arg_val)
        num_parsed_elements = num_parsed_elements + 1
        if (num_parsed_elements > num_elements) then
          write (6, *) 'ERROR: the number of parsed elements exceeds the size of the array (', num_elements, ')'
          success = .false.
          return
        end if
        int_array(num_parsed_elements) = arg_val
      end do
      if(num_parsed_elements /= num_elements) then
        write (6, *) 'ERROR: only ', num_parsed_elements,' integers read while ', num_elements, ' were expected'
        success = .false.
        return
      endif
    end if
  end subroutine parse_int_array

  ! check if inconsistencies in input parameters
  subroutine check_execute(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    use mod_command, only: k_post_action_read_new_line
    use mod_hiiolib1, only: genchk

    class(check_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action
    UNUSED_DUMMY(this)
    UNUSED_DUMMY(statement_parser)
    call genchk()
    post_action = k_post_action_read_new_line

  end subroutine check_execute

  subroutine eadiab_execute(this, statement_parser, post_action)
    ! adiabatic energy calculation
    ! eadiab,jobfile,[nch1],[nch2]
    ! with eadiab,jobfile:
    !  nchmin = 1
    !  nchmax = 10
    ! with eadiab,jobfile,nch1:
    !  nchmin = 1
    !  nchmax = nch1
    ! with eadiab,jobfile,nch1,nch2:
    !  nchmin = nch1
    !  nchmax = nch2

    use mod_statement_parser, only: statement_parser_type
    use mod_command, only: k_post_action_read_new_line
    use mod_hibrid4, only: eadiab1
    use mod_file, only: jobnam
    use mod_hiutil, only: lower, upper

    class(eadiab_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    integer :: nchmin, nchmax
    character(len=:), allocatable :: fnam1
    character(len=:), allocatable :: arg

    UNUSED_DUMMY(this)

    fnam1 = statement_parser%get_token(equal_is_delimiter=.false.)  ! disable-warnings:maybe-uninitialized (fnam1)
    call lower(fnam1)
    call upper(fnam1(1:1))
    if(fnam1 .eq. ' ') fnam1 = jobnam
    arg = statement_parser%get_token(equal_is_delimiter=.false.)  ! disable-warnings:maybe-uninitialized (arg)
    if (arg .eq. ' ') then
      nchmin = 1
      nchmax = 10
    else
      read (arg, *, err=2860, end=2860) nchmax
      arg = statement_parser%get_token(equal_is_delimiter=.false.)
      if (arg .eq. ' ') then
        nchmin = 1
      else
        nchmin = nchmax
        read (arg, *, err=2860, end=2860) nchmax
      end if
    end if
    call eadiab1(fnam1, nchmin, nchmax)
2860 write (6, *) 'Parameters to EADIAB cannot be recognized'

    post_action = k_post_action_read_new_line


  end subroutine eadiab_execute

  subroutine flux_execute(this, statement_parser, post_action)
    !  flux calculation,jobfile,mchannel,iflux,thresh,iprint

    use mod_statement_parser, only: statement_parser_type
    use mod_command, only: k_post_action_read_new_line
    use mod_hibrid4, only: psi
    use mod_file, only: jobnam
    use mod_hiutil, only: lower, upper, assignment_parse

    class(flux_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    character(len=:), allocatable :: fnam1
    character(len=:), allocatable :: arg
    real(8) :: a(7)
    character*8 empty_var_list(0)
    integer :: i, j

    UNUSED_DUMMY(this)

    fnam1 = statement_parser%get_token(equal_is_delimiter=.false.)  ! disable-warnings:maybe-uninitialized (argument)
    call lower(fnam1)
    call upper(fnam1(1:1))
    if(fnam1 .eq. ' ') fnam1 = jobnam
    do i = 1, 7
      a(i) = 0
      if(.not. statement_parser%statement_end_reached()) then
        arg = statement_parser%get_token(equal_is_delimiter=.false.)
        call assignment_parse(arg, empty_var_list,j, a(i))
      end if
    end do
    call psi(fnam1,a)

    post_action = k_post_action_read_new_line
  end subroutine flux_execute

  subroutine hypxsc_execute(this, statement_parser, post_action)
    !  hyperfine xcs routine (originally written by j. klos,
    !  rewritten by p.j. dagdigian
    !  hypxsc,jobfile, ienerg ,nucspin, j1, j2
    use mod_statement_parser, only: statement_parser_type
    use mod_hiutil, only: get_token, lower, upper, assignment_parse
    use mod_file, only: jobnam
    use mod_command, only: k_post_action_read_new_line
    use mod_hypxsc, only: hypxsc

    class(hypxsc_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    character(len=40) :: fnam1  ! job file name
    real(8) :: a(4)  ! 4 real arguments:
    ! 1: ienerg
    ! 2: nucspin
    ! 3: j1
    ! 4: j2

    character*8 empty_var_list(0)
    character(len=40) :: code
    integer :: i, j

    UNUSED_DUMMY(this)
    fnam1 = statement_parser%get_token(equal_is_delimiter=.false.)
    if(fnam1 .eq. ' ') fnam1 = jobnam
    call lower(fnam1)
    call upper(fnam1(1:1))
    do 2013 i = 1,4
       a(i) = 0.d0
       if(statement_parser%statement_end_reached()) goto 2013
       code = statement_parser%get_token(equal_is_delimiter=.false.)
       call assignment_parse(code,empty_var_list,j,a(i))
    2013   continue
    call hypxsc(fnam1,a)
    post_action = k_post_action_read_new_line
  end subroutine hypxsc_execute




  ! indout values
  ! specify indout values in the form
  ! indout,niout,indout(1),...,indout(niout)
  ! terminate the string with a semicolon if other parameters will follow
  ! on the same card, e.g. indout,2,1,-1;energ=e1,e2,e3;jtot1=0,jtot2=2....
  subroutine indout_execute(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    use mod_command, only: k_post_action_interpret_next_statement
    use mod_hiutil, only: get_token, assignment_parse
    use mod_coiout, only: niout, indout

    class(indout_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    character(len=:), allocatable :: argument
    integer :: i, j
    real(8) :: arg_val
    character*8 empty_var_list(0)
    logical :: parse_succeeded

    UNUSED_DUMMY(this)
    
    call parse_int_array(statement_parser, 'indout', indout, parse_succeeded)
    if (parse_succeeded) then
      niout = size(indout)
    end if

    post_action = k_post_action_interpret_next_statement
  end subroutine indout_execute

  ! input, output, label and job file names
  ! input=infile, output=outfile, job=jobfile
  ! input, output, and label are now lower case:  mha 6.6.91
  subroutine input_execute(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    use mod_command, only: k_post_action_exit_hibridon, k_post_action_read_new_line
    use mod_hinput_state, only: batch
    use mod_file, only: input
    use mod_hiutil, only: get_token, lower, upper

    class(input_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    character(len=40) :: argument
    integer :: len
    class(command_type), pointer :: command
    logical :: existf
    UNUSED_DUMMY(this)
    argument = statement_parser%get_token(equal_is_delimiter=.false.)
    input = argument
    call lower(input)
    call upper(input(1:1))
    inquire (file=input, exist=existf)
    if (.not. existf) then
      len = index (input,' ')
      write (6, 901) input(1:len)
  901     format (' *** INPUT FILE ',(a),' DOES NOT EXIST')
      if(batch) then
        post_action = k_post_action_exit_hibridon
      else
        post_action = k_post_action_read_new_line
      end if
    end if
    command => command_mgr%get_command('READ')
    call command%execute(statement_parser, post_action)
  end subroutine input_execute

  subroutine intcrs_execute(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    use mod_hibrid5, only :intcrs
    !.....integral cross sections
    !  intcrs,jobfile,in1,in2,ienerg,maxjtot
    use mod_hiutil, only: assignment_parse
    use mod_command, only: k_post_action_read_new_line
    use mod_file, only: jobnam
    use mod_hiutil, only: get_token, lower, upper

    implicit none
    class(intcrs_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action
    character(len=40) :: code


    character(len=40) :: fnam1

    integer :: i, j
    integer, parameter :: k_num_args = 4
    real(8) :: a(k_num_args)  ! in1,in2,ienerg,maxjtot
    character*8 empty_var_list(0)

    UNUSED_DUMMY(this)
    fnam1 = statement_parser%get_token(equal_is_delimiter=.false.)
    if(fnam1 .eq. ' ') fnam1 = jobnam
    call lower(fnam1)
    call upper(fnam1(1:1))
    do i = 1, k_num_args
      a(i) = 0
      if(statement_parser%statement_end_reached()) cycle
      code = statement_parser%get_token(equal_is_delimiter=.false.)
      call assignment_parse(code,empty_var_list,j,a(i))
    end do
    call intcrs(fnam1,a)
    post_action = k_post_action_read_new_line
  end subroutine intcrs_execute

  ! j1j2 values
  ! specify j1j2 values in the form
  ! j1j2,numj,j1j2(1),...,j1j2(numj)
  ! terminate the string with a semicolon if other parameters will follow
  ! on the same card, e.g. j1j2,2,00,10;energ=e1,e2,e3;jtot1=0,jtot2=2....
  subroutine j1j2_execute(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    use mod_command, only: k_post_action_interpret_next_statement
    use mod_hiutil, only: get_token, assignment_parse
    use mod_two, only: nj1j2
    use mod_par, only: lpar
    use lpar_enum, only: LPAR_TWOMOL

    class(j1j2_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    character(len=:), allocatable :: argument
    integer :: i, j
    real(8) :: arg_val
    character*8 empty_var_list(0)
    logical :: parse_succeeded

    UNUSED_DUMMY(this)

    if (.not. lpar(LPAR_TWOMOL)) then
      write (6, 465)
    465   format(' ** NUMJ CAN ONLY BE DEFINED IF TWOMOL = .TRUE.')
      post_action = k_post_action_interpret_next_statement
      return
    endif
    
    call parse_int_array(statement_parser, 'nj1j2', nj1j2, parse_succeeded)
    if (.not. parse_succeeded) then
      post_action = k_post_action_interpret_next_statement
      return
    end if

    post_action = k_post_action_interpret_next_statement
  end subroutine j1j2_execute

  ! input, output, label and job file names
  ! input=infile, output=outfile, job=jobfile
  ! input, output, and label are now lower case:  mha 6.6.91
  subroutine job_execute(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    use mod_command, only: k_post_action_interpret_next_statement
    use mod_file, only: jobnam
    use mod_hiutil, only: get_token, lower, upper

    class(job_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    character(len=40) :: code
    UNUSED_DUMMY(this)
    code = statement_parser%get_token(equal_is_delimiter=.false.)
    jobnam = code
    call lower(jobnam)
    call upper(jobnam(1:1))
    post_action = k_post_action_interpret_next_statement
  end subroutine job_execute

  ! input, output, label and job file names
  ! input=infile, output=outfile, job=jobfile
  ! input, output, and label are now lower case:  mha 6.6.91
  subroutine label_execute(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    use mod_command, only: k_post_action_interpret_next_statement
    use mod_parpot, only: label => pot_label
    use mod_hiutil, only: lenstr

    class(label_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    integer :: l

    integer :: low, lend

    UNUSED_DUMMY(this)

    low = index (statement_parser%current_line, '=') + 1
    lend = index (statement_parser%current_line, ';') - 1
    if(lend.lt.0) then
      lend=lenstr(statement_parser%current_line)
      l=0
    else
      l=lend+2
    end if
    label = statement_parser%current_line(low:lend)

    statement_parser%current_pos = l
    post_action = k_post_action_interpret_next_statement
  end subroutine label_execute

  ! input, output, label and job file names
  ! input=infile, output=outfile, job=jobfile
  ! input, output, and label are now lower case:  mha 6.6.91
  subroutine output_execute(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    use mod_command, only: k_post_action_interpret_next_statement
    use funit
    use mod_file, only: output
    use mod_hiutil, only: get_token, lower, upper
    use mod_hiiolib1, only: openf

    class(output_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    character(len=40) :: argument

    logical :: openfl

    UNUSED_DUMMY(this)
    argument = statement_parser%get_token(equal_is_delimiter=.false.)
    output = argument
    call lower(output)
    call upper(output(1:1))
    inquire(FUNIT_OUT, opened=openfl)
    if(openfl) then
      endfile(FUNIT_OUT)
      close(FUNIT_OUT)
    end if
    call openf(FUNIT_OUT, output, 'sf', 0)
    post_action = k_post_action_interpret_next_statement
  end subroutine output_execute

  subroutine partc_execute(this, statement_parser, post_action)
    ! print selected partial cross sections from pcs file
    use mod_statement_parser, only: statement_parser_type
    use mod_command, only: k_post_action_read_new_line
    use mod_file, only: jobnam
    use mod_hiutil, only: get_token, lower, upper, assignment_parse
    use mod_hibrid5, only : readpc
    use mod_codim, only: nmax => mmax
    use mod_coamat, only: scmat => toto ! scmat(1)
    class(partc_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action
    character*8 empty_var_list(0)
    character(len=:), allocatable :: arg
    character(len=:), allocatable :: fnam1

    integer :: i, j
    integer, parameter :: k_num_args = 8
    real(8) :: a(k_num_args)

    UNUSED_DUMMY(this)

    fnam1 = statement_parser%get_token(equal_is_delimiter=.false.)  ! disable-warnings:maybe-uninitialized (fnam1)
    if(fnam1 .eq. ' ') fnam1 = jobnam
    call lower(fnam1)
    call upper(fnam1(1:1))
    do i = 1, k_num_args
      a(i) = 0
      if(.not. statement_parser%statement_end_reached()) then
        arg = statement_parser%get_token(equal_is_delimiter=.false.)
        call assignment_parse(arg, empty_var_list, j, a(i))
      end if
    end do
    call readpc(fnam1, a, scmat, nmax)
    post_action = k_post_action_read_new_line
  end subroutine partc_execute

  subroutine prsbr_execute(this, statement_parser, post_action)
    ! pressure broadening cross sections - added by p. dagdigian
    use mod_statement_parser, only: statement_parser_type
    use mod_hiutil, only: assignment_parse
    use mod_command, only: k_post_action_read_new_line
    use mod_file, only: jobnam
    use mod_hiutil, only: get_token, lower, upper
    use mod_hiprsbr, only: prsbr
    class(prsbr_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action
    character*8 empty_var_list(0)
    character*40 :: code
    character(len=40) :: fnam1
    character(len=40) :: fnam2

    integer :: i, j
    integer, parameter :: k_num_args = 12
    real(8) :: a(k_num_args)

    UNUSED_DUMMY(this)

    fnam1 = statement_parser%get_token(equal_is_delimiter=.false.)
    if(fnam1 .eq. ' ') fnam1 = jobnam
    call lower(fnam1)
    call upper(fnam1(1:1))
    ! get iener for 1st smt file
    a(1) = 0.d0
    if(statement_parser%statement_end_reached()) goto 3205
    code = statement_parser%get_token(equal_is_delimiter=.false.)
    call assignment_parse(code,empty_var_list,j,a(1))
    3205 fnam2 = statement_parser%get_token(equal_is_delimiter=.false.)
    if(fnam2 .eq. ' ') fnam2 = jobnam
    call lower(fnam2)
    call upper(fnam2(1:1))
    ! get iener for 2nd smt file
    a(2) = 0.d0
    if(statement_parser%statement_end_reached()) goto 3210
    code = statement_parser%get_token(equal_is_delimiter=.false.)
    call assignment_parse(code,empty_var_list,j,a(2))
    ! get k, j1, in1, j2, in2, diag, j1p, in1p, j2p, in2p
    3210 do 3220 i = 3, 12
      a(i) = 0.d0
      if(statement_parser%statement_end_reached()) goto 3220
      code = statement_parser%get_token(equal_is_delimiter=.false.)
      call assignment_parse(code,empty_var_list,j,a(i))
    3220 continue
    call prsbr(fnam1,fnam2,a)
    post_action = k_post_action_read_new_line
  end subroutine prsbr_execute

  subroutine read_execute(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    use mod_hinput_state, only: batch
    use mod_par, only: lpar
    use mod_si_params, only: icode, set_param_names
    use mod_hinput_state, only: irpot, irinp
    use mod_command, only: k_post_action_interpret_next_statement
    use lpar_enum
    use mod_par, only: pcod
    use mod_hiiolib1, only: gendat
    use mod_hisystem, only: sysdat

    class(read_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    integer :: ione
    UNUSED_DUMMY(this)
    UNUSED_DUMMY(statement_parser)

    ! read
    call set_param_names(lpar(LPAR_BOUNDC),pcod,icode)
    call gendat
    ione = 1
    call sysdat(irpot, lpar(LPAR_READPT),ione)
    irinp = 1
    if (batch) lpar(LPAR_BATCH) = .true.
    post_action = k_post_action_interpret_next_statement
  end subroutine read_execute

  subroutine run_execute(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    use mod_par, only: ipar, lpar, rpar
    use mod_si_params, only: iicode, ircode, icode, set_param_names
    use ipar_enum
    use lpar_enum
    use mod_coiout, only: niout, indout
    use mod_cosout, only: nnout, jout
    use mod_coener, only: energ
    use mod_hibrid2, only: enord
    use mod_hinput_state, only: irpot, irinp
    use mod_command, only: k_post_action_write_cr_and_exit, k_post_action_exit_hinput, k_post_action_exit_hibridon, k_post_action_read_new_line
    use rpar_enum
    use mod_parpot, only: potnam=>pot_name, label=>pot_label
    use mod_hiiolib1, only: genchk
    !
    ! label:execute_run
    !
    ! start execution,run
    use mod_command, only: k_post_action_read_new_line
    use mod_file, only: output
    use mod_par, only: pcod
    use mod_sav, only: ixpar, rxpar
    use mod_hiutil, only: get_token, lower, upper
    use mod_hiiolib1, only: openf

    implicit none
    class(run_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    integer :: i, j
    integer :: nerg

    UNUSED_DUMMY(this)
    UNUSED_DUMMY(statement_parser)

    call update_nu_params()
    call set_param_names(lpar(LPAR_BOUNDC),pcod,icode)
    if (lpar(LPAR_CSFLAG).and.ipar(IPAR_NUD).ne.1) lpar(LPAR_NUCROS)=.true.
    nerg=ipar(IPAR_NERG)
    ! check to see if flags are ok if wavefunction desired or
    ! photodissociation calculation
    call genchk
    call enord(energ,nerg)
    do i = 1,ircode
      rxpar(i) = rpar(i)
    end do
    do i = 1,iicode
      ixpar(i) = ipar(i)
    end do
    if(irinp.eq.0) then
      write(6,505)
      505 format (/,' ** SAVE DEFAULT VARIABLES OR SPECIFY INPUT', &
              ' FILE WITH INP = filename')
      if(lpar(LPAR_BATCH)) then
        post_action = k_post_action_exit_hibridon
        return
      else
        post_action = k_post_action_read_new_line
        return
      end if
    end if
    if(rpar(RPAR_SCAT_TOLAI).eq.0) then  ! graffy: todo : shouldn't it be RPAR_XMU instead of RPAR_SCAT_TOLAI here ?
      write(6,507)
      507   format(/,' ** SPECIFY COLLISION REDUCED MASS WITH XMU = mass')
      post_action = k_post_action_read_new_line
      return
    end if
    if(irpot.ne.0.or..not.lpar(LPAR_READPT)) then
    ! open output file
    ! first make sure it is lower case
      call lower(output)
      call upper(output(1:1))
      call openf(9,output,'sf',0)
    !     write input data to file outpt
      write (9, 508) label
      508   format (1x,a)
      write (9, 508) potnam
      write (9, 240)
      240   format(1h ,30('='))
      write(9,710) 'Parameters:',(pcod(j),ipar(j),j = 1,iicode)
      write(9,720) (pcod(iicode+j),rpar(j),j = 1,ircode)
      call print_main_params(out_unit=9)
      710 format(5x,'*** ',(a)/(4(1x,a7,'=',i4,7x)))
      720 format(4(1x,a7,'=',1pg11.4))

      736   format(1x,(a),10i4,/,9x,10i4)
      737   format(1x,(a),i2,(a),20(i3,1x))
      739   format (1x,(a),i2,(a),20(2i1,'  ') )
      call enord(energ,ipar(IPAR_NERG))
      if(ipar(IPAR_NERG).gt.0) write(9,740) (energ(j),j = 1,ipar(IPAR_NERG))
      740 format(1x,'** Energies:',(t15,5f15.6))
      if(nnout.ne.0) then
        if (.not.lpar(LPAR_TWOMOL) ) then
          write (9,737) 'NOUT: ',nnout, &
               '; JOUT:',(jout(j), j=1,iabs(nnout) )
        else
          write (9,739) 'NOUT: ',nnout, &
             '; J1/J2-OUT: ', &
            (jout(j)/10, mod (jout(j),10), j = 1,iabs(nnout))
        end if
      end if
      if(niout.ne.0) write (9,736) 'INDOUT: ',(indout(j), j=1,niout)
      write (9, 240)
      post_action = k_post_action_exit_hinput
      return
    else
      write(6,510)
    510   format(' Potential not yet defined!')
      post_action = k_post_action_read_new_line
      return
    end if

    post_action = k_post_action_write_cr_and_exit
  end subroutine run_execute

  subroutine show_execute(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    use mod_par, only: ipar, lpar, rpar
    use mod_hinput_state, only: batch
    use mod_si_params, only: iicode, ircode, icode, set_param_names
    use ipar_enum
    use lpar_enum
    use mod_coiout, only: niout, indout
    use mod_conlam, only: nlammx
    use mod_codim, only: nmax => mmax
    use mod_cosout, only: nnout, jout
    use mod_coener, only: energ
    use mod_command, only: k_post_action_interpret_next_statement, k_post_action_read_new_line
    use mod_hibrid2, only: enord
    use mod_parpot, only: potnam=>pot_name, label=>pot_label
    use mod_file, only: input, output, jobnam
    use mod_par, only: pcod
    use mod_hiutil, only: get_token

    ! show all parameters and flags
    ! show
    class(show_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    integer :: j
    integer :: leninp, lenjob, lenout
    logical :: jtrunc
    character*40 :: code
    character(len=K_MAX_USER_LINE_LENGTH) :: answer

    UNUSED_DUMMY(this)
    call set_param_names(lpar(LPAR_BOUNDC),pcod,icode)
    code = statement_parser%get_token(equal_is_delimiter=.false.)

    if (.not.lpar(LPAR_BOUNDC)) then
      write(6,710) &
       'Parameters (scattering):',(pcod(j),ipar(j),j = 1,iicode)
    else
      write(6,710) &
       'Parameters (bound-state):',(pcod(j),ipar(j),j = 1,iicode)
    endif
    write(6,720)   (pcod(iicode+j),rpar(j),j = 1,ircode-1)
    write(6,1720)  pcod(iicode+ircode), rpar(ircode)
    if(nnout.ne.0) then
      if (.not.lpar(LPAR_TWOMOL) ) then
        write (6,737) 'NOUT: ',nnout, &
               '; JOUT:',(jout(j), j=1,iabs(nnout) )
      else
        write (6,739) 'NOUT: ',nnout, &
             '; J1/J2-OUT: ', &
            (jout(j)/10, mod (jout(j),10), j = 1,iabs(nnout))
      end if
    end if
    if(niout.ne.0) write (6,701) 'INDOUT:',(indout(j), j=1,niout)
    10 format((a))
    701   format(1x,(a),10i5,/,8x,10i5,/,5x,10i5)
    call print_main_params(out_unit=6)
    write(6,731)  nmax, nlammx
    if (ipar(IPAR_LSCREEN) .le. 24 .and. .not. batch) then
      write (6, 703)
    703   format (6x,'enter <return> to continue,', &
                 ' <q> for prompt')
      read (5, 10) answer
      if (answer(1:1) .eq. 'q' .or. answer(1:1) .eq. 'q') then  ! fixme: both 'or' conditions are identical
        post_action = k_post_action_read_new_line
        return
      end if
    end if
    call enord(energ,ipar(IPAR_NERG))
    if(ipar(IPAR_NERG).gt.0) write(6,740) (energ(j),j = 1,ipar(IPAR_NERG))
    710 format(5x,'*** ',(a)/(4(1x,a7,'=',i4,7x)))
    720 format(4(1x,a7,'=',1pg11.4))
    1720 format(1x,a7,'=',f10.5)
    731 format(1x,'** Maximum Channels: ', i4, '; ', &
      'Maximum Anisotropic Terms: ',i5)
    737   format(1x,(a),i2,(a),20(i3,1x))
    739   format (1x,(a),i2,(a),20(2i1,'  ') )
    740 format(1x,'** Energies:',(t15,5f15.6))
    lenout = index(output,' ')-1
    lenjob = index(jobnam,' ')-1
    jtrunc = .false.
    if (lenjob.gt.8) then
       jobnam=jobnam(1:8)
       jtrunc = .true.
    endif
    leninp=index(input,' ')-1
    write (6, 751) label
    751 format (1x,'** Label:      ',(a))
    write (6, 752) potnam
    752 format (1x,'** Pot name:      ',(a))
    if (.not. jtrunc) then
       write (6, 753) input(1:leninp), &
         output(1:lenout), jobnam(1:lenjob)
    753 format(1x,'** Input File:  ',(a),/ &
           1x,'** Output file: ',(a),/,1x,'** Jobname:     ',(a))
    else
       write (6, 754) input(1:leninp), &
         output(1:lenout), jobnam(1:lenjob)
    754 format(1x,'** Input File:  ',(a),/ &
           1x,'** Output file: ',(a),/,1x,'** Jobname:     ',(a), &
           ' (** TRUNCATED TO 8 CHARACTERS **)')
    endif
    post_action = k_post_action_interpret_next_statement
  end subroutine show_execute

  subroutine stmix_execute(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    use mod_histmix, only: stmix
    ! singlet-triplet collisional mixing - added by p. dagdigian
    use mod_hiutil, only: get_token, lower, upper, assignment_parse
    use mod_file, only: jobnam
    use mod_command, only: k_post_action_read_new_line

    class(stmix_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    character(len=40) :: fnam1  ! 1st smt file name
    character(len=40) :: fnam2  ! 2nd smt file name
    real(8) :: a(7)  ! 7 real arguments:
    ! 1: ienerg1: ienerg for 1st smt file
    ! 2: ienerg2: ienerg for 2nd smt file
    ! 3: dele
    ! 4: emax
    ! 5: istats
    ! 6: istatt
    ! 7: wso

    character*8 empty_var_list(0)
    character(len=40) :: code
    integer :: i, j

    UNUSED_DUMMY(this)
    fnam1 = statement_parser%get_token(equal_is_delimiter=.false.)
    if(fnam1 .eq. ' ') fnam1 = jobnam
    call lower(fnam1)
    call upper(fnam1(1:1))
    ! get iener for 1st smt file
    a(1) = 0.d0
    if(statement_parser%statement_end_reached()) goto 3005
    code = statement_parser%get_token(equal_is_delimiter=.false.)
    call assignment_parse(code,empty_var_list,j,a(1))
    3005 fnam2 = statement_parser%get_token(equal_is_delimiter=.false.)
    if(fnam2 .eq. ' ') fnam2 = jobnam
    call lower(fnam2)
    call upper(fnam2(1:1))
    ! get iener for 2nd smt file
    a(2) = 0.d0
    if(statement_parser%statement_end_reached()) goto 3010
    code = statement_parser%get_token(equal_is_delimiter=.false.)
    call assignment_parse(code,empty_var_list,j,a(2))
    ! get dele, emax, istats, istatt, wso
    3010 do 3020 i = 3, 7
      a(i) = 0.d0
      if(statement_parser%statement_end_reached()) goto 3020
      code = statement_parser%get_token(equal_is_delimiter=.false.)
      call assignment_parse(code,empty_var_list,j,a(i))
    3020 continue
    call stmix(fnam1,fnam2,a)
    post_action = k_post_action_read_new_line
  end subroutine stmix_execute

  subroutine sysconf_execute(this, statement_parser, post_action)
    !  print out system parameters
    use mod_statement_parser, only: statement_parser_type
    use mod_command, only: k_post_action_read_new_line
    use mod_hiutil, only: sys_conf

    class(sysconf_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    UNUSED_DUMMY(this)
    UNUSED_DUMMY(statement_parser)
    call sys_conf()
    post_action = k_post_action_read_new_line
  end subroutine sysconf_execute

  subroutine trnprt_execute(this, statement_parser, post_action)
    ! transport cross sections - added by p. dagdigian
    ! TRNPRT,JOB,IENERG,IN1,IN2,JTOTMX,JMIN,JMAX
    use mod_statement_parser, only: statement_parser_type
    use mod_hiutil, only: get_token, lower, upper, assignment_parse
    use mod_file, only: jobnam
    use mod_hitrnprt, only: trnprt
    use mod_command, only: k_post_action_read_new_line

    class(trnprt_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    real(8) :: a(6)  ! 6 real arguments
    character*8 empty_var_list(0)
    character(len=40) :: fnam1
    character(len=40) :: code
    integer :: j

    UNUSED_DUMMY(this)
    fnam1 = statement_parser%get_token(equal_is_delimiter=.false.)
    if(fnam1 .eq. ' ') fnam1 = jobnam
    call lower(fnam1)
    call upper(fnam1(1:1))
    ! get iener for 1st smt file
    a(1) = 0.d0
    if(.not. statement_parser%statement_end_reached()) then
      code = statement_parser%get_token(equal_is_delimiter=.false.)
      call assignment_parse(code,empty_var_list,j,a(1))
      ! todo: get in1, in2, jtotmx, join, jmax
    end if
    call trnprt(fnam1, a)
    post_action = k_post_action_read_new_line
  end subroutine trnprt_execute

  subroutine turn_execute(this, statement_parser, post_action)
    ! determine turning point from isotropic potential
    use mod_statement_parser, only: statement_parser_type
    use mod_hiutil, only: get_token, lower, upper, assignment_parse
    use mod_hibrid1, only: find_turning_point
    use mod_hinput_state, only: irpot
    use mod_command, only: k_post_action_read_new_line, k_post_action_interpret_next_statement

    class(turn_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action

    real(8) :: e
    real(8) :: r
    integer :: i
    logical :: success

    UNUSED_DUMMY(this)
    UNUSED_DUMMY(statement_parser)

    if(irpot .eq. 0) then
      write(6,510)  ! potentiel not yet defined
      510   format(' Potential not yet defined!')
      post_action = k_post_action_read_new_line
      return
    end if

    e = get_max_energy(success)
    if (.not. success) then
      post_action = k_post_action_read_new_line
      return
    end if

    r = find_turning_point(e)
    write(6,1620) r
    1620 format(' Turning point for isotropic potential at R = ', f5.2, ' bohr')
    post_action = k_post_action_interpret_next_statement
  end subroutine turn_execute

  subroutine init()
    use mod_assert, only: fassert
    class(command_type), allocatable :: com
    if (.not. associated(command_mgr)) then
      allocate(command_mgr)
      command_mgr%num_commands = 0
    end if

    com = check_command_type()
    call command_mgr%register_command('CHECK', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = eadiab_command_type()
    call command_mgr%register_command('EADIAB', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = flux_command_type()
    call command_mgr%register_command('FLUX', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = hypxsc_command_type()
    call command_mgr%register_command('HYPXSC', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = indout_command_type()
    call command_mgr%register_command('INDOUT', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = input_command_type()
    call command_mgr%register_command('INPUT', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = intcrs_command_type()
    call command_mgr%register_command('INTCRS', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = j1j2_command_type()
    call command_mgr%register_command('J1J2', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = job_command_type()
    call command_mgr%register_command('JOB', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = label_command_type()
    call command_mgr%register_command('LABEL', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = output_command_type()
    call command_mgr%register_command('OUTPUT', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = partc_command_type()
    call command_mgr%register_command('PARTC', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = prsbr_command_type()
    call command_mgr%register_command('PRSBR', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = read_command_type()
    call command_mgr%register_command('READ', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = run_command_type()
    call command_mgr%register_command('RUN', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = show_command_type()
    call command_mgr%register_command('SHOW', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = stmix_command_type()
    call command_mgr%register_command('STMIX', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = sysconf_command_type()
    call command_mgr%register_command('SYSCONF', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    com = trnprt_command_type()
    call command_mgr%register_command('TRNPRT', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    ASSERT(associated(command_mgr))
  end subroutine init

end module mod_hicommands
