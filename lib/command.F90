#include "assert.h"
#include "command.inc.F90"

! code inspired from https://stackoverflow.com/questions/58749579/what-is-the-canonical-way-to-allocate-and-construct-polymorphic-objects-in-fortr

module mod_statement_parser

  type, public :: statement_parser_type
    character(len=K_MAX_USER_LINE_LENGTH) current_line ! current statement line being parsed
    integer        :: current_pos = 0 ! points to first character to be searched in current_line
    contains

    procedure                  :: initialize => statement_parser_type_constructor
    final                      :: statement_parser_type_destructor

    procedure                  :: get_token => statement_parser_type_get_token
    procedure                  :: statement_end_reached => statement_parser_statement_end_reached
    procedure                  :: prev_char_is => statement_parser_prev_char_is
  end type statement_parser_type

  ! interface to create an instance of statement_parser_type using a construct familiar to other languages:
  ! g1 = statement_parser_type()
  ! interface statement_parser_type
  !   module procedure create_statement_parser_type
  ! end interface statement_parser_type

contains

! function create_statement_parser_type() result(statement_parser)

!   class(statement_parser_type), allocatable :: statement_parser
!   allocate(statement_parser)
!   statement_parser%current_line = ''
!   statement_parser%current_pos = 0
! end function

! constructor for statement_parser_type
subroutine statement_parser_type_constructor(this)
  class(statement_parser_type) :: this
  this%current_line = ''
  this%current_pos = 0
end subroutine statement_parser_type_constructor

! destructor for statement_parser_type
subroutine statement_parser_type_destructor(this)
  type(statement_parser_type) :: this
end subroutine statement_parser_type_destructor

function statement_parser_type_get_token(this, equal_is_delimiter) result(token)
  use mod_hiutil, only: get_token
  use mod_assert, only: fassert
  implicit none
  class(statement_parser_type), intent(inout) :: this
  logical, intent(in) :: equal_is_delimiter

  character(len=:), allocatable :: token

  character(len=40) :: tmp_token
  integer :: token_length
  integer :: current_pos

  current_pos = this%current_pos
  if (equal_is_delimiter) then
    current_pos = -current_pos
  end if

  call get_token(this%current_line, current_pos, tmp_token, token_length)
  allocate(character(len=token_length) :: token)
  token = tmp_token(1:token_length)
  ASSERT(current_pos >= 0)
  this%current_pos = current_pos
end function

function statement_parser_statement_end_reached(this)
  class(statement_parser_type), intent(in) :: this
  logical :: statement_parser_statement_end_reached

  statement_parser_statement_end_reached = (this%current_pos == 0)
end function

function statement_parser_prev_char_is(this, char)
  use mod_assert, only: fassert
  class(statement_parser_type), intent(in) :: this
  character(len=1), intent(in) :: char

  logical :: statement_parser_prev_char_is

  integer :: prev_char_index

  ASSERT(this%current_pos > 1)
  prev_char_index = this%current_pos-1
  statement_parser_prev_char_is = (this%current_line(prev_char_index:prev_char_index) .eq. char)
end function

end module mod_statement_parser


module mod_command
  implicit none

  integer, parameter :: k_max_num_commands = 100

  ! post_actions indicate what should be done after a command has been executed
  integer, parameter :: k_post_action_read_new_line = 1  ! ask the user a new statements line then interpret it
  integer, parameter :: k_post_action_interpret_next_statement = 2  ! interprets the next statement in statments line
  integer, parameter :: k_post_action_write_cr_and_exit = 3  ! write carriage return and exit hinput function
  integer, parameter :: k_post_action_exit_hinput = 4  ! exit hinput function
  integer, parameter :: k_post_action_exit_hibridon = 5  ! terminate hibridon program
  
  
  type, abstract :: command_type
  contains
    procedure(command_execute_interface), deferred :: execute
  end type command_type

  abstract interface
    subroutine command_execute_interface(this, statement_parser, post_action)
    use mod_statement_parser, only: statement_parser_type
    import command_type
      class(command_type) :: this
      class(statement_parser_type), intent(inout) :: statement_parser
      ! intent(in) :: statement_parser%current_line  ! the string containing one line of a user input eg 'jtot=42;run;prints,myjob,1,3,2,4,7'
      ! intent(in) :: statement_parser%current_pos  ! beginning of command argument (index in statements  of the first character of the first argument of the current command)
      ! intent(out) :: statement_parser%current_pos  ! index in statements of the beginning of next statement, in case post_action == k_post_action_interpret_next_statement
      integer, intent(out) :: post_action  ! indicates what the caller is expected to do once the command has been executed
    end subroutine
  end interface

  type command_item_type
    character(len=8) :: codex
    class(command_type), allocatable :: item
  end type command_item_type

  type command_mgr_type
    integer :: k_max_num_commands = 100
    integer :: num_commands
    type(command_item_type) :: commands(k_max_num_commands)
  contains
    procedure :: register_command => commands_register_command
    procedure :: execute_command => commands_execute_command
    procedure :: get_command => commands_get_command
  end type command_mgr_type

  ! interface to create an instance of command_mgr_type using a construct familiar to other languages:
  ! g1 = command_mgr_type()
  ! interface command_mgr_type
  !   module procedure command_mgr_ctor
  ! end interface command_mgr_type

contains

  subroutine commands_register_command(this, codex, command)
    use mod_assert, only: fassert
    class(command_mgr_type), intent(inout) :: this
    character(len=*), intent(in) :: codex
    class(command_type), intent(in), allocatable :: command
    this%num_commands = this%num_commands + 1
    ASSERT(this%num_commands <= this%k_max_num_commands)
    ! allocate(this%commands(this%num_commands))
    this%commands(this%num_commands)%codex = codex
    this%commands(this%num_commands)%item = command
  end subroutine commands_register_command

  subroutine commands_execute_command(this, command_statement, post_action)
    use mod_statement_parser, only: statement_parser_type
    class(command_mgr_type), intent(inout) :: this
    character(len=*), intent(in) :: command_statement
    integer, intent(out) :: post_action
    integer :: i
    type(statement_parser_type) :: statement_parser

    call statement_parser%initialize()
    statement_parser%current_line = command_statement
    do i=1, this%num_commands
      if (this%commands(i)%codex == command_statement) then
        call this%commands(i)%item%execute(statement_parser, post_action=post_action)
        return
      end if
    end do
  end subroutine commands_execute_command

  function commands_get_command(this, command_name) result(command)
    use mod_assert, only: fassert
    class(command_mgr_type), intent(inout), target :: this
    character(len=*), intent(in) :: command_name
    class(command_type), pointer :: command
    integer :: i
    do i=1, this%num_commands
      if (this%commands(i)%codex == command_name) then
        command => this%commands(i)%item
        return
      end if
    end do
    ASSERT(.false.)  ! unable to find the given command
    command => null()
  end function commands_get_command

  ! function command_mgr_ctor() result(command_mgr)
  !   class(command_mgr_type), pointer :: command_mgr
  !   class(command_mgr_type), allocatable, target :: mgr
  !   command_mgr => mgr  ! disable-warnings:target-lifetime

  ! end function command_mgr_ctor

end module mod_command


