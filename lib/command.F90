#include "assert.h"
#include "command.inc.F90"

! code inspired from https://stackoverflow.com/questions/58749579/what-is-the-canonical-way-to-allocate-and-construct-polymorphic-objects-in-fortr

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
    subroutine command_execute_interface(this, statements, bofargs, next_statement, post_action)
    import command_type
      class(command_type) :: this
      character(len=K_MAX_USER_LINE_LENGTH), intent(in) :: statements  ! the string containing one line of a user input eg 'jtot=42;run;prints,myjob,1,3,2,4,7'
      integer, intent(in) :: bofargs  ! beginning of command argument (index in statements  of the first character of the first argument of the current command)
      integer, intent(out) :: next_statement  ! index in statements of the beginning of next statement, in case post_action == k_post_action_interpret_next_statement
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
  interface command_mgr_type
    module procedure command_mgr_ctor
  end interface command_mgr_type

contains

  subroutine commands_register_command(this, codex, command)
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
    class(command_mgr_type), intent(inout) :: this
    character(len=*), intent(in) :: command_statement
    integer, intent(out) :: post_action
    integer :: i
    character(len=K_MAX_USER_LINE_LENGTH) :: line 
    integer :: bofargs = 1
    integer :: next_statement
    do i=1, this%num_commands
      if (this%commands(i)%codex == command_statement) then
        call this%commands(i)%item%execute(statements=line, bofargs=bofargs, next_statement=next_statement, post_action=post_action)
        return
      end if
    end do
  end subroutine commands_execute_command

  function commands_get_command(this, command_name) result(command)
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

  function command_mgr_ctor() result(command_mgr)
    class(command_mgr_type), pointer :: command_mgr
    class(command_mgr_type), allocatable, target :: mgr
    command_mgr => mgr

  end function command_mgr_ctor

end module mod_command


