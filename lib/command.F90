#include "assert.h"
#include "command.inc.F90"

! code inspired from https://stackoverflow.com/questions/58749579/what-is-the-canonical-way-to-allocate-and-construct-polymorphic-objects-in-fortr

module mod_command
  implicit none

  integer, parameter :: k_max_num_commands = 100

  type, abstract :: command_type
  contains
    procedure(command_execute_interface), deferred :: execute
  end type command_type

  abstract interface
    subroutine command_execute_interface(this, user_input_line, boca)
    import command_type
      class(command_type) :: this
      character(len=K_MAX_USER_LINE_LENGTH), intent(in) :: user_input_line  ! the string containing one line of a user input eg 'jtot=42;run;prints,myjob,1,3,2,4,7'
      integer, intent(in) :: boca  ! beginning of command argument (index inuser_input_line  of the first character of the first argument of the current command)
    end subroutine
  end interface

  type command_item_type
    character(len=8) :: codex
    class(command_type), allocatable :: item
  end type command_item_type

  type command_mgr_type
    integer :: num_commands
    type(command_item_type) :: commands(3)
  contains
    procedure :: register_command => commands_register_command
    procedure :: execute_command => commands_execute_command
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
    ! allocate(this%commands(this%num_commands))
   this%commands(this%num_commands)%codex = codex
    this%commands(this%num_commands)%item = command
  end subroutine commands_register_command

  subroutine commands_execute_command(this, command_statement)
    class(command_mgr_type), intent(inout) :: this
    character(len=*), intent(in) :: command_statement
    integer :: i
    character(len=K_MAX_USER_LINE_LENGTH) :: line 
    integer :: boca = 1
    do i=1, this%num_commands
      if (this%commands(i)%codex == command_statement) then
        call this%commands(i)%item%execute(user_input_line=line, boca=boca)
      end if
    end do
  end subroutine commands_execute_command

  function command_mgr_ctor() result(command_mgr)
    class(command_mgr_type), pointer :: command_mgr
    class(command_mgr_type), allocatable, target :: mgr
    command_mgr => mgr

  end function command_mgr_ctor

end module mod_command


