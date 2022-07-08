#include "assert.h"
#define K_MAX_USER_LINE_LENGTH 80

subroutine driver()
write (6,*) 'coucou : executing subroutine driver'
end subroutine driver

module mod_command
  implicit none

  integer, parameter :: k_max_num_commands = 100

  type, abstract :: command_type
    character(len=8) :: codex
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
    class(command_type), allocatable :: item
  end type command_item_type

  type command_mgr_type
    integer :: num_commands
    type(command_item_type) :: commands(3)
  contains
    procedure :: register_command => commands_register_command
  end type command_mgr_type

  ! interface to create an instance of command_mgr_type using a construct familiar to other languages:
  ! g1 = command_mgr_type()
  interface command_mgr_type
    module procedure create_command_mgr_type
  end interface command_mgr_type

contains

  subroutine commands_register_command(this, command)
    class(command_mgr_type), intent(inout) :: this
    class(command_type), intent(in), allocatable :: command
  !   ! class(command_type), pointer :: command_ptr
  !   integer :: num_commands
    this%num_commands = this%num_commands + 1
    this%commands(this%num_commands)%item = command
  !   ASSERT( this%num_commands <= k_max_num_commands )
  !   !this%toto => command
  !   num_commands = this%num_commands
  !   ! command_ptr => this%commands(num_commands)%
  !   write(6,*) 'command codex : ',command%codex
  !   ! allocate(this%toto)
  !   ! this%toto%p => command
  !   !allocate(this%commands(num_commands))
  !   !this%commands(num_commands) = command
  end subroutine commands_register_command

  function create_command_mgr_type() result(command_mgr)
    class(command_mgr_type), pointer :: command_mgr
    class(command_mgr_type), allocatable, target :: mgr
    command_mgr => mgr

    write (6,*) 'coucou'
  end function create_command_mgr_type

end module mod_command

module mod_hicommands
use mod_command, only: command_type, command_mgr_type
implicit none
  class(command_mgr_type), pointer :: command_mgr  ! singleton

  ! showpot
  type, extends(command_type) :: showpot_command_type
  contains
    procedure :: execute => showpot_execute
  end type showpot_command_type

  interface showpot_command_type
    module procedure :: showpot_command_constructor
  end interface showpot_command_type

  ! dummyc1
  type, extends(command_type) :: dummyc1_command_type
  contains
    procedure :: execute => dummyc1_execute
  end type dummyc1_command_type

  interface dummyc1_command_type
    module procedure :: dummyc1_command_constructor
  end interface dummyc1_command_type

contains

  function showpot_command_constructor()
    type(showpot_command_type) :: showpot_command_constructor
    write(6,*) 'showpot_command_constructor'
    showpot_command_constructor%codex='SHOWPOT'
  end function showpot_command_constructor


  subroutine showpot_execute(this, user_input_line, boca)
    class(showpot_command_type) :: this
    character(len=K_MAX_USER_LINE_LENGTH), intent(in) :: user_input_line  ! the string containing one line of a user input eg 'jtot=42;run;prints,myjob,1,3,2,4,7'
    integer, intent(in) :: boca  ! beginning of command argument (index inuser_input_line  of the first character of the first argument of the current command)
    write(6,*) "************************************************************"
    write(6,*) "Entering the DRIVER subroutine of the potential"
    write(6,*) "Press Ctrl+D to go back to Hibridon's console"
    write(6,*) "************************************************************"
    call driver
  end subroutine

  function dummyc1_command_constructor()
    type(dummyc1_command_type) :: dummyc1_command_constructor
    write(6,*) 'dummyc1_command_constructor'
    dummyc1_command_constructor%codex='DUMMYC1'
  end function dummyc1_command_constructor


  subroutine dummyc1_execute(this, user_input_line, boca)
    class(dummyc1_command_type) :: this
    character(len=K_MAX_USER_LINE_LENGTH), intent(in) :: user_input_line  ! the string containing one line of a user input eg 'jtot=42;run;prints,myjob,1,3,2,4,7'
    integer, intent(in) :: boca  ! beginning of command argument (index inuser_input_line  of the first character of the first argument of the current command)
    write (6,*) 'executing dummyc1'
  end subroutine

  subroutine init()
    class(command_type), allocatable :: com
    class(command_type), allocatable :: com1
    if (.not. associated(command_mgr)) then
      allocate(command_mgr)
      command_mgr%num_commands = 0
      write (6,*)  'coucou init: manager has been allocated'
    end if
    com = showpot_command_type()
    !allocate(showpot_command_type :: command_mgr%commands(1))
    call command_mgr%register_command(com)
    com1 = dummyc1_command_type()
    call command_mgr%register_command(com1)
    !command_mgr%commands(1)%item = showpot_command_type()
    ! call command_mgr%register_command( command )
    write(6,*) 'after register num_commands=', command_mgr%num_commands
    write(6,*) 'after register, 1st codex is ', command_mgr%commands(1)%item%codex
    ASSERT(associated(command_mgr))
  end subroutine init

end module mod_hicommands

subroutine test_commands()
  use mod_hicommands, only: init, command_mgr
  character(len=K_MAX_USER_LINE_LENGTH) :: line
  call init()


  write (6, *) 'after init, 1st codex is ', command_mgr%commands(1)%item%codex
  write (6, *) 'after init, 2nd codex is ', command_mgr%commands(2)%item%codex
  call command_mgr%commands(1)%item%execute(user_input_line=line, boca=1)
  call command_mgr%commands(2)%item%execute(user_input_line=line, boca=1)
end subroutine test_commands

program unittest_commands

call test_commands()
        
end program unittest_commands