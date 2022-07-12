#include "command.inc.F90"

module mod_unitcommands
  use mod_command, only: command_type

  ! dummyc1
  type, extends(command_type) :: dummyc1_command_type
  contains
    procedure :: execute => dummyc1_execute
  end type dummyc1_command_type

  ! dummyc2
  type, extends(command_type) :: dummyc2_command_type
  contains
    procedure :: execute => dummyc2_execute
  end type dummyc2_command_type

contains

  subroutine dummyc1_execute(this, statements, bofargs, next_statement, post_action)
    class(dummyc1_command_type) :: this
    character(len=K_MAX_USER_LINE_LENGTH), intent(in) :: statements
    integer, intent(in) :: bofargs
    integer, intent(out) :: next_statement
    integer, intent(out) :: post_action
    write (6,*) 'executing dummyc1'
  end subroutine

  subroutine dummyc2_execute(this, statements, bofargs, next_statement, post_action)
    class(dummyc2_command_type) :: this
    character(len=K_MAX_USER_LINE_LENGTH), intent(in) :: statements
    integer, intent(in) :: bofargs
    integer, intent(out) :: next_statement
    integer, intent(out) :: post_action
    write (6,*) 'executing dummyc2'
  end subroutine

end module mod_unitcommands


subroutine test_commands()
  use mod_unitcommands, only: dummyc2_command_type, dummyc1_command_type
  use mod_command, only: command_type, command_mgr_type

  character(len=K_MAX_USER_LINE_LENGTH) :: line
  class(command_mgr_type), allocatable :: command_mgr
  class(command_type), allocatable :: com
  integer :: next_statement
  integer :: post_action

  allocate(command_mgr)
  command_mgr%num_commands = 0

  com = dummyc2_command_type()
  call command_mgr%register_command('DUMMYC2', com)

  deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free
  com = dummyc1_command_type()
  call command_mgr%register_command('DUMMYC1', com)

  write(6,*) 'after register num_commands=', command_mgr%num_commands
  write(6,*) 'after register, 1st codex is ', command_mgr%commands(1)%codex

  write (6, *) 'after init, 1st codex is ', command_mgr%commands(1)%codex
  write (6, *) 'after init, 2nd codex is ', command_mgr%commands(2)%codex
  call command_mgr%commands(1)%item%execute(statements=line, bofargs=1, next_statement=next_statement, post_action=post_action)
  call command_mgr%commands(2)%item%execute(statements=line, bofargs=1, next_statement=next_statement, post_action=post_action)


  call command_mgr%execute_command('DUMMYC1', post_action)
  call command_mgr%execute_command('DUMMYC2', post_action)
end subroutine test_commands

program unittest_commands

call test_commands()
        
end program unittest_commands