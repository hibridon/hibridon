#include "command.inc.F90"
#include "unused.h"

module mod_unitcommands
  use mod_command, only: command_type
  use mod_statement_parser, only: statement_parser_type

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

  subroutine dummyc1_execute(this, statement_parser, post_action)
    use mod_command, only: k_post_action_read_new_line
    class(dummyc1_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action
    UNUSED_DUMMY(this)
    UNUSED_DUMMY(statement_parser)
    write (6,*) 'executing dummyc1'
    post_action = k_post_action_read_new_line
  end subroutine

  subroutine dummyc2_execute(this, statement_parser, post_action)
    use mod_command, only: k_post_action_read_new_line
    class(dummyc2_command_type) :: this
    class(statement_parser_type), intent(inout) :: statement_parser
    integer, intent(out) :: post_action
    UNUSED_DUMMY(this)
    UNUSED_DUMMY(statement_parser)
    write (6,*) 'executing dummyc2'
    post_action = k_post_action_read_new_line
  end subroutine

end module mod_unitcommands

module mod_command_unit_test
contains
subroutine test_commands()
  use mod_unitcommands, only: dummyc2_command_type, dummyc1_command_type
  use mod_command, only: command_type, command_mgr_type
  use mod_statement_parser, only: statement_parser_type

  type(statement_parser_type) :: statement_parser
  class(command_mgr_type), allocatable :: command_mgr
  class(command_type), allocatable :: com
  integer :: post_action

  allocate(command_mgr)
  command_mgr%num_commands = 0

  com = dummyc2_command_type()
  call command_mgr%register_command('DUMMYC2', com)

  deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free
  com = dummyc1_command_type()
  call command_mgr%register_command('DUMMYC1', com)

  call statement_parser%initialize()

  write(6,*) 'after register num_commands=', command_mgr%num_commands
  write(6,*) 'after register, 1st codex is ', command_mgr%commands(1)%codex

  write (6, *) 'after init, 1st codex is ', command_mgr%commands(1)%codex
  write (6, *) 'after init, 2nd codex is ', command_mgr%commands(2)%codex
  call command_mgr%commands(1)%item%execute(statement_parser=statement_parser, post_action=post_action)
  call command_mgr%commands(2)%item%execute(statement_parser=statement_parser, post_action=post_action)


  call command_mgr%execute_command('DUMMYC1', post_action)
  call command_mgr%execute_command('DUMMYC2', post_action)
end subroutine test_commands

end module mod_command_unit_test


program unittest_commands
use mod_command_unit_test, only: test_commands
call test_commands()
        
end program unittest_commands