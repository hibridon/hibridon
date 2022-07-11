#include "assert.h"
#include "command.inc.F90"

module mod_hicommands
use mod_command, only: command_type, command_mgr_type
implicit none
  class(command_mgr_type), pointer :: command_mgr  ! singleton

  ! showpot
  type, extends(command_type) :: showpot_command_type
  contains
    procedure :: execute => showpot_execute
  end type showpot_command_type

contains

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

  subroutine init()
    class(command_type), allocatable :: com
    if (.not. associated(command_mgr)) then
      allocate(command_mgr)
      command_mgr%num_commands = 0
    end if
    com = showpot_command_type()
    call command_mgr%register_command('SHOWPOT', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    ASSERT(associated(command_mgr))
  end subroutine init

end module mod_hicommands
