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
    !procedure :: register_command => commands_register_command
  end type command_mgr_type

  ! interface to create an instance of command_mgr_type using a construct familiar to other languages:
  ! g1 = command_mgr_type()
  interface command_mgr_type
    module procedure create_command_mgr_type
  end interface command_mgr_type

contains

  ! subroutine commands_register_command(this, command)
  !   class(command_mgr_type), intent(inout) :: this
  !   class(command_type), intent(in), allocatable :: command
  !   ! class(command_type), pointer :: command_ptr
  !   integer :: num_commands
  !   this%num_commands = this%num_commands + 1
  !   ASSERT( this%num_commands <= k_max_num_commands )
  !   !this%toto => command
  !   num_commands = this%num_commands
  !   ! command_ptr => this%commands(num_commands)%
  !   write(6,*) 'command codex : ',command%codex
  !   ! allocate(this%toto)
  !   ! this%toto%p => command
  !   !allocate(this%commands(num_commands))
  !   !this%commands(num_commands) = command
  ! end subroutine commands_register_command

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

  subroutine init()
    if (.not. associated(command_mgr)) then
      allocate(command_mgr)
      command_mgr%num_commands = 0
      write (6,*)  'coucou init: manager has been allocated'
    end if

    !allocate(showpot_command_type :: command_mgr%commands(1))
    command_mgr%commands(1)%item = showpot_command_type()
    ! call command_mgr%register_command( command )
    write(6,*) 'after register num_commands=', command_mgr%num_commands
    write(6,*) 'after register, 1st codex is ', command_mgr%commands(1)%item%codex
    ASSERT(associated(command_mgr))
  end subroutine init

end module mod_hicommands

subroutine test_commands()
  use mod_hicommands, only: init
  call init()
end subroutine test_commands

subroutine test_polyarray()
  implicit none

  type basetype
  end type basetype

  type, extends(basetype) :: exttype1
  end type exttype1

  type, extends(exttype1) :: exttype2
  end type exttype2

  type arraytype
    class(basetype), allocatable :: comp
  end type arraytype

  type(arraytype), dimension(3) :: ary

  integer :: i

  allocate (basetype::ary(1)%comp)
  allocate (exttype1::ary(2)%comp)
  allocate (exttype2::ary(3)%comp)

  do i=1,3
  select type (st=>ary(i)%comp)
  type is (basetype)
  print 101,i,"basetype"
  type is (exttype1)
  print 101,i,"exttype1"
  type is (exttype2)
  print 101,i,"exttype2"
  class default
  print 101,i,"unknown"
  end select
  101    format ("The dynamic type of ary(",i0,")%comp is ",A)
  end do
end subroutine test_polyarray

program PolyArray

! call test_polyarray()
call test_commands()
        
end program PolyArray