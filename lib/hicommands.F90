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

  ! prsbr
  type, extends(command_type) :: prsbr_command_type
  contains
    procedure :: execute => prsbr_execute
  end type prsbr_command_type

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

  subroutine prsbr_execute(this, user_input_line, boca)
    use mod_hiutil, only: assignment_parse
    class(prsbr_command_type) :: this
    character(len=K_MAX_USER_LINE_LENGTH), intent(in) :: user_input_line  ! the string containing one line of a user input eg 'jtot=42;run;prints,myjob,1,3,2,4,7'
    integer, intent(in) :: boca  ! beginning of command argument (index inuser_input_line  of the first character of the first argument of the current command)
    ! pressure broadening cross sections - added by p. dagdigian
    integer :: l
    character(len=40) :: fnam1, fnam2
    character*8 empty_var_list(0)
    character*40 :: code
    integer :: i, j, lc
    integer, parameter :: k_num_args = 12
    real(8) :: a(k_num_args)

    common /cofile/ input, output, jobnam, savfil
    character*40 :: input
    character*40 :: output
    character*40 :: jobnam
    character*40 :: savfil


    l = boca

    call get_token(user_input_line,l,fnam1,lc)
    if(fnam1 .eq. ' ') fnam1 = jobnam
    call lower(fnam1)
    call upper(fnam1(1:1))
    ! get iener for 1st smt file
    a(1) = 0.d0
    if(l .eq. 0) goto 3205
    call get_token(user_input_line,l,code,lc)
    call assignment_parse(code(1:lc),empty_var_list,j,a(1))
    3205 call get_token(user_input_line,l,fnam2,lc)
    if(fnam2 .eq. ' ') fnam2 = jobnam
    call lower(fnam2)
    call upper(fnam2(1:1))
    ! get iener for 2nd smt file
    a(2) = 0.d0
    if(l .eq. 0) goto 3210
    call get_token(user_input_line,l,code,lc)
    call assignment_parse(code(1:lc),empty_var_list,j,a(2))
    ! get k, j1, in1, j2, in2, diag, j1p, in1p, j2p, in2p
    3210 do 3220 i = 3, 12
      a(i) = 0.d0
      if(l .eq. 0) goto 3220
      call get_token(user_input_line,l,code,lc)
      call assignment_parse(code(1:lc),empty_var_list,j,a(i))
    3220 continue
    call prsbr(fnam1,fnam2,a)
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

    com = prsbr_command_type()
    call command_mgr%register_command('PRSBR', com)
    deallocate(com)  ! without deallocation, address sanitizer would detect a heap-use-after-free if com is reused

    ASSERT(associated(command_mgr))
  end subroutine init

end module mod_hicommands
