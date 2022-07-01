! investigating solutions to replace log code such as :

!   write (6,16) rmu * xmconv, ered * econv, jtot, jlpar
!   write (9,16) rmu * xmconv, ered * econv, jtot, jlpar
! 16     format(/' **  2/2 ATOM ON UNCORRUGATED SURFACE ** RMU=',f9.4, &
!       '             E=', f9.2, /, ' JTOT=', i5, 2x,' JLPAR=',i2)

! problems:
!   1. the labels pollute the view, making the code structure harder to read
!   2. the arguments list is duplicated, making the source code harder to maintain
!   3. the arguments list is duplicated, potentially impacting performance (same expressions is evaluated more than once)
!   4. the logging code is not as compact as it could, thus polluting the view by making source code fatter: at least 3 lines are needed while 1 could be sufficent

! solution1
!   benefits:
!   - (1) removes the labels, thus making code easier to read
!   - (4) 2 lines of source code instead of 3
!   inconvenients:
!   - COMMA instead of ',' makes argument list harder to read
!   - duplicates the format in memory
!   - duplicates the format and arg list in source code, making it harder to maintain
!
! solution2
!
!   writes to a string array (because a string would not accept the '/' format), then sends the string to the output
!
!   benefits:
!   - (1) removes the labels, thus making code easier to read
!   - (3) the values are only evaluated once (increase performance over the existing code)
!   - (4) 1 lines of source code instead of 3
!   inconvenients:
!   - COMMA instead of ',' makes argument list harder to read
!   - need to declare USE_LOG2 in each subroutine
!
! solution3
!
!   simply embeds the original log messages into a preprocessing macro which just hides the label
!
!   benefits:
!   - (1) removes the labels from the view (they still exist), thus making code easier to read
!   - (4) 1 lines of source code instead of 3
!   inconvenients:
!   - COMMA instead of ',' makes argument list harder to read
!   - COMMA instead of ',' makes also format harder to read
!   - a bit less generic regarding units : needs one macro for each length of unit array
!
! solution4
!
!   solution similar to solution2 that doesn't make use of fortran preprocessor
!   Uses the special dt formatter mechanism inspired by https://stackoverflow.com/questions/19906446/how-to-wrap-the-fortran-write-statement
!   but this mechanism doesn't bring much benefit (we could have a similar simpler solution without it)
!
!   benefits:
!   - (1) removes the labels, thus making code easier to read
!   - (3) the values are only evaluated once (increase performance over the existing code)
!   - (4) 2 lines of source code instead of 3
!   inconvenients:
!   - need to declare a message variable in each subroutine
!
#define COMMA ,
#define __NL__  NEWLINE
#define MYWRITE1(UNIT, FMT, X) write(UNIT,FMT) X

#define USE_LOG2 use mod_write, only: tmp_message, log, std_log_units
#define LOG2(UNITS, FMT, X) tmp_message(:)(1:1) = char(0); write(tmp_message,FMT) X ; call log(UNITS, tmp_message)
#define LOG3_1(UNIT1, FMT, X) __LINE__ format(FMT); write(UNIT1,__LINE__) X
#define LOG3_2(UNIT1, UNIT2, FMT, X) __LINE__ format(FMT); write(UNIT1,__LINE__) X; write(UNIT2,__LINE__) X

module mod_write
  integer, allocatable :: std_log_units(:)
  ! writing to an unallocated string is not supported https://community.intel.com/t5/Intel-Fortran-Compiler/internal-formatted-write-to-deferred-length-allocatable/td-p/751843
  character(len=:), allocatable :: tmp_message(:)
  integer, parameter :: max_lines = 10
  contains
  subroutine init()
    allocate(integer :: std_log_units(2))
    std_log_units(1) = 6
    std_log_units(2) = 9
    allocate(character(2048) :: tmp_message(max_lines))
  end subroutine init

  subroutine log(log_units, message)
    implicit none
    integer, allocatable, intent(in) :: log_units(:)
    character(len=:), allocatable, intent(in) :: message(:)
    integer :: num_lines
    integer :: char
    integer :: iline, iunit

    num_lines = 0
    do iline = 1, max_lines
      char = iachar(message(iline)(1:1))
      ! write(6,'(i3)') char 
      if (char /= 0) then
        num_lines = num_lines + 1
      end if
    end do

    do iunit = 1, size(log_units)
      do iline = 1, num_lines
        !write(6,*) iline
          !write(6,'(a)',advance="no") trim(message(iline))
        write(log_units(iunit),'(a)') trim(message(iline))
      end do
    end do
  end subroutine
end module mod_write

module logger4_mod

  integer, parameter :: max_lines = 10

  type message_t
    character(len=2048) :: lines(max_lines)
  contains
    procedure message_write
    generic :: write(formatted) => message_write
  end type message_t

  logical :: debug_verbose = .true.

contains

  subroutine log(units, message)
    implicit none
    integer, intent(in) :: units(2)
    class(message_t) :: message
    integer :: iunit
    do iunit = 1, size(units)
      ! write(6,*) 'writing to unit', units(iunit)
      write(unit=units(iunit), fmt='(dt)') message
    end do
    message%lines(:)(1:1) = char(0) ! empty all lines in the message because the next write into it might leave old trailing lines
  end subroutine

  subroutine message_write(message, unit, iotype, v_list, iostat, iomsg)
    implicit none
    class(message_t), intent(in)  :: message
    integer, intent(in)  :: unit
    character(len=*), intent(in)  :: iotype
    integer, intent(in), dimension(:)  :: v_list
    integer, intent(out)  :: iostat
    character(len=*), intent(inout)  :: iomsg
    integer :: char

    integer :: iline, num_lines

    if (debug_verbose) then
      num_lines = 0
      do iline = 1, max_lines
        char = iachar(message%lines(iline)(1:1))
        ! write(6,'(i3)') char 
        if (char /= 0) then
          num_lines = num_lines + 1
        end if
      end do

      do iline = 1, num_lines
        write(unit,'(a,/)', iostat=iostat, iomsg=iomsg, advance='yes') trim(message%lines(iline))  ! graffy: no idea why I need to add a newline here '/'
      end do
    else
      write(unit, '()', advance='no')
    end if

  end subroutine message_write

end module logger4_mod

subroutine test1_original()
  real(8) :: rmu, xmconv, ered, econv
  integer :: jtot, jlpar
  rmu = 2.0
  xmconv = 3.0
  ered = 4.0
  econv = 5.5
  jtot = 13
  jlpar = 14

  write (6,16) rmu * xmconv, ered * econv, jtot, jlpar
  write (9,16) rmu * xmconv, ered * econv, jtot, jlpar
16     format(/' **  2/2 ATOM ON UNCORRUGATED SURFACE ** RMU=',f9.4, &
      '             E=', f9.2, /, ' JTOT=', i5, 2x,' JLPAR=',i2)

end subroutine test1_original

subroutine test1_solution1()
  real(8) :: rmu, xmconv, ered, econv
  integer :: jtot, jlpar
  rmu = 2.0
  xmconv = 3.0
  ered = 4.0
  econv = 5.5
  jtot = 13
  jlpar = 14
  

  MYWRITE1(6, '(/," **  2/2 ATOM ON UNCORRUGATED SURFACE ** RMU=",f9.4,"             E=",f9.2,/," JTOT=",i5,2x," JLPAR=",i2)', rmu * xmconv COMMA ered * econv COMMA jtot COMMA jlpar)
  MYWRITE1(9, '(/," **  2/2 ATOM ON UNCORRUGATED SURFACE ** RMU=",f9.4,"             E=",f9.2,/," JTOT=",i5,2x," JLPAR=",i2)', rmu * xmconv COMMA ered * econv COMMA jtot COMMA jlpar)

end subroutine test1_solution1

subroutine test1_solution2()
  USE_LOG2
  real(8) :: rmu, xmconv, ered, econv
  integer :: jtot, jlpar
  rmu = 2.0
  xmconv = 3.0
  ered = 4.0
  econv = 5.5
  jtot = 13
  jlpar = 14

  LOG2(std_log_units, '(/," **  2/2 ATOM ON UNCORRUGATED SURFACE ** RMU=",f9.4,"             E=",f9.2,/," JTOT=",i5,2x," JLPAR=",i2)', rmu * xmconv COMMA ered * econv COMMA jtot COMMA jlpar)

end subroutine test1_solution2

subroutine test1_solution3()
  real(8) :: rmu, xmconv, ered, econv
  integer :: jtot, jlpar
  rmu = 2.0
  xmconv = 3.0
  ered = 4.0
  econv = 5.5
  jtot = 13
  jlpar = 14

  LOG3_2(6, 9, / COMMA ' **  2/2 ATOM ON UNCORRUGATED SURFACE ** RMU=' COMMA f9.4 COMMA '             E=' COMMA  f9.2 COMMA  / COMMA  ' JTOT=' COMMA  i5 COMMA  2x COMMA ' JLPAR=' COMMA i2, rmu * xmconv COMMA ered * econv COMMA jtot COMMA jlpar)

end subroutine test1_solution3

subroutine test1_solution4()
  use logger4_mod, only : message_t, log
  implicit none
  type(message_t) message
  real(8) :: rmu, xmconv, ered, econv
  integer :: jtot, jlpar
  rmu = 2.0
  xmconv = 3.0
  ered = 4.0
  econv = 5.5
  jtot = 13
  jlpar = 14

  write(message%lines, '(/," **  2/2 ATOM ON UNCORRUGATED SURFACE ** RMU=",f9.4,"             E=",f9.2,/," JTOT=",i5,2x," JLPAR=",i2)') rmu * xmconv, ered * econv, jtot, jlpar 
  call log((/6, 9/), message)

end subroutine test1_solution4


program test_write
use mod_write, only: init
implicit none

call init()
open(unit=9, file='/tmp/toto.txt', access='append')

call test1_original()
call test1_solution1()
call test1_solution2()
call test1_solution3()
call test1_solution4()

close(9)

end program test_write


