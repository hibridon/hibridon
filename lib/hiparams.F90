!************************************************************************
!                                                                       *
!                  module param_goup                                    *
!                                                                       *
!************************************************************************
!  handles group of parameters                                          *
!************************************************************************
#include "assert.h"

! declare param_array_type_i
#define T i
#define TT integer
#include "hiparam_array.F90.inc"

! declare param_array_type_r
#define T r
#define TT real(8)
#include "hiparam_array.F90.inc"

! declare param_array_type_l
#define T l
#define TT logical
#include "hiparam_array.F90.inc"

module mod_param_group
  use mod_param_array_type_i, only: param_array_type_i => param_array_type, iparam_type => param_type, iroparam_type => roparam_type
  use mod_param_array_type_r, only: param_array_type_r => param_array_type, rparam_type => param_type, rroparam_type => roparam_type
  use mod_param_array_type_l, only: param_array_type_l => param_array_type, lparam_type => param_type, lroparam_type => roparam_type
  implicit none
  integer, parameter :: kmaxpar = 150

  type :: param_group_type
    type(param_array_type_i) :: iparams
    type(param_array_type_r) :: rparams
    type(param_array_type_l) :: lparams
  !type(param_array_type_i(kmaxpar)) :: iparams
    !type(param_array_type_r(kmaxpar)) :: rparams
    !type(param_array_type_l(kmaxpar)) :: lparams

  contains 
    procedure, public :: create_iparam => param_group_create_iparam
    procedure, public :: create_rparam => param_group_create_rparam
    procedure, public :: create_lparam => param_group_create_lparam
    procedure, public :: get_ivalue => param_group_get_ivalue
    procedure, public :: get_iparam => param_group_get_iparam
    procedure, public :: get_iroparam => param_group_get_iroparam
  end type param_group_type

  type(param_group_type) :: basis_params
 
contains

  function param_group_create_iparam(this, param_name)
    class(param_group_type), intent(inout) :: this
    character(len=*), intent(in) :: param_name
    type(iparam_type) :: param_group_create_iparam

    param_group_create_iparam = this%iparams%create_param(param_name)
  end function

  function param_group_create_rparam(this, param_name)
    class(param_group_type), intent(inout) :: this
    character(len=*), intent(in) :: param_name
    type(rparam_type) :: param_group_create_rparam

    param_group_create_rparam = this%rparams%create_param(param_name)
  end function

  function param_group_create_lparam(this, param_name)
    class(param_group_type), intent(inout) :: this
    character(len=*), intent(in) :: param_name
    type(lparam_type) :: param_group_create_lparam

    param_group_create_lparam = this%lparams%create_param(param_name)
  end function

  function param_group_get_ivalue(this, param_name)
    class(param_group_type), intent(in) :: this
    character(len=*), intent(in) :: param_name
    integer :: param_group_get_ivalue
    type(iroparam_type) :: toto
    toto = this%iparams%get_roparam(param_name)
    param_group_get_ivalue = toto%get_value()
  end function param_group_get_ivalue

  function param_group_get_iparam(this, param_name)
    class(param_group_type), intent(inout) :: this
    character(len=*), intent(in) :: param_name
    type(iparam_type) :: param_group_get_iparam
    param_group_get_iparam = this%iparams%get_param(param_name)
  end function param_group_get_iparam

  function param_group_get_iroparam(this, param_name)
    class(param_group_type), intent(in) :: this
    character(len=*), intent(in) :: param_name
    type(iroparam_type) :: param_group_get_iroparam
    param_group_get_iroparam = this%iparams%get_roparam(param_name)
  end function param_group_get_iroparam

  ! logical function basis_is_twomol(this)
  !   implicit none
  !   class(param_group_type), intent(in) :: this
  !   integer :: ibasty
  !   ibasty = this%id
  !   basis_is_twomol = is_twomol(ibasty)
  ! end function basis_is_twomol

end module mod_param_group

#define TESTING_PARAM_GOUP
#ifdef TESTING_PARAM_GOUP

subroutine test_read_only(iparam)
  use mod_param_array_type_i, only: iparam_type => param_type
  implicit none
  type(iparam_type), intent(in) :: iparam
  write (*,*) 'test_read_only : value of iparam = ', iparam%get_value()
  ! call iparam%set_value(7) ! this would not compile, which is good !
end subroutine test_read_only

subroutine test_param_group_read_only(param_group)
  use mod_param_group, only: param_group_type, iroparam_type
  implicit none
  type(param_group_type), intent(in) :: param_group
  type(iroparam_type) :: iparam
  write (*,*) 'test_param_group_read_only : testing that this function cannot modify the parameters'
  iparam = param_group%get_iroparam('JTOT2')
  ! call iparam%value(13)  ! compiler doesn't allow that, which is what we want ! Error: ‘set_value’ at (1) is not a member of the ‘roparam_type’ structure; did you mean ‘value’?
  ! iparam%value = 13  ! compiler doesn't allow that, which is what we want : Error: Component ‘value’ at (1) is a PRIVATE component of ‘roparam_type’
  write (*,*) 'test_param_group_read_only : value of VVMIN = ', param_group%get_ivalue('VVMIN')  

end subroutine test_param_group_read_only
  
program test
  use mod_param_group, only: param_group_type, iparam_type
  integer, target :: toto
  integer, pointer :: jtot2
  integer :: jtot2_index
  type(iparam_type) :: vmin_par

  type(param_group_type), target :: group

  write (*,*) 'populating a group of integer parameters'
  call group%iparams%append('JTOT1')
  call group%iparams%append('JTOT2')
  call group%iparams%append('JMIN')
  call group%iparams%append('JMAX')
  vmin_par = group%create_iparam('VMIN')
  call vmin_par%set_value(66)

  call group%iparams%set('JTOT2', 42)
  write (*,*) 'value of JTOT2 = ', group%iparams%get('JTOT2')

  write (*,*) 'current contents of the integer parameters'
  call group%iparams%print()

  write (*,*) 'setting VMIN to 666'
  call vmin_par%set_value(666)

  write (*,*) 'current contents of the integer parameters'
  call group%iparams%print()

  write (*,*) 'renaming VMIN as VVMIN'
  call vmin_par%set_name('VVMIN')

  write (*,*) 'current contents of the integer parameters'
  call group%iparams%print()

  call test_read_only(vmin_par)

  write (*,*) 'current contents of the integer parameters'
  call group%iparams%print()

  call test_param_group_read_only(group)

  write (*,*) 'current contents of the integer parameters'
  call group%iparams%print()

  jtot2 => toto
  jtot2 => group%iparams%values(2)

  write (*,*) 'setting JTOT2 to 47'
  jtot2 = group%iparams%get('JTOT2')
  jtot2 = 47
  write (*,*) 'current contents of the integer parameters'
  call group%iparams%print()

  call group%rparams%append('PI')
  call group%rparams%set('PI', 3.14d0)
  write (*,*) group%rparams%get('PI')
  write (*,*) 'current contents of the real parameters'
  call group%rparams%print()

  call group%lparams%append('IS_TRUE')
  call group%lparams%set('IS_TRUE', .true.)
  write (*,*) group%lparams%get('IS_TRUE')

  write (*,*) 'current contents of the logical parameters'
  call group%lparams%print()

end program

#endif