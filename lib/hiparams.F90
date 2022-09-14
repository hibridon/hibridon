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