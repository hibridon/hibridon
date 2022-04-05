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
  use mod_param_array_type_i, only: param_array_type_i => param_array_type
  use mod_param_array_type_r, only: param_array_type_r => param_array_type
  use mod_param_array_type_l, only: param_array_type_l => param_array_type
  implicit none
  integer, parameter :: kmaxpar = 150

  type :: param_group_type
    integer :: toto
    type(param_array_type_i) :: iparams
    type(param_array_type_r) :: rparams
    type(param_array_type_l) :: lparams
  !type(param_array_type_i(kmaxpar)) :: iparams
    !type(param_array_type_r(kmaxpar)) :: rparams
    !type(param_array_type_l(kmaxpar)) :: lparams

  contains 
    ! procedure, public :: is_twomol => basis_is_twomol
  end type param_group_type

  type(param_group_type) :: basis_params
 
contains

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

program test
  use mod_param_group, only: param_group_type
  integer, target :: toto
  integer, pointer :: jtot2
  integer :: jtot2_index

  type(param_group_type), target :: group
  call group%iparams%append('JTOT1')
  call group%iparams%append('JTOT2')
  call group%iparams%append('JMIN')
  call group%iparams%append('JMAX')
  call group%iparams%set('JTOT2', 42)
  write (*,*) group%iparams%get('JTOT2')
  call group%iparams%print()

  jtot2 => toto
  jtot2 => group%toto
  jtot2 => group%iparams%values(2)
  jtot2 = group%iparams%get('JTOT2')
  jtot2 = 47
  call group%iparams%print()

  call group%rparams%append('PI')
  call group%rparams%set('PI', 3.14d0)
  write (*,*) group%rparams%get('PI')
  call group%rparams%print()

  call group%lparams%append('IS_TRUE')
  call group%lparams%set('IS_TRUE', .true.)
  write (*,*) group%lparams%get('IS_TRUE')
  call group%lparams%print()

end program

#endif