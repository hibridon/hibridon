#include "assert.h"
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

#endif
program test
  use mod_param_group, only: param_group_type, iparam_type
  integer, target :: toto
  integer, pointer :: jtot2
  integer :: jtot2_index
  integer, pointer :: int_ptr
  integer, pointer :: int_ptr2
  type(iparam_type) :: vmin_par
  character(len=64) :: param_value_as_string

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

  ! testing get_vref
  toto = 777
  write (*,*) 'toto : ', toto
  int_ptr => toto
  write (*,*) 'toto : ', toto
  write (*,*) 'int_ptr : ', int_ptr
  int_ptr = 888
  write (*,*) 'toto : ', toto
  write (*,*) 'int_ptr : ', int_ptr

  

  toto = vmin_par%get_vref()
  write (*,*) 'toto : ', toto
  int_ptr = vmin_par%get_vref()
  write (*,*) 'int_ptr : ', int_ptr
  vmin_par%get_vref() = 43
  write (*,'("vmin_par value:",i5)') vmin_par%get_value()
  ASSERT(vmin_par%get_value() == 43)

  param_value_as_string = '42'
  read (param_value_as_string, *) vmin_par%get_vref()
  write (*,'("vmin_par value:",i5)') vmin_par%get_value()
  ASSERT(vmin_par%get_value() == 42)

end program