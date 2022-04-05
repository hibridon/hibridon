#define PASTE(X)         X
#define PASTE2(X)        PASTE(X)_
#define CONCATHELP(X, Y) PASTE2(X)Y
#define CONCAT(X, Y)     CONCATHELP(X,Y)

!#define PARAM_ARRAY_TYPE      CONCAT(param_array_type,T)
#define MODULE_NAME           CONCAT(mod_param_array_type,T)
#define MAX_PARAMS 10
module MODULE_NAME

  type :: param_array_type
    !integer, len :: max_params = 10
    TT, dimension(MAX_PARAMS) :: values
    character(len=8), dimension(MAX_PARAMS) :: names
    integer :: num_params = 0
  contains
    procedure, public :: append => param_array_append
    procedure, public :: set    => param_array_set
    procedure, public :: get    => param_array_get
    procedure, public :: get_ptr    => param_array_get_ptr
    procedure, public :: print  => param_array_print
  end type

contains
  subroutine param_array_append(this, param_name)
    class(param_array_type), intent(inout) :: this
    character(len=*), intent(in) :: param_name
    ! write (*,*) 'param_array_append : appending ', param_name
    this%num_params = this%num_params + 1
    this%names( this%num_params ) = param_name
  end subroutine

  subroutine param_array_set(this, param_name, param_value)
    class(param_array_type), intent(inout) :: this
    character(len=*), intent(in) :: param_name
    TT, intent(in) :: param_value

    integer :: i
    ! write (*,*) 'param_array_set : ', param_name, ' <- ', param_value
    do i = 1, this%num_params
      ! write (*,*) 'param_array_set : ', this%names(i)
      if ( this%names(i) == param_name ) then
        this%values(i) = param_value
        return
      end if

    end do
    ASSERT( .false. )
  end subroutine

  function param_array_get(this, param_name)
    class(param_array_type), intent(in) :: this
    character(len=*), intent(in) :: param_name
    TT :: param_array_get
    integer :: i
    do i = 1, this%num_params
      if ( this%names(i) == param_name ) then
        param_array_get = this%values(i)
        return
      end if
    end do
    ASSERT( .false. )
  end function

  function param_array_get_ptr(this, param_name)
    class(param_array_type), intent(in) :: this
    character(len=*), intent(in) :: param_name
    TT, pointer :: param_array_get
    integer :: i
    do i = 1, this%num_params
      if ( this%names(i) == param_name ) then
        param_array_get = this%values(i)
        return
      end if
    end do
    ASSERT( .false. )
  end function

  subroutine param_array_print(this)
    class(param_array_type), intent(in) :: this
    integer :: i
    do i = 1, this%num_params
      write (*, *) this%names(i), this%values(i)
    end do
  end subroutine

end module MODULE_NAME

#undef T
#undef TT


#undef PASTE
#undef PASTE2
#undef CONCATHELP
#undef CONCAT
! #undef PARAM_ARRAY_TYPE
#undef MODULE_NAME
