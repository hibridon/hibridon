#define PASTE(X)         X
#define PASTE2(X)        PASTE(X)_
#define CONCATHELP(X, Y) PASTE2(X)Y
#define CONCAT(X, Y)     CONCATHELP(X,Y)

!#define PARAM_ARRAY_TYPE      CONCAT(param_array_type,T)
#define MODULE_NAME           CONCAT(mod_param_array_type,T)
#define MAX_PARAMS 10
module MODULE_NAME

  type :: param_type
  private
    character(len=8), pointer :: name
    TT, pointer :: value
  contains
  procedure, public :: get_value => param_get_value
  procedure, public :: set_value => param_set_value
  procedure, public :: get_name => param_get_name
  procedure, public :: set_name => param_set_name
  end type

  type :: roparam_type
  private
    character(len=8), pointer :: name
    TT, pointer :: value
  contains
  procedure, public :: get_value => roparam_get_value
  procedure, public :: get_name => roparam_get_name
  end type

  type :: param_array_type
    !integer, len :: max_params = 10
    TT, dimension(MAX_PARAMS) :: values
    character(len=8), dimension(MAX_PARAMS) :: names
    integer :: num_params = 0
  contains
    procedure, public :: append => param_array_append
    procedure, public :: create => param_array_create
    procedure, public :: create_param => param_array_create_param
    procedure, public :: set    => param_array_set
    procedure, public :: get    => param_array_get
    procedure, public :: get_ptr    => param_array_get_ptr
    procedure, public :: get_param    => param_array_get_param
    procedure, public :: get_roparam    => param_array_get_roparam
    procedure, public :: print  => param_array_print
  end type

contains
  function param_get_value(this)
    class(param_type), intent(in) :: this
    TT :: param_get_value
    param_get_value = this%value
  end function param_get_value

  subroutine param_set_value(this, value)
    class(param_type), intent(inout) :: this
    TT, intent(in) :: value
    this%value = value
  end subroutine param_set_value

  function param_get_name(this)
    class(param_type), intent(in) :: this
    character(len=8) :: param_get_name
    param_get_name = this%name
  end function param_get_name

  subroutine param_set_name(this, name)
    class(param_type), intent(inout) :: this
    character(len=*), intent(in) :: name
    this%name = name
  end subroutine param_set_name

  function roparam_get_value(this)
    class(roparam_type), intent(in) :: this
    TT :: roparam_get_value
    roparam_get_value = this%value
  end function roparam_get_value

  function roparam_get_name(this)
    class(roparam_type), intent(in) :: this
    character(len=8) :: roparam_get_name
    roparam_get_name = this%name
  end function roparam_get_name


  subroutine param_array_append(this, param_name)
    class(param_array_type), intent(inout) :: this
    character(len=*), intent(in) :: param_name
    ! write (*,*) 'param_array_append : appending ', param_name
    if(this%num_params == MAX_PARAMS) stop 'not enough room for storing a new parameter: please increase MAX_PARAMS'
    this%num_params = this%num_params + 1
    this%names( this%num_params ) = param_name
  end subroutine

  function param_array_create(this, param_name)
    class(param_array_type), intent(inout) :: this
    character(len=*), intent(in) :: param_name
    TT, pointer :: param_array_create
    ! write (*,*) 'param_array_append : appending ', param_name
    if (this%num_params == MAX_PARAMS) stop 'not enough room for storing a new parameter: please increase MAX_PARAMS'
    this%num_params = this%num_params + 1
    this%names( this%num_params ) = param_name
    param_array_create = this%values( this%num_params )
  end function

  function param_array_create_param(this, param_name)
    class(param_array_type), intent(inout), target :: this
    character(len=*), intent(in) :: param_name
    type(param_type) :: param_array_create_param
    ! write (*,*) 'param_array_append : appending ', param_name
    this%num_params = this%num_params + 1
    this%names( this%num_params ) = param_name
    param_array_create_param%name => this%names( this%num_params )
    param_array_create_param%value => this%values( this%num_params )
  end function

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
    TT, pointer :: param_array_get_ptr
    integer :: i
    do i = 1, this%num_params
      if ( this%names(i) == param_name ) then
        param_array_get_ptr = this%values(i)
        return
      end if
    end do
    ASSERT( .false. )
  end function

  function param_array_get_param(this, param_name)
    class(param_array_type), intent(in), target :: this
    character(len=*), intent(in) :: param_name
    type(param_type) :: param_array_get_param
    integer :: i
    do i = 1, this%num_params
      if ( this%names(i) == param_name ) then
        param_array_get_param%value => this%values(i)
        param_array_get_param%name => this%names(i)
        return
      end if
    end do
    ASSERT( .false. )
  end function

  function param_array_get_roparam(this, param_name)
    class(param_array_type), intent(in), target :: this
    character(len=*), intent(in) :: param_name
    type(roparam_type) :: param_array_get_roparam
    integer :: i
    do i = 1, this%num_params
      if ( this%names(i) == param_name ) then
        param_array_get_roparam%value => this%values(i)
        param_array_get_roparam%name => this%names(i)
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
