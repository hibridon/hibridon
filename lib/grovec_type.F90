! grovec is a growable vector : it automatically grows when elements are appended
! the type of elements contained in grovec can be chosen by the user: real, integer, or logical
! the grovec code uses a template mechanism to allow the support of any element type. As a result, this
! file is not expected to be compiled as is; instead it should be included from a proper fortran source file,
! after having defined the following preprocessing variables:
! - GROVEC_CLASS_NAME : the name chosen by the user for the resulting grovec type. A good name for a type representing a growable array of integers could be igrovec_type ('i' for integer)
! - GROVEC_ELEMENT_TYPE : the type of the elements stored in the resulting grovec type. Eg 'real(8)'

  ! <GROVEC_CLASS_NAME>_block : a block of elements 
  type, public                 :: CONCATENATE(GROVEC_CLASS_NAME,_block)
    integer                    :: block_size = 0   ! number of cells in this block (each cell can store an element of the growable array)
    GROVEC_ELEMENT_TYPE, allocatable           :: p(:) ! the actual array that is used as a storage for the elements of the growable array
  contains
    final                      :: CONCATENATE(GROVEC_CLASS_NAME,_block_destroy)  ! destructor of the bock
  end type CONCATENATE(GROVEC_CLASS_NAME,_block)

  ! <GROVEC_CLASS_NAME> : represents a growable array
  ! each instance of a growable array stores its elements in blocks (all blocks are expected to have the same size)
  type, public                 :: GROVEC_CLASS_NAME
    integer                    :: num_elements = 0   ! the current number of elements stored in this growable array
    integer                    :: num_blocks = 0     ! the number 
    integer                    :: num_allocated_blocks = 0  ! the current number of allocated blocks
    integer                    :: block_size = 0     ! the size of the blocks used by this growable array (in other words, the maximum number of elements contained in a block)
    type(CONCATENATE(GROVEC_CLASS_NAME,_block)), allocatable       :: blocks(:)
  contains
    final                      :: CONCATENATE(GROVEC_CLASS_NAME,_destroy)
    procedure                  :: set_element => CONCATENATE(GROVEC_CLASS_NAME,_set_element)
    procedure                  :: get_element => CONCATENATE(GROVEC_CLASS_NAME,_get_element)
    procedure                  :: append => CONCATENATE(GROVEC_CLASS_NAME,_append)
    procedure                  :: empty => CONCATENATE(GROVEC_CLASS_NAME,_empty)

  end type GROVEC_CLASS_NAME

  ! interface to create an instance of GROVEC_CLASS_NAME using a construct familiar to other languages:
  ! g1 = GROVEC_CLASS_NAME()
  interface GROVEC_CLASS_NAME
    module procedure CONCATENATE(create_,GROVEC_CLASS_NAME)
  end interface GROVEC_CLASS_NAME

