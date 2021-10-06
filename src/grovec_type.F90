  type                         :: CONCATENATE(GROVEC_CLASS_NAME,_block)
    integer                    :: block_size = 0
    GROVEC_ELEMENT_TYPE, allocatable           :: p(:)
  contains
    final                      :: CONCATENATE(GROVEC_CLASS_NAME,_block_destroy)
  end type CONCATENATE(GROVEC_CLASS_NAME,_block)
  type, public                 :: GROVEC_CLASS_NAME
    integer*8                  :: num_elements = 0
    integer                    :: num_blocks = 0
    integer                    :: num_allocated_blocks = 0
    integer                    :: block_size = 0
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

