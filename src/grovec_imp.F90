  subroutine CONCATENATE(GROVEC_CLASS_NAME,_block_destroy)(this)
    type(CONCATENATE(GROVEC_CLASS_NAME,_block))             :: this
    if (allocated(this%p)) then
      ! write (6,*) 'deallocating grovec_block%p'
      deallocate(this%p)
    end if
  end subroutine

  function CONCATENATE(create_,GROVEC_CLASS_NAME)(block_size, num_blocks) Result(g)
    integer, intent(in) :: block_size
    integer, intent(in) :: num_blocks
    type(GROVEC_CLASS_NAME) :: g
    g%num_elements = 0
    g%num_blocks = num_blocks
    g%num_allocated_blocks = 0
    g%block_size = block_size
    allocate(g%blocks(0:g%num_blocks-1))
  end function CONCATENATE(create_,GROVEC_CLASS_NAME)

  subroutine CONCATENATE(GROVEC_CLASS_NAME,_destroy)(this)
    type(GROVEC_CLASS_NAME)             :: this
    if (allocated(this%blocks)) then
      ! write (6,*) 'deallocating grovec%blocks'
      deallocate(this%blocks)
    end if
  end subroutine

  subroutine CONCATENATE(GROVEC_CLASS_NAME,_ensure_block_is_available)(this, block_index)
    class(GROVEC_CLASS_NAME)             :: this
    integer, intent(in)        :: block_index

    if (block_index >= this%num_blocks) then
      stop 'block_index exceeds the number of blocks'
    end if
    do while (this%num_allocated_blocks <= block_index)
      allocate(this%blocks(this%num_allocated_blocks)%p(0:this%block_size-1))
      this%blocks(this%num_allocated_blocks)%block_size = this%block_size
      this%num_allocated_blocks = this%num_allocated_blocks + 1
    end do
  end subroutine

  subroutine CONCATENATE(GROVEC_CLASS_NAME,_set_element)(this, el_index, el_value)
    class(GROVEC_CLASS_NAME)             :: this
    integer(8), intent(in)      :: el_index
    GROVEC_ELEMENT_TYPE, intent(in)        :: el_value

    integer                    :: block_index
    integer                    :: el_index_in_block

    block_index =  (el_index-1) / this%block_size
    el_index_in_block = el_index - 1 - (block_index * this%block_size)
    
    call CONCATENATE(GROVEC_CLASS_NAME,_ensure_block_is_available)(this, block_index)

    this%blocks(block_index)%p(el_index_in_block) = el_value
    this%num_elements = max(this%num_elements, el_index)
  end subroutine

  function CONCATENATE(GROVEC_CLASS_NAME,_get_element)(this, el_index)
    class(GROVEC_CLASS_NAME)             :: this
    integer(8), intent(in)        :: el_index
    GROVEC_ELEMENT_TYPE                    :: CONCATENATE(GROVEC_CLASS_NAME,_get_element)

    integer                    :: block_index
    integer                    :: el_index_in_block

    block_index =  (el_index-1) / this%block_size
    el_index_in_block = el_index - 1 - (block_index * this%block_size)
    CONCATENATE(GROVEC_CLASS_NAME,_get_element) = this%blocks(block_index)%p(el_index_in_block)
  end function

  subroutine CONCATENATE(GROVEC_CLASS_NAME,_append)(this, el_value)
    class(GROVEC_CLASS_NAME)             :: this
    GROVEC_ELEMENT_TYPE, intent(in)        :: el_value

    call this%set_element(this%num_elements+1, el_value)
  end subroutine CONCATENATE(GROVEC_CLASS_NAME,_append)

  subroutine CONCATENATE(GROVEC_CLASS_NAME,_empty)(this)
    class(GROVEC_CLASS_NAME)             :: this
    integer :: block_index = 0
    do block_index = 0, this%num_allocated_blocks-1
      deallocate(this%blocks(block_index)%p)
    end do
    this%num_elements = 0
    this%num_allocated_blocks = 0
  end subroutine CONCATENATE(GROVEC_CLASS_NAME,_empty)
