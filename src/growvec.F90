module mod_growvec
  implicit none
  private
  type                         :: block
    integer                    :: block_size = 1000000
    real(8), pointer           :: p(:) => null()
  end type block
  type, public                 :: growvec
    integer(8)                 :: num_elements = 0
    integer                    :: num_blocks = 1000
    integer                    :: num_allocated_blocks = 0
    integer                    :: block_size = 1000000
    type(block), pointer       :: blocks(:) => null()
  contains
    procedure                  :: set_element => growvec_set_element
    procedure                  :: get_element => growvec_get_element

  end type growvec

  interface growvec
    module procedure create_growvec
  end interface growvec


contains
  function create_growvec() Result(g)
    !integer(8) :: block_size    
    type(growvec) :: g
    write (6,*) 'create_growvec'
    g%num_elements = 0
    g%num_blocks = 1024
    g%num_allocated_blocks = 0
    g%block_size = 1024*1024
    allocate(g%blocks(0:g%num_blocks-1))
  end function create_growvec

  subroutine growvec_ensure_block_is_available(this, block_index)
    class(growvec)             :: this
    integer, intent(in)        :: block_index

    write (6,*) 'growvec_ensure_block_is_available: block_index =', block_index
    do while (this%num_allocated_blocks <= block_index)
      write (6,*) 'allocating block', this%num_allocated_blocks
      write (6,*) 'this%block_size = ', this%block_size
      allocate(this%blocks(this%num_allocated_blocks)%p(0:this%block_size-1))
      write (6,*) 'block allocated'
      this%num_allocated_blocks = this%num_allocated_blocks + 1
      write (6,*) 'coucou'
    end do
  end subroutine

  subroutine growvec_set_element(this, el_index, el_value)
    class(growvec)             :: this
    integer, intent(in)        :: el_index
    real(8), intent(in)        :: el_value

    integer                    :: block_index
    integer                    :: el_index_in_block

    write (6,*) 'el_index=', el_index
    write (6,*) 'this%block_size=', this%block_size
    block_index =  (el_index-1) / this%block_size
    write (6,*) 'block_index=', block_index
    el_index_in_block = el_index - (block_index * this%block_size)
    write (6,*) 'el_index_in_block=', el_index_in_block
    
    call growvec_ensure_block_is_available(this, block_index)

    this%blocks(block_index)%p(el_index_in_block) = el_value
  end subroutine

  function growvec_get_element(this, el_index)
    class(growvec)             :: this
    integer, intent(in)        :: el_index
    real(8)                    :: growvec_get_element

    integer                    :: block_index
    integer                    :: el_index_in_block

    block_index =  (el_index-1) / this%block_size
    !write (6,*) 'growvec_get_element: block_index=', block_index
    el_index_in_block = el_index - (block_index * this%block_size)
    !write (6,*) 'growvec_get_element: el_index_in_block=', el_index_in_block
    growvec_get_element = this%blocks(block_index)%p(el_index_in_block)
  end function

end module mod_growvec

program test_growvec
use mod_growvec, only:growvec
type(growvec) :: g1
g1 = growvec() ! (block_size=1024*1024)
!g1 = create_growvec(block_size=1024*1024)
call g1%set_element(1, 3.d0)
call g1%set_element(1024, 1024.d0)
call g1%set_element(1024*1024+1, 1025.d0)
write (6,*) 'g1(1)=', g1%get_element(1)
write (6,*) 'g1(1024)=', g1%get_element(1024)
write (6,*) 'g1(1024*1024+1)=', g1%get_element(1024*1024+1)
end program test_growvec