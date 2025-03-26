! interfaces to blas and lapack subroutines
module mod_hiblas
  interface

    function ddot(n, dx, incx, dy, incy)
      integer, intent(in) :: n
      real(8), intent(in) :: dx(*)
      integer, intent(in) :: incx
      real(8), intent(in) :: dy(*)
      integer, intent(in) :: incy
      real(8) :: ddot
    end function

  end interface
end module mod_hiblas
