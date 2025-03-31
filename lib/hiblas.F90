! interfaces to blas and lapack subroutines
module mod_hiblas
  interface

    subroutine dcopy(n, dx, incx, dy, incy)
      integer, intent(in) :: n
      real(8), intent(in) :: dx(*)
      integer, intent(in) :: incx
      real(8), intent(out) :: dy(*)
      integer, intent(in) :: incy
    end subroutine

    function ddot(n, dx, incx, dy, incy)
      integer, intent(in) :: n
      real(8), intent(in) :: dx(*)
      integer, intent(in) :: incx
      real(8), intent(in) :: dy(*)
      integer, intent(in) :: incy
      real(8) :: ddot
    end function

    subroutine dscal(n, da, dx, incx)
      integer, intent(in) :: n
      real(8), intent(in) :: da
      real(8), intent(inout) :: dx(*)
      integer, intent(in) :: incx
    end subroutine

  end interface
end module mod_hiblas
