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

    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character, intent(in) :: transa
      character, intent(in) :: transb
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, intent(in) :: k
      real(8), intent(in) :: alpha
      real(8), dimension(lda,*), intent(in) :: a
      integer, intent(in) :: lda
      real(8), dimension(ldb,*), intent(in) :: b
      integer, intent(in) :: ldb
      real(8), intent(in) :: beta
      real(8), dimension(ldc,*), intent(inout) :: c
      integer, intent(in) :: ldc
    end subroutine

    subroutine dscal(n, da, dx, incx)
      integer, intent(in) :: n
      real(8), intent(in) :: da
      real(8), intent(inout) :: dx(*)
      integer, intent(in) :: incx
    end subroutine

    subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
      character, intent(in) :: jobz
      character, intent(in) :: uplo
      integer, intent(in) :: n
      real(8), intent(inout) :: a(lda, *)
      integer, intent(in) :: lda
      real(8), intent(out) :: w(*)
      real(8), intent(out) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine

    subroutine dsyevx(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
      character, intent(in) :: jobz
      character, intent(in) :: range
      character, intent(in) :: uplo
      integer, intent(in) :: n
      real(8), intent(inout) :: a(lda, *)
      integer, intent(in) :: lda
      real(8), intent(in) :: vl
      real(8), intent(in) :: vu
      integer, intent(in) :: il
      integer, intent(in) :: iu
      real(8), intent(in) :: abstol
      integer, intent(out) :: m
      real(8), intent(out) :: w(*)
      real(8), intent(out) :: z(ldz, *)
      integer, intent(in) :: ldz
      real(8), intent(out) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: iwork(*)
      integer, intent(out) :: ifail(*)
      integer, intent(out) :: info
    end subroutine

    subroutine dsygvx(itype, jobz, range, uplo, n, a, lda, b, ldb, vl, vu, il, iu, abstol, m, w, z, ldz, work, lwork, iwork, ifail, info)
      integer, intent(in) :: itype
      character, intent(in) :: jobz
      character, intent(in) :: range
      character, intent(in) :: uplo
      integer, intent(in) :: n
      real(8), intent(inout) :: a(lda, *)
      integer, intent(in) :: lda
      real(8), intent(inout) :: b(ldb, *)
      integer, intent(in) :: ldb
      real(8), intent(in) :: vl
      real(8), intent(in) :: vu
      integer, intent(in) :: il
      integer, intent(in) :: iu
      real(8), intent(in) :: abstol
      integer, intent(out) :: m
      real(8), intent(out) :: w(*)
      real(8), intent(out) :: z(ldz, *)
      integer, intent(in) :: ldz
      real(8), intent(out) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: iwork(*)
      integer, intent(out) :: ifail(*)
      integer, intent(out) :: info
    end subroutine

  end interface
end module mod_hiblas
