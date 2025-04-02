#include "assert.h"
! interfaces to blas and lapack subroutines
module mod_hiblas
  interface

    subroutine daxpy(n, da, dx, incx, dy, incy)
      integer, intent(in) :: n
      real(8), intent(in) :: da
      real(8), intent(in) :: dx(*)
      integer, intent(in) :: incx
      real(8), intent(inout) :: dy(*)
      integer, intent(in) :: incy
    end subroutine

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

    subroutine dgetrf(m, n, a, lda, ipiv, info)
      integer, intent(in) :: m
      integer, intent(in) :: n
      real(8), intent(inout) :: a(lda,*)
      integer, intent(in) :: lda
      integer, intent(out) :: ipiv(*)
      integer, intent(out) :: info
    end subroutine

    subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
      integer, intent(in) :: n
      real(8), intent(inout) :: a(lda,*)
      integer, intent(in) :: lda
      integer, intent(in) :: ipiv(*)
      real(8), intent(out) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine

    subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
      character, intent(in) :: trans
      integer, intent(in) :: n
      integer, intent(in) :: nrhs
      real(8), intent(in) :: a(lda,*)
      integer, intent(in) :: lda
      integer, intent(in) :: ipiv(*)
      real(8), intent(inout) :: b(ldb,*)
      integer, intent(in) :: ldb
      integer, intent(out) :: info
    end subroutine

    subroutine drot(n, dx, incx, dy, incy, c, s)
      integer, intent(in) :: n
      real(8), intent(inout) :: dx(*)
      integer, intent(in) :: incx
      real(8), intent(inout) :: dy(*)
      integer, intent(in) :: incy
      real(8), intent(in) :: c
      real(8), intent(in) :: s
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

contains

subroutine daxpy_wrapper(n, da, dx, incx, dy, incy)
implicit none
integer, intent(in)                       :: n
real(8), intent(in)                       :: da
real(8), intent(in), dimension(1+(n+1)*incx)    :: dx
integer, intent(in)                       :: incx
real(8), intent(inout), dimension(1+(n+1)*incx)    :: dy
integer, intent(in)                       :: incy

real(8) :: sum_for_nan_detection
integer :: i

sum_for_nan_detection = dx(1)
do i=1,n
      sum_for_nan_detection = sum_for_nan_detection + dx(1+((i-1)*incx))
end do
sum_for_nan_detection = dx(1)
sum_for_nan_detection = sum(dx(1:n:incx))

call daxpy(n, da, dx, incx, dy, incy)
end subroutine
!     ------------------------------------------------------------------
subroutine dgetri_wrapper(n, a, lda, ipiv, work, lwork, info)
implicit none
integer, intent(in)                       :: n
real(8), intent(inout), dimension(lda, n) :: a
integer, intent(in)                       :: lda
integer, intent(in),    dimension(n)      :: ipiv
real(8), intent(out),   dimension(lwork)  :: work
integer, intent(in)                       :: lwork
integer, intent(out)                      :: info

real(8) :: sum_for_nan_detection

#if defined(FLOATING_POINT_ERRORS_CAUSE_EXCEPTIONS)
! for some reason, dgetri triggers a floating point exception 
! if work array contains NaNs
work = 0.0
#endif       
info = 0
sum_for_nan_detection = sum(a(1:n,:))
sum_for_nan_detection = sum(ipiv)

call dgetri(n, a, lda, ipiv, work, lwork, info)
end subroutine
!     ------------------------------------------------------------------
subroutine dsyevr_wrapper(jobz, range, uplo, n, a, lda, vl, vu, &
      il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)
use mod_assert, only: fassert
implicit none
#ifdef FLOATING_POINT_ERRORS_CAUSE_EXCEPTIONS
! use dsyevx instead of dsyevr if floating point exceptions are on
! because xxxevr routines may create NaNs according to 
! http://www.democritos.it/activities/IT-MC/ACML-4.0.1/LAPACK_002dIEEE.html
! Also, we found that dsyevr triggers floating point exceptions
! in its call to ilaenv :
! http://www.netlib.org/lapack/explore-html/df/d3b/dsyevr_8f_source.html
! ieeeok = ilaenv( 10, 'DSYEVR', 'N', 1, 2, 3, 4 )

! http://www.democritos.it/activities/IT-MC/ACML-4.0.1/LAPACK_002dIEEE.html
!CALL ilaenvset(10,'X','X',0,0,0,0,1,INFO)
!CALL ilaenvset(11,'X','X',0,0,0,0,1,INFO)

! so, to disable these false positives, we replace dsyevr with dsyevx, which is slower but does the same thing.
#define USE_DSYEVX
#endif
character, intent(in) :: jobz
character, intent(in) :: range
character, intent(in) :: uplo
integer, intent(in) :: n
real(8), dimension(lda, n), intent(inout) :: a
integer, intent(in) :: lda
real(8), intent(in) :: vl
real(8), intent(in) :: vu
integer, intent(in) :: il
integer, intent(in) :: iu
real(8), intent(in) :: abstol
integer, intent(out) :: m
real(8), dimension(n), intent(out) :: w
real(8), dimension(ldz, n), intent(out) :: z
integer, intent(in) :: ldz
integer, dimension(2*n), intent(out) :: isuppz
real(8), dimension(lwork), intent(out) :: work
integer, intent(in) :: lwork
integer, dimension(liwork), intent(out) :: iwork
integer, intent(in) :: liwork
integer, intent(out) :: info

#ifdef USE_DSYEVX
integer :: dsyevx_lwork
real(8), dimension(8*n) :: dsyevx_work
integer, dimension(5*n) :: dsyevx_iwork
integer, dimension(n) :: ifail
#endif
! real(8), dimension(lda, n) :: eye
real(8) :: sum_for_nan_detection
! real(8) :: dummy
! integer :: i

! eye = a(1:n,:)
! eye = 0
! do i = 1, n
!       eye(i,i) = 1.0
! end do
! dummy = 1.0 / sum_for_nan_detection
! eye = 0.0
#ifdef DEBUG_DSYEV_WRAPPER
write(6, *) 'dsyevr_wrapper: jobz = ', jobz
write(6, *) 'dsyevr_wrapper: range = ', range
write(6, *) 'dsyevr_wrapper: uplo = ', uplo
write(6, *) 'dsyevr_wrapper: n = ', n
write(6, *) 'dsyevr_wrapper: lda = ', lda
write(6, *) 'dsyevr_wrapper: vl = ', vl
write(6, *) 'dsyevr_wrapper: vu = ', vu
write(6, *) 'dsyevr_wrapper: il = ', il
write(6, *) 'dsyevr_wrapper: iu = ', iu
write(6, *) 'dsyevr_wrapper: abstol = ', abstol
write(6, *) 'dsyevr_wrapper: ldz = ', ldz
write(6, *) 'dsyevr_wrapper: lwork = ', lwork
write(6, *) 'dsyevr_wrapper: liwork = ', liwork
write(6, *) 'dsyevr_wrapper: a(1,1) = ', a(1,1)
#endif
!sum_for_nan_detection = a(1,1) + a(1,2)
sum_for_nan_detection = sum(a(1:n,:))
! write(6, *) 'dsyevr_wrapper: eye = ', eye
! write(6, *) 'dsyevr_wrapper: a = ', a

! m = 0
! w = 0.0
! if ( jobz == 'V' ) then
!      z = 0.0
! end if
! isuppz = 0
! straing: it seems necessay to initialize work, otherwise
! floating point exceptions can happen (at least on gfortran
! debug build with fpe on)
work = 0 
!iwork = 0
!info = 0


#ifdef USE_DSYEVX
dsyevx_lwork = size(dsyevx_work)
dsyevx_work = 0.0
dsyevx_iwork = 0
isuppz = 0 ! to remove warning #6843 because of unused out argument
iwork = 0 ! to remove warning #6843 because of unused out argument
#ifdef DEBUG_DSYEV_WRAPPER
write(6, *) 'dsyevr_wrapper: dsyevx_lwork = ', dsyevx_lwork
write(6, *) 'dsyevr_wrapper: size(dsyevx_iwork) = ', size(dsyevx_iwork)
#endif
ifail = 0
call dsyevx(jobz, range, uplo, n, a, lda, vl, vu, &
      il, iu, abstol, m, w, z, ldz, dsyevx_work, dsyevx_lwork, dsyevx_iwork, ifail, info)
#else
call dsyevr(jobz, range, uplo, n, a, lda, vl, vu, &
      il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)
#endif
end subroutine

!#define DEBUG_DSYEV_WRAPPER
subroutine dsyev_wrapper(jobz, uplo, n, a, lda, w, work, lwork, info)
use mod_assert, only: fassert
implicit none
character, intent(in) :: jobz
character, intent(in) :: uplo
integer, intent(in) :: n
real(8), dimension(lda, n), intent(inout) :: a
integer, intent(in) :: lda
real(8), dimension(n), intent(out) :: w
real(8), dimension(lwork), intent(out) :: work
integer, intent(in) :: lwork
integer, intent(out) :: info

! real(8), dimension(lda, n) :: eye
real(8) :: sum_for_nan_detection
! real(8) :: dummy
! integer :: i

! eye = a(1:n,:)
! eye = 0
! do i = 1, n
!       eye(i,i) = 1.0
! end do
! dummy = 1.0 / sum_for_nan_detection
! eye = 0.0
#ifdef DEBUG_DSYEV_WRAPPER
write(6, *) 'dsyev_wrapper: jobz = ', jobz
write(6, *) 'dsyev_wrapper: uplo = ', uplo
write(6, *) 'dsyev_wrapper: n = ', n
write(6, *) 'dsyev_wrapper: lda = ', lda
write(6, *) 'dsyev_wrapper: lwork = ', lwork
write(6, *) 'dsyev_wrapper: a(1,1) = ', a(1,1)
#endif
!sum_for_nan_detection = a(1,1) + a(1,2)
sum_for_nan_detection = sum(a(1:n,:))
! write(6, *) 'dsyevr_wrapper: eye = ', eye
! write(6, *) 'dsyevr_wrapper: a = ', a

! m = 0
! w = 0.0
! if ( jobz == 'V' ) then
!      z = 0.0
! end if
! isuppz = 0
! straing: it seems necessay to initialize work, otherwise
! floating point exceptions can happen (at least on gfortran
! debug build with fpe on)
work = 0
!iwork = 0
!info = 0

call dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
if (info /= 0) then
  stop 43
end if
ASSERT(info == 0)
end subroutine

end module mod_hiblas
