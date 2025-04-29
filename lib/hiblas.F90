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

    subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
      character, intent(in) :: trans
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, intent(in) :: nrhs
      real(8), intent(inout) :: a(lda, n)
      integer, intent(in) :: lda
      real(8), intent(inout) :: b(ldb, nrhs)
      integer, intent(in) :: ldb
      real(8), intent(out) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: info
    end subroutine

    subroutine dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, intent(in) :: nrhs
      real(8), intent(inout) :: a(lda, n)
      integer, intent(in) :: lda
      real(8), intent(inout) :: b(ldb, nrhs)
      integer, intent(in) :: ldb
      real(8), intent(out) :: s(*)
      real(8), intent(in) :: rcond
      integer, intent(out) :: rank
      real(8), intent(out) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: iwork(*)
      integer, intent(out) :: info
    end subroutine

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

    function dnrm2(n, x, incx)
      real(8) :: dnrm2
      integer, intent(in) :: n
      real(8), intent(in) :: x(*)
      integer, intent(in) :: incx
    end function

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

    subroutine dswap(n, dx, incx, dy, incy)
      integer, intent(in) :: n
      real(8), intent(inout) :: dx(*)
      integer, intent(in) :: incx
      real(8), intent(inout) :: dy(*)
      integer, intent(in) :: incy
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

    subroutine dsyevr(jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, m, w, z, ldz, isuppz, work, lwork, iwork, liwork, info)
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
      integer, intent(out) :: isuppz(*)
      real(8), intent(out) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: iwork(*)
      integer, intent(in) :: liwork
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

    subroutine dsygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info)
      integer, intent(in) :: itype
      character, intent(in) :: jobz
      character, intent(in) :: uplo
      integer, intent(in) :: n
      real(8), intent(inout) :: a(lda, *)
      integer, intent(in) :: lda
      real(8), intent(inout) :: b(ldb, *)
      integer, intent(in) :: ldb
      real(8), intent(out) :: w(*)
      real(8), intent(out) :: work(*)
      integer, intent(in) :: lwork
      integer, intent(out) :: iwork(*)
      integer, intent(in) :: liwork
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

    function idamax(n, dx, incx)
      integer :: idamax
      integer, intent(in) :: n
      real(8), intent(in) :: dx(*)
      integer, intent(in) :: incx
    end function


  end interface

contains

#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_IBM) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
integer          function ilaenv( ispec, name, opts, n1, n2, n3, &
                 n4 )
!
!  -- LAPACK auxiliary routine (preliminary version) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 20, 1992
!
!     .. Scalar Arguments ..
character*( * )    name, opts
integer            ispec, n1, n2, n3, n4
!     ..
!
!  Purpose
!  =======
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
logical            cname, sname
character*1        c1
character*2        c2, c4
character*3        c3
character*6        subnam
integer            i, ic, iz, nb, nbmin, nx
!     ..
!     .. Intrinsic Functions ..
intrinsic          char, ichar, int, min, real
!     ..
!     .. Executable Statements ..
!
go to ( 100, 100, 100, 400, 500, 600, 700, 800 ) ispec
!
!     Invalid value for ISPEC
!
ilaenv = -1
return
!
100 continue
!
!     Convert NAME to upper case if the first character is lower case.
!
ilaenv = 1
subnam = name
ic = ichar( subnam( 1:1 ) )
iz = ichar( 'Z' )
if( iz.eq.90 .or. iz.eq.122 ) then
!
!        ASCII character set
!
   if( ic.ge.97 .and. ic.le.122 ) then
      subnam( 1:1 ) = char( ic-32 )
      do 10 i = 2, 6
         ic = ichar( subnam( i:i ) )
         if( ic.ge.97 .and. ic.le.122 ) &
            subnam( i:i ) = char( ic-32 )
10       continue
   end if
!
else if( iz.eq.233 .or. iz.eq.169 ) then
!
!        EBCDIC character set
!
   if( ( ic.ge.129 .and. ic.le.137 ) .or. &
       ( ic.ge.145 .and. ic.le.153 ) .or. &
       ( ic.ge.162 .and. ic.le.169 ) ) then
      subnam( 1:1 ) = char( ic+64 )
      do 20 i = 2, 6
         ic = ichar( subnam( i:i ) )
         if( ( ic.ge.129 .and. ic.le.137 ) .or. &
             ( ic.ge.145 .and. ic.le.153 ) .or. &
             ( ic.ge.162 .and. ic.le.169 ) ) &
            subnam( i:i ) = char( ic+64 )
20       continue
   end if
!
else if( iz.eq.218 .or. iz.eq.250 ) then
!
!        Prime machines:  ASCII+128
!
   if( ic.ge.225 .and. ic.le.250 ) then
      subnam( 1:1 ) = char( ic-32 )
      do 30 i = 2, 6
         ic = ichar( subnam( i:i ) )
         if( ic.ge.225 .and. ic.le.250 ) &
            subnam( i:i ) = char( ic-32 )
30       continue
   end if
end if
!
c1 = subnam( 1:1 )
sname = c1.eq.'S' .or. c1.eq.'D'
cname = c1.eq.'C' .or. c1.eq.'Z'
if( .not.( cname .or. sname ) ) &
   return
c2 = subnam( 2:3 )
c3 = subnam( 4:6 )
c4 = c3( 2:3 )
!
go to ( 110, 200, 300 ) ispec
!
110 continue
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
nb = 1
!
if( c2.eq.'GE' ) then
   if( c3.eq.'TRF' ) then
      if( sname ) then
         nb = 64
      else
         nb = 64
      end if
   else if( c3.eq.'QRF' .or. c3.eq.'RQF' .or. c3.eq.'LQF' .or. &
            c3.eq.'QLF' ) then
      if( sname ) then
         nb = 32
      else
         nb = 32
      end if
   else if( c3.eq.'HRD' ) then
      if( sname ) then
         nb = 32
      else
         nb = 32
      end if
   else if( c3.eq.'BRD' ) then
      if( sname ) then
         nb = 32
      else
         nb = 32
      end if
   else if( c3.eq.'TRI' ) then
      if( sname ) then
         nb = 64
      else
         nb = 64
      end if
   end if
else if( c2.eq.'PO' ) then
   if( c3.eq.'TRF' ) then
      if( sname ) then
         nb = 64
      else
         nb = 64
      end if
   end if
else if( c2.eq.'SY' ) then
   if( c3.eq.'TRF' ) then
      if( sname ) then
         nb = 64
      else
         nb = 64
      end if
   else if( sname .and. c3.eq.'TRD' ) then
      nb = 1
   else if( sname .and. c3.eq.'GST' ) then
      nb = 64
   end if
else if( cname .and. c2.eq.'HE' ) then
   if( c3.eq.'TRF' ) then
      nb = 64
   else if( c3.eq.'TRD' ) then
      nb = 1
   else if( c3.eq.'GST' ) then
      nb = 64
   end if
else if( sname .and. c2.eq.'OR' ) then
   if( c3( 1:1 ).eq.'G' ) then
      if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or. &
          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or. &
          c4.eq.'BR' ) then
         nb = 32
      end if
   else if( c3( 1:1 ).eq.'M' ) then
      if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or. &
          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or. &
          c4.eq.'BR' ) then
         nb = 32
      end if
   end if
else if( cname .and. c2.eq.'UN' ) then
   if( c3( 1:1 ).eq.'G' ) then
      if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or. &
          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or. &
          c4.eq.'BR' ) then
         nb = 32
      end if
   else if( c3( 1:1 ).eq.'M' ) then
      if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or. &
          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or. &
          c4.eq.'BR' ) then
         nb = 32
      end if
   end if
else if( c2.eq.'GB' ) then
   if( c3.eq.'TRF' ) then
      if( sname ) then
         if( n4.le.64 ) then
            nb = 1
         else
            nb = 32
         end if
      else
         if( n4.le.64 ) then
            nb = 1
         else
            nb = 32
         end if
      end if
   end if
else if( c2.eq.'PB' ) then
   if( c3.eq.'TRF' ) then
      if( sname ) then
         if( n2.le.64 ) then
            nb = 1
         else
            nb = 32
         end if
      else
         if( n2.le.64 ) then
            nb = 1
         else
            nb = 32
         end if
      end if
   end if
else if( c2.eq.'TR' ) then
   if( c3.eq.'TRI' ) then
      if( sname ) then
         nb = 64
      else
         nb = 64
      end if
   end if
else if( c2.eq.'LA' ) then
   if( c3.eq.'UUM' ) then
      if( sname ) then
         nb = 64
      else
         nb = 64
      end if
   end if
else if( sname .and. c2.eq.'ST' ) then
   if( c3.eq.'EBZ' ) then
      nb = 1
   end if
end if
ilaenv = nb
return
!
200 continue
!
!     ISPEC = 2:  minimum block size
!
nbmin = 2
if( c2.eq.'GE' ) then
   if( c3.eq.'QRF' .or. c3.eq.'RQF' .or. c3.eq.'LQF' .or. &
       c3.eq.'QLF' ) then
      if( sname ) then
         nbmin = 2
      else
         nbmin = 2
      end if
   else if( c3.eq.'HRD' ) then
      if( sname ) then
         nbmin = 2
      else
         nbmin = 2
      end if
   else if( c3.eq.'BRD' ) then
      if( sname ) then
         nbmin = 2
      else
         nbmin = 2
      end if
   else if( c3.eq.'TRI' ) then
      if( sname ) then
         nbmin = 2
      else
         nbmin = 2
      end if
   end if
else if( c2.eq.'SY' ) then
   if( c3.eq.'TRF' ) then
      if( sname ) then
         nbmin = 2
      else
         nbmin = 2
      end if
   else if( sname .and. c3.eq.'TRD' ) then
      nbmin = 2
   end if
else if( cname .and. c2.eq.'HE' ) then
   if( c3.eq.'TRD' ) then
      nbmin = 2
   end if
else if( sname .and. c2.eq.'OR' ) then
   if( c3( 1:1 ).eq.'G' ) then
      if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or. &
          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or. &
          c4.eq.'BR' ) then
         nbmin = 2
      end if
   else if( c3( 1:1 ).eq.'M' ) then
      if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or. &
          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or. &
          c4.eq.'BR' ) then
         nbmin = 2
      end if
   end if
else if( cname .and. c2.eq.'UN' ) then
   if( c3( 1:1 ).eq.'G' ) then
      if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or. &
          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or. &
          c4.eq.'BR' ) then
         nbmin = 2
      end if
   else if( c3( 1:1 ).eq.'M' ) then
      if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or. &
          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or. &
          c4.eq.'BR' ) then
         nbmin = 2
      end if
   end if
end if
ilaenv = nbmin
return
!
300 continue
!
!     ISPEC = 3:  crossover point
!
nx = 0
if( c2.eq.'GE' ) then
   if( c3.eq.'QRF' .or. c3.eq.'RQF' .or. c3.eq.'LQF' .or. &
       c3.eq.'QLF' ) then
      if( sname ) then
         nx = 128
      else
         nx = 128
      end if
   else if( c3.eq.'HRD' ) then
      if( sname ) then
         nx = 128
      else
         nx = 128
      end if
   else if( c3.eq.'BRD' ) then
      if( sname ) then
         nx = 128
      else
         nx = 128
      end if
   end if
else if( c2.eq.'SY' ) then
   if( sname .and. c3.eq.'TRD' ) then
      nx = 1
   end if
else if( cname .and. c2.eq.'HE' ) then
   if( c3.eq.'TRD' ) then
      nx = 1
   end if
else if( sname .and. c2.eq.'OR' ) then
   if( c3( 1:1 ).eq.'G' ) then
      if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or. &
          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or. &
          c4.eq.'BR' ) then
         nx = 128
      end if
   end if
else if( cname .and. c2.eq.'UN' ) then
   if( c3( 1:1 ).eq.'G' ) then
      if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or. &
          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or. &
          c4.eq.'BR' ) then
         nx = 128
      end if
   end if
end if
ilaenv = nx
return
!
400 continue
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
ilaenv = 6
return
!
500 continue
!
!     ISPEC = 5:  minimum column dimension (not used)
!
ilaenv = 2
return
!
600 continue
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
ilaenv = int( real( min( n1, n2 ) )*1.6e0 )
return
!
700 continue
!
!     ISPEC = 7:  number of processors (not used)
!
ilaenv = 1
return
!
800 continue
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
ilaenv = 50
return
!
!     End of ILAENV
!
end
subroutine dlaev2( a, b, c, rt1, rt2, cs1, sn1 )
!
!  -- LAPACK auxiliary routine (version 1.0b) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
double precision   a, b, c, cs1, rt1, rt2, sn1
!     ..
!
!  Purpose
!  =======
!
!  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
!  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
!  eigenvector for RT1, giving the decomposition
!
!     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
!     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
!
!  Arguments
!  =========
!
!  A       (input) DOUBLE PRECISION
!          The (1,1) entry of the 2-by-2 matrix.
!
!  B       (input) DOUBLE PRECISION
!          The (1,2) entry and the conjugate of the (2,1) entry of the
!          2-by-2 matrix.
!
!  C       (input) DOUBLE PRECISION
!          The (2,2) entry of the 2-by-2 matrix.
!
!  RT1     (output) DOUBLE PRECISION
!          The eigenvalue of larger absolute value.
!
!  RT2     (output) DOUBLE PRECISION
!          The eigenvalue of smaller absolute value.
!
!  CS1     (output) DOUBLE PRECISION
!  SN1     (output) DOUBLE PRECISION
!          The vector (CS1, SN1) is a unit right eigenvector for RT1.
!
!  Further Details
!  ===============
!
!  RT1 is accurate to a few ulps barring over/underflow.
!
!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.
!
!  CS1 and SN1 are accurate to a few ulps barring over/underflow.
!
!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.
!
! =====================================================================
!
!     .. Parameters ..
double precision   one
parameter          ( one = 1.0d0 )
double precision   two
parameter          ( two = 2.0d0 )
double precision   zero
parameter          ( zero = 0.0d0 )
double precision   half
parameter          ( half = 0.5d0 )
!     ..
!     .. Local Scalars ..
integer            sgn1, sgn2
double precision   ab, acmn, acmx, acs, adf, cs, ct, df, rt, sm, &
                   tb, tn
!     ..
!     .. Intrinsic Functions ..
intrinsic          abs, sqrt
!     ..
!     .. Executable Statements ..
!
!     Compute the eigenvalues
!
sm = a + c
df = a - c
adf = abs( df )
tb = b + b
ab = abs( tb )
if( abs( a ).gt.abs( c ) ) then
   acmx = a
   acmn = c
else
   acmx = c
   acmn = a
end if
if( adf.gt.ab ) then
   rt = adf*sqrt( one+( ab / adf )**2 )
else if( adf.lt.ab ) then
   rt = ab*sqrt( one+( adf / ab )**2 )
else
!
!        Includes case AB=ADF=0
!
   rt = ab*sqrt( two )
end if
if( sm.lt.zero ) then
   rt1 = half*( sm-rt )
   sgn1 = -1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
   rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
else if( sm.gt.zero ) then
   rt1 = half*( sm+rt )
   sgn1 = 1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
   rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
else
!
!        Includes case RT1 = RT2 = 0
!
   rt1 = half*rt
   rt2 = -half*rt
   sgn1 = 1
end if
!
!     Compute the eigenvector
!
if( df.ge.zero ) then
   cs = df + rt
   sgn2 = 1
else
   cs = df - rt
   sgn2 = -1
end if
acs = abs( cs )
if( acs.gt.ab ) then
   ct = -tb / cs
   sn1 = one / sqrt( one+ct*ct )
   cs1 = ct*sn1
else
   if( ab.eq.zero ) then
      cs1 = one
      sn1 = zero
   else
      tn = -cs / tb
      cs1 = one / sqrt( one+tn*tn )
      sn1 = tn*cs1
   end if
end if
if( sgn1.eq.sgn2 ) then
   tn = cs1
   cs1 = -sn1
   sn1 = tn
end if
return
!
!     End of DLAEV2
!
end
subroutine dlasyf( uplo, n, nb, kb, a, lda, ipiv, w, ldw, info )
use mod_hiblas, only: dscal, dcopy, dgemm, dswap, idamax
!
!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
character          uplo
integer            info, kb, lda, ldw, n, nb
!     ..
!     .. Array Arguments ..
integer            ipiv( * )
double precision   a( lda, * ), w( ldw, * )
!     ..
!
!  Purpose
!  =======
!
!  DLASYF computes a partial factorization of a real symmetric matrix A
!  using the Bunch-Kaufman diagonal pivoting method. The partial
!  factorization has the form:
!
!  A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:
!        ( 0  U22 ) (  0   D  ) ( U12' U22' )
!
!  A  =  ( L11  0 ) (  D   0  ) ( L11' L21' )  if UPLO = 'L'
!        ( L21  I ) (  0  A22 ) (  0    I   )
!
!  where the order of D is at most NB. The actual order is returned in
!  the argument KB, and is either NB or NB-1, or N if N <= NB.
!
!  DLASYF is an auxiliary routine called by DSYTRF. It uses blocked code
!  (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
!  A22 (if UPLO = 'L').
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NB      (input) INTEGER
!          The maximum number of columns of the matrix A that should be
!          factored.  NB should be at least 2 to allow for 2-by-2 pivot
!          blocks.
!
!  KB      (output) INTEGER
!          The number of columns of A that were actually factored.
!          KB is either NB-1 or NB, or N if N <= NB.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, A contains details of the partial factorization.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If UPLO = 'U', only the last KB elements of IPIV are set;
!          if UPLO = 'L', only the first KB elements are set.
!
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!
!  W       (workspace) DOUBLE PRECISION array, dimension (LDW,NB)
!
!  LDW     (input) INTEGER
!          The leading dimension of the array W.  LDW >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!               has been completed, but the block diagonal matrix D is
!               exactly singular.
!
!  =====================================================================
!
!     .. Parameters ..
double precision   zero, one
parameter          ( zero = 0.0d+0, one = 1.0d+0 )
double precision   eight, sevten
parameter          ( eight = 8.0d+0, sevten = 17.0d+0 )
!     ..
!     .. Local Scalars ..
integer            imax, j, jb, jj, jmax, jp, k, kk, kkw, kp, &
                   kstep, kw
double precision   absakk, alpha, colmax, d11, d21, d22, r1, &
                   rowmax, t
!     ..
!     .. External Functions ..
logical            lsame
external           lsame
!     ..
!     .. External Subroutines ..
external           dgemv
!     ..
!     .. Intrinsic Functions ..
intrinsic          abs, max, min, sqrt
!     ..
!     .. Executable Statements ..
!
info = 0
!
!     Initialize ALPHA for use in choosing pivot block size.
!
alpha = ( one+sqrt( sevten ) ) / eight
!
if( lsame( uplo, 'U' ) ) then
!
!        Factorize the trailing columns of A using the upper triangle
!        of A and working backwards, and compute the matrix W = U12*D
!        for use in updating A11
!
!        K is the main loop index, decreasing from N in steps of 1 or 2
!
!        KW is the column of W which corresponds to column K of A
!
   k = n
10    continue
   kw = nb + k - n
!
!        Exit from loop
!
   if( ( k.le.n-nb+1 .and. nb.lt.n ) .or. k.lt.1 ) &
      go to 30
!
!        Copy column K of A to column KW of W and update it
!
   call dcopy( k, a( 1, k ), 1, w( 1, kw ), 1 )
   if( k.lt.n ) &
      call dgemv( 'No transpose', k, n-k, -one, a( 1, k+1 ), lda, &
                  w( k, kw+1 ), ldw, one, w( 1, kw ), 1 )
!
   kstep = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
   absakk = abs( w( k, kw ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
   if( k.gt.1 ) then
      imax = idamax( k-1, w( 1, kw ), 1 )
      colmax = abs( w( imax, kw ) )
   else
      colmax = zero
   end if
!
   if( max( absakk, colmax ).eq.zero ) then
!
!           Column K is zero: set INFO and continue
!
      if( info.eq.0 ) &
         info = k
      kp = k
   else
      if( absakk.ge.alpha*colmax ) then
!
!              no interchange, use 1-by-1 pivot block
!
         kp = k
      else
!
!              Copy column IMAX to column KW-1 of W and update it
!
         call dcopy( imax, a( 1, imax ), 1, w( 1, kw-1 ), 1 )
         call dcopy( k-imax, a( imax, imax+1 ), lda, &
                     w( imax+1, kw-1 ), 1 )
         if( k.lt.n ) &
            call dgemv( 'No transpose', k, n-k, -one, a( 1, k+1 ), &
                        lda, w( imax, kw+1 ), ldw, one, &
                        w( 1, kw-1 ), 1 )
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
         jmax = imax + idamax( k-imax, w( imax+1, kw-1 ), 1 )
         rowmax = abs( w( jmax, kw-1 ) )
         if( imax.gt.1 ) then
            jmax = idamax( imax-1, w( 1, kw-1 ), 1 )
            rowmax = max( rowmax, abs( w( jmax, kw-1 ) ) )
         end if
!
         if( absakk.ge.alpha*colmax*( colmax / rowmax ) ) then
!
!                 no interchange, use 1-by-1 pivot block
!
            kp = k
         else if( abs( w( imax, kw-1 ) ).ge.alpha*rowmax ) then
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
            kp = imax
!
!                 copy column KW-1 of W to column KW
!
            call dcopy( k, w( 1, kw-1 ), 1, w( 1, kw ), 1 )
         else
!
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!
            kp = imax
            kstep = 2
         end if
      end if
!
      kk = k - kstep + 1
      kkw = nb + kk - n
!
!           Updated column KP is already stored in column KKW of W
!
      if( kp.ne.kk ) then
!
!              Copy non-updated column KK to column KP
!
         a( kp, k ) = a( kk, k )
         call dcopy( k-1-kp, a( kp+1, kk ), 1, a( kp, kp+1 ), &
                     lda )
         call dcopy( kp, a( 1, kk ), 1, a( 1, kp ), 1 )
!
!              Interchange rows KK and KP in last KK columns of A and W
!
         call dswap( n-kk+1, a( kk, kk ), lda, a( kp, kk ), lda )
         call dswap( n-kk+1, w( kk, kkw ), ldw, w( kp, kkw ), &
                     ldw )
      end if
!
      if( kstep.eq.1 ) then
!
!              1-by-1 pivot block D(k): column KW of W now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
!              Store U(k) in column k of A
!
         call dcopy( k, w( 1, kw ), 1, a( 1, k ), 1 )
         r1 = one / a( k, k )
         call dscal( k-1, r1, a( 1, k ), 1 )
      else
!
!              2-by-2 pivot block D(k): columns KW and KW-1 of W now
!              hold
!
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
         if( k.gt.2 ) then
!
!                 Store U(k) and U(k-1) in columns k and k-1 of A
!
            d21 = w( k-1, kw )
            d11 = w( k, kw ) / d21
            d22 = w( k-1, kw-1 ) / d21
            t = one / ( d11*d22-one )
            d21 = t / d21
            do 20 j = 1, k - 2
               a( j, k-1 ) = d21*( d11*w( j, kw-1 )-w( j, kw ) )
               a( j, k ) = d21*( d22*w( j, kw )-w( j, kw-1 ) )
20             continue
         end if
!
!              Copy D(k) to A
!
         a( k-1, k-1 ) = w( k-1, kw-1 )
         a( k-1, k ) = w( k-1, kw )
         a( k, k ) = w( k, kw )
      end if
   end if
!
!        Store details of the interchanges in IPIV
!
   if( kstep.eq.1 ) then
      ipiv( k ) = kp
   else
      ipiv( k ) = -kp
      ipiv( k-1 ) = -kp
   end if
!
!        Decrease K and return to the start of the main loop
!
   k = k - kstep
   go to 10
!
30    continue
!
!        Update the upper triangle of A11 (= A(1:k,1:k)) as
!
!        A11 := A11 - U12*D*U12' = A11 - U12*W'
!
!        computing blocks of NB columns at a time
!
   do 50 j = ( ( k-1 ) / nb )*nb + 1, 1, -nb
      jb = min( nb, k-j+1 )
!
!           Update the upper triangle of the diagonal block
!
      do 40 jj = j, j + jb - 1
         call dgemv( 'No transpose', jj-j+1, n-k, -one, &
                     a( j, k+1 ), lda, w( jj, kw+1 ), ldw, one, &
                     a( j, jj ), 1 )
40       continue
!
!           Update the rectangular superdiagonal block
!
      call dgemm( 'No transpose', 'Transpose', j-1, jb, n-k, -one, &
                  a( 1, k+1 ), lda, w( j, kw+1 ), ldw, one, &
                  a( 1, j ), lda )
50    continue
!
!        Put U12 in standard form by partially undoing the interchanges
!        in columns k+1:n
!
   j = k + 1
60    continue
   jj = j
   jp = ipiv( j )
   if( jp.lt.0 ) then
      jp = -jp
      j = j + 1
   end if
   j = j + 1
   if( jp.ne.jj .and. j.le.n ) &
      call dswap( n-j+1, a( jp, j ), lda, a( jj, j ), lda )
   if( j.le.n ) &
      go to 60
!
!        Set KB to the number of columns factorized
!
   kb = n - k
!
else
!
!        Factorize the leading columns of A using the lower triangle
!        of A and working forwards, and compute the matrix W = L21*D
!        for use in updating A22
!
!        K is the main loop index, increasing from 1 in steps of 1 or 2
!
   k = 1
70    continue
!
!        Exit from loop
!
   if( ( k.ge.nb .and. nb.lt.n ) .or. k.gt.n ) &
      go to 90
!
!        Copy column K of A to column K of W and update it
!
   call dcopy( n-k+1, a( k, k ), 1, w( k, k ), 1 )
   call dgemv( 'No transpose', n-k+1, k-1, -one, a( k, 1 ), lda, &
               w( k, 1 ), ldw, one, w( k, k ), 1 )
!
   kstep = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
   absakk = abs( w( k, k ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
   if( k.lt.n ) then
      imax = k + idamax( n-k, w( k+1, k ), 1 )
      colmax = abs( w( imax, k ) )
   else
      colmax = zero
   end if
!
   if( max( absakk, colmax ).eq.zero ) then
!
!           Column K is zero: set INFO and continue
!
      if( info.eq.0 ) &
         info = k
      kp = k
   else
      if( absakk.ge.alpha*colmax ) then
!
!              no interchange, use 1-by-1 pivot block
!
         kp = k
      else
!
!              Copy column IMAX to column K+1 of W and update it
!
         call dcopy( imax-k, a( imax, k ), lda, w( k, k+1 ), 1 )
         call dcopy( n-imax+1, a( imax, imax ), 1, w( imax, k+1 ), &
                     1 )
         call dgemv( 'No transpose', n-k+1, k-1, -one, a( k, 1 ), &
                     lda, w( imax, 1 ), ldw, one, w( k, k+1 ), 1 )
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
         jmax = k - 1 + idamax( imax-k, w( k, k+1 ), 1 )
         rowmax = abs( w( jmax, k+1 ) )
         if( imax.lt.n ) then
            jmax = imax + idamax( n-imax, w( imax+1, k+1 ), 1 )
            rowmax = max( rowmax, abs( w( jmax, k+1 ) ) )
         end if
!
         if( absakk.ge.alpha*colmax*( colmax / rowmax ) ) then
!
!                 no interchange, use 1-by-1 pivot block
!
            kp = k
         else if( abs( w( imax, k+1 ) ).ge.alpha*rowmax ) then
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
            kp = imax
!
!                 copy column K+1 of W to column K
!
            call dcopy( n-k+1, w( k, k+1 ), 1, w( k, k ), 1 )
         else
!
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!
            kp = imax
            kstep = 2
         end if
      end if
!
      kk = k + kstep - 1
!
!           Updated column KP is already stored in column KK of W
!
      if( kp.ne.kk ) then
!
!              Copy non-updated column KK to column KP
!
         a( kp, k ) = a( kk, k )
         call dcopy( kp-k-1, a( k+1, kk ), 1, a( kp, k+1 ), lda )
         call dcopy( n-kp+1, a( kp, kk ), 1, a( kp, kp ), 1 )
!
!              Interchange rows KK and KP in first KK columns of A and W
!
         call dswap( kk, a( kk, 1 ), lda, a( kp, 1 ), lda )
         call dswap( kk, w( kk, 1 ), ldw, w( kp, 1 ), ldw )
      end if
!
      if( kstep.eq.1 ) then
!
!              1-by-1 pivot block D(k): column k of W now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
!              Store L(k) in column k of A
!
         call dcopy( n-k+1, w( k, k ), 1, a( k, k ), 1 )
         if( k.lt.n ) then
            r1 = one / a( k, k )
            call dscal( n-k, r1, a( k+1, k ), 1 )
         end if
      else
!
!              2-by-2 pivot block D(k): columns k and k+1 of W now hold
!
!              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
!
!              where L(k) and L(k+1) are the k-th and (k+1)-th columns
!              of L
!
         if( k.lt.n-1 ) then
!
!                 Store L(k) and L(k+1) in columns k and k+1 of A
!
            d21 = w( k+1, k )
            d11 = w( k+1, k+1 ) / d21
            d22 = w( k, k ) / d21
            t = one / ( d11*d22-one )
            d21 = t / d21
            do 80 j = k + 2, n
               a( j, k ) = d21*( d11*w( j, k )-w( j, k+1 ) )
               a( j, k+1 ) = d21*( d22*w( j, k+1 )-w( j, k ) )
80             continue
         end if
!
!              Copy D(k) to A
!
         a( k, k ) = w( k, k )
         a( k+1, k ) = w( k+1, k )
         a( k+1, k+1 ) = w( k+1, k+1 )
      end if
   end if
!
!        Store details of the interchanges in IPIV
!
   if( kstep.eq.1 ) then
      ipiv( k ) = kp
   else
      ipiv( k ) = -kp
      ipiv( k+1 ) = -kp
   end if
!
!        Increase K and return to the start of the main loop
!
   k = k + kstep
   go to 70
!
90    continue
!
!        Update the lower triangle of A22 (= A(k:n,k:n)) as
!
!        A22 := A22 - L21*D*L21' = A22 - L21*W'
!
!        computing blocks of NB columns at a time
!
   do 110 j = k, n, nb
      jb = min( nb, n-j+1 )
!
!           Update the lower triangle of the diagonal block
!
      do 100 jj = j, j + jb - 1
         call dgemv( 'No transpose', j+jb-jj, k-1, -one, &
                     a( jj, 1 ), lda, w( jj, 1 ), ldw, one, &
                     a( jj, jj ), 1 )
100       continue
!
!           Update the rectangular subdiagonal block
!
      if( j+jb.le.n ) &
         call dgemm( 'No transpose', 'Transpose', n-j-jb+1, jb, &
                     k-1, -one, a( j+jb, 1 ), lda, w( j, 1 ), ldw, &
                     one, a( j+jb, j ), lda )
110    continue
!
!        Put L21 in standard form by partially undoing the interchanges
!        in columns 1:k-1
!
   j = k - 1
120    continue
   jj = j
   jp = ipiv( j )
   if( jp.lt.0 ) then
      jp = -jp
      j = j - 1
   end if
   j = j - 1
   if( jp.ne.jj .and. j.ge.1 ) &
      call dswap( j, a( jp, 1 ), lda, a( jj, 1 ), lda )
   if( j.ge.1 ) &
      go to 120
!
!        Set KB to the number of columns factorized
!
   kb = k - 1
!
end if
return
!
!     End of DLASYF
!
end
#endif

#if defined(HIB_NONE)
subroutine dsytri( uplo, n, a, lda, ipiv, work, info )
use mod_hiblas, only: dswap
!
!  -- LAPACK routine (version 1.0b) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
character          uplo
integer            info, lda, n
!     ..
!     .. Array Arguments ..
integer            ipiv( * )
double precision   a( lda, * ), work( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTRI computes the inverse of a real symmetric indefinite matrix
!  A using the factorization A = U*D*U' or A = L*D*L' computed by
!  DSYTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the details of the factorization are stored
!          as an upper or lower triangular matrix.
!          = 'U':  Upper triangular (form is A = U*D*U')
!          = 'L':  Lower triangular (form is A = L*D*L')
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the block diagonal matrix D and the multipliers
!          used to obtain the factor U or L as computed by DSYTRF.
!
!          On exit, if INFO = 0, the (symmetric) inverse of the original
!          matrix.  If UPLO = 'U', the upper triangular part of the
!          inverse is formed and the part of A below the diagonal is not
!          referenced; if UPLO = 'L' the lower triangular part of the
!          inverse is formed and the part of A above the diagonal is
!          not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D
!          as determined by DSYTRF.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, D(k,k) = 0; the matrix is singular and its
!               inverse could not be computed.
!
!  =====================================================================
!
!     .. Parameters ..
double precision   one, zero
parameter          ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. Local Scalars ..
logical            upper
integer            k, kp, kstep
double precision   ak, akkp1, akp1, d, t, temp
!     ..
!     .. External Functions ..
logical            lsame
double precision   ddot
external           lsame, ddot
!     ..
!     .. External Subroutines ..
external           dcopy, dsymv, xerbla
!     ..
!     .. Intrinsic Functions ..
intrinsic          abs, max
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
info = 0
upper = lsame( uplo, 'U' )
if( .not.upper .and. .not.lsame( uplo, 'L' ) ) then
   info = -1
else if( n.lt.0 ) then
   info = -2
else if( lda.lt.max( 1, n ) ) then
   info = -4
end if
if( info.ne.0 ) then
   call xerbla( 'DSYTRI', -info )
   return
end if
!
!     Quick return if possible
!
if( n.eq.0 ) &
   return
!
!     Check that the diagonal matrix D is nonsingular.
!
if( upper ) then
!
!        Upper triangular storage: examine D from bottom to top
!
   do 10 info = n, 1, -1
      if( ipiv( info ).gt.0 .and. a( info, info ).eq.zero ) &
         return
10    continue
else
!
!        Lower triangular storage: examine D from top to bottom.
!
   do 20 info = 1, n
      if( ipiv( info ).gt.0 .and. a( info, info ).eq.zero ) &
         return
20    continue
end if
info = 0
!
if( upper ) then
!
!        Compute inv(A) from the factorization A = U*D*U'.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
   k = 1
30    continue
!
!        If K > N, exit from loop.
!
   if( k.gt.n ) &
      go to 40
!
   if( ipiv( k ).gt.0 ) then
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
      a( k, k ) = one / a( k, k )
!
!           Compute column K of the inverse.
!
      if( k.gt.1 ) then
         call dcopy( k-1, a( 1, k ), 1, work, 1 )
         call dsymv( uplo, k-1, -one, a, lda, work, 1, zero, &
                     a( 1, k ), 1 )
         a( k, k ) = a( k, k ) - ddot( k-1, work, 1, a( 1, k ), &
                     1 )
      end if
      kstep = 1
   else
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
      t = abs( a( k, k+1 ) )
      ak = a( k, k ) / t
      akp1 = a( k+1, k+1 ) / t
      akkp1 = a( k, k+1 ) / t
      d = t*( ak*akp1-one )
      a( k, k ) = akp1 / d
      a( k+1, k+1 ) = ak / d
      a( k, k+1 ) = -akkp1 / d
!
!           Compute columns K and K+1 of the inverse.
!
      if( k.gt.1 ) then
         call dcopy( k-1, a( 1, k ), 1, work, 1 )
         call dsymv( uplo, k-1, -one, a, lda, work, 1, zero, &
                     a( 1, k ), 1 )
         a( k, k ) = a( k, k ) - ddot( k-1, work, 1, a( 1, k ), &
                     1 )
         a( k, k+1 ) = a( k, k+1 ) - &
                       ddot( k-1, a( 1, k ), 1, a( 1, k+1 ), 1 )
         call dcopy( k-1, a( 1, k+1 ), 1, work, 1 )
         call dsymv( uplo, k-1, -one, a, lda, work, 1, zero, &
                     a( 1, k+1 ), 1 )
         a( k+1, k+1 ) = a( k+1, k+1 ) - &
                         ddot( k-1, work, 1, a( 1, k+1 ), 1 )
      end if
      kstep = 2
   end if
!
   kp = abs( ipiv( k ) )
   if( kp.ne.k ) then
!
!           Interchange rows and columns K and KP in the leading
!           submatrix A(1:k+1,1:k+1)
!
      call dswap( kp-1, a( 1, k ), 1, a( 1, kp ), 1 )
      call dswap( k-kp-1, a( kp+1, k ), 1, a( kp, kp+1 ), lda )
      temp = a( k, k )
      a( k, k ) = a( kp, kp )
      a( kp, kp ) = temp
      if( kstep.eq.2 ) then
         temp = a( k, k+1 )
         a( k, k+1 ) = a( kp, k+1 )
         a( kp, k+1 ) = temp
      end if
   end if
!
   k = k + kstep
   go to 30
40    continue
!
else
!
!        Compute inv(A) from the factorization A = L*D*L'.
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2, depending on the size of the diagonal blocks.
!
   k = n
50    continue
!
!        If K < 1, exit from loop.
!
   if( k.lt.1 ) &
      go to 60
!
   if( ipiv( k ).gt.0 ) then
!
!           1 x 1 diagonal block
!
!           Invert the diagonal block.
!
      a( k, k ) = one / a( k, k )
!
!           Compute column K of the inverse.
!
      if( k.lt.n ) then
         call dcopy( n-k, a( k+1, k ), 1, work, 1 )
         call dsymv( uplo, n-k, -one, a( k+1, k+1 ), lda, work, 1, &
                     zero, a( k+1, k ), 1 )
         a( k, k ) = a( k, k ) - ddot( n-k, work, 1, a( k+1, k ), &
                     1 )
      end if
      kstep = 1
   else
!
!           2 x 2 diagonal block
!
!           Invert the diagonal block.
!
      t = abs( a( k, k-1 ) )
      ak = a( k-1, k-1 ) / t
      akp1 = a( k, k ) / t
      akkp1 = a( k, k-1 ) / t
      d = t*( ak*akp1-one )
      a( k-1, k-1 ) = akp1 / d
      a( k, k ) = ak / d
      a( k, k-1 ) = -akkp1 / d
!
!           Compute columns K-1 and K of the inverse.
!
      if( k.lt.n ) then
         call dcopy( n-k, a( k+1, k ), 1, work, 1 )
         call dsymv( uplo, n-k, -one, a( k+1, k+1 ), lda, work, 1, &
                     zero, a( k+1, k ), 1 )
         a( k, k ) = a( k, k ) - ddot( n-k, work, 1, a( k+1, k ), &
                     1 )
         a( k, k-1 ) = a( k, k-1 ) - &
                       ddot( n-k, a( k+1, k ), 1, a( k+1, k-1 ), &
                       1 )
         call dcopy( n-k, a( k+1, k-1 ), 1, work, 1 )
         call dsymv( uplo, n-k, -one, a( k+1, k+1 ), lda, work, 1, &
                     zero, a( k+1, k-1 ), 1 )
         a( k-1, k-1 ) = a( k-1, k-1 ) - &
                         ddot( n-k, work, 1, a( k+1, k-1 ), 1 )
      end if
      kstep = 2
   end if
!
   kp = abs( ipiv( k ) )
   if( kp.ne.k ) then
!
!           Interchange rows and columns K and KP in the trailing
!           submatrix A(k-1:n,k-1:n)
!
      if( kp.lt.n ) &
         call dswap( n-kp, a( kp+1, k ), 1, a( kp+1, kp ), 1 )
      call dswap( kp-k-1, a( k+1, k ), 1, a( kp, k+1 ), lda )
      temp = a( k, k )
      a( k, k ) = a( kp, kp )
      a( kp, kp ) = temp
      if( kstep.eq.2 ) then
         temp = a( k, k-1 )
         a( k, k-1 ) = a( kp, k-1 )
         a( kp, k-1 ) = temp
      end if
   end if
!
   k = k - kstep
   go to 50
60    continue
end if
!
return
!
!     End of DSYTRI
!
end
subroutine dsytrf( uplo, n, a, lda, ipiv, work, lwork, info )
!
!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
character          uplo
integer            info, lda, lwork, n
!     ..
!     .. Array Arguments ..
integer            ipiv( * )
double precision   a( lda, * ), work( lwork )
!     ..
!
!  Purpose
!  =======
!
!  DSYTRF computes the factorization of a real symmetric matrix A using
!  the Bunch-Kaufman diagonal pivoting method:
!
!     A = U*D*U'  or  A = L*D*L'
!
!  where U (or L) is a product of permutation and unit upper (lower)
!  triangular matrices, U' is the transpose of U, and D is symmetric and
!  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, the block diagonal matrix D and the multipliers used
!          to obtain the factor U or L (see below for further details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
!          If INFO returns 0, then WORK(1) returns N*NB, the minimum
!          value of LWORK required to use the optimal blocksize.
!
!  LWORK   (input) INTEGER
!          The length of WORK.  LWORK should be >= N*NB, where NB is the
!          block size returned by ILAENV.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!               has been completed, but the block diagonal matrix D is
!               exactly singular, and division by zero will occur if it
!               is used to solve a system of equations.
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', then A = U*D*U', where
!     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    v    0   )   k-s
!     U(k) =  (   0    I    0   )   s
!             (   0    0    I   )   n-k
!                k-s   s   n-k
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!
!  If UPLO = 'L', then A = L*D*L', where
!     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    0     0   )  k-1
!     L(k) =  (   0    I     0   )  s
!             (   0    v     I   )  n-k-s+1
!                k-1   s  n-k-s+1
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!
!  =====================================================================
!
!     .. Local Scalars ..
logical            upper
integer            iinfo, iws, j, k, kb, ldwork, nb, nbmin
!     ..
!     .. External Functions ..
logical            lsame
integer            ilaenv
external           lsame, ilaenv
!     ..
!     .. External Subroutines ..
external           dlasyf, dsytf2, xerbla
!     ..
!     .. Intrinsic Functions ..
intrinsic          max
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
info = 0
upper = lsame( uplo, 'U' )
if( .not.upper .and. .not.lsame( uplo, 'L' ) ) then
   info = -1
else if( n.lt.0 ) then
   info = -2
else if( lda.lt.max( 1, n ) ) then
   info = -4
else if( lwork.lt.1 ) then
   info = -7
end if
if( info.ne.0 ) then
   call xerbla( 'DSYTRF', -info )
   return
end if
!
!     Determine the block size
!
nb = ilaenv( 1, 'DSYTRF', uplo, n, -1, -1, -1 )
nbmin = 2
ldwork = n
if( nb.gt.1 .and. nb.lt.n ) then
   iws = ldwork*nb
   if( lwork.lt.iws ) then
      nb = max( lwork / ldwork, 1 )
      nbmin = max( 2, ilaenv( 2, 'DSYTRF', uplo, n, -1, -1, -1 ) )
   end if
else
   iws = 1
end if
if( nb.lt.nbmin ) &
   nb = n
!
if( upper ) then
!
!        Factorize A as U*D*U' using the upper triangle of A
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        KB, where KB is the number of columns factorized by DLASYF;
!        KB is either NB or NB-1, or K for the last block
!
   k = n
10    continue
!
!        If K < 1, exit from loop
!
   if( k.lt.1 ) &
      go to 40
!
   if( k.gt.nb ) then
!
!           Factorize columns k-kb+1:k of A and use blocked code to
!           update columns 1:k-kb
!
      call dlasyf( uplo, k, nb, kb, a, lda, ipiv, work, ldwork, &
                   iinfo )
   else
!
!           Use unblocked code to factorize columns 1:k of A
!
      call dsytf2( uplo, k, a, lda, ipiv, iinfo )
      kb = k
   end if
!
!        Set INFO on the first occurrence of a zero pivot
!
   if( info.eq.0 .and. iinfo.gt.0 ) &
      info = iinfo
!
!        Decrease K and return to the start of the main loop
!
   k = k - kb
   go to 10
!
else
!
!        Factorize A as L*D*L' using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        KB, where KB is the number of columns factorized by DLASYF;
!        KB is either NB or NB-1, or N-K+1 for the last block
!
   k = 1
20    continue
!
!        If K > N, exit from loop
!
   if( k.gt.n ) &
      go to 40
!
   if( k.le.n-nb ) then
!
!           Factorize columns k:k+kb-1 of A and use blocked code to
!           update columns k+kb:n
!
      call dlasyf( uplo, n-k+1, nb, kb, a( k, k ), lda, ipiv( k ), &
                   work, ldwork, iinfo )
   else
!
!           Use unblocked code to factorize columns k:n of A
!
      call dsytf2( uplo, n-k+1, a( k, k ), lda, ipiv( k ), iinfo )
      kb = n - k + 1
   end if
!
!        Set INFO on the first occurrence of a zero pivot
!
   if( info.eq.0 .and. iinfo.gt.0 ) &
      info = iinfo + k - 1
!
!        Adjust IPIV
!
   do 30 j = k, k + kb - 1
      if( ipiv( j ).gt.0 ) then
         ipiv( j ) = ipiv( j ) + k - 1
      else
         ipiv( j ) = ipiv( j ) - k + 1
      end if
30    continue
!
!        Increase K and return to the start of the main loop
!
   k = k + kb
   go to 20
!
end if
!
40 continue
work( 1 ) = iws
return
!
!     End of DSYTRF
!
end
subroutine dsytf2( uplo, n, a, lda, ipiv, info )
use mod_hiblas, only: dscal, drot, dswap, idamax
!
!  -- LAPACK routine (version 1.0b) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
character          uplo
integer            info, lda, n
!     ..
!     .. Array Arguments ..
integer            ipiv( * )
double precision   a( lda, * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTF2 computes the factorization of a real symmetric matrix A using
!  the Bunch-Kaufman diagonal pivoting method:
!
!     A = U*D*U'  or  A = L*D*L'
!
!  where U (or L) is a product of permutation and unit upper (lower)
!  triangular matrices, U' is the transpose of U, and D is symmetric and
!  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!
!          On exit, the block diagonal matrix D and the multipliers used
!          to obtain the factor U or L (see below for further details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (output) INTEGER array, dimension (N)
!          Details of the interchanges and the block structure of D.
!          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!          interchanged and D(k,k) is a 1-by-1 diagonal block.
!          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
!               has been completed, but the block diagonal matrix D is
!               exactly singular, and division by zero will occur if it
!               is used to solve a system of equations.
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', then A = U*D*U', where
!     U = P(n)*U(n)* ... *P(k)U(k)* ...,
!  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
!  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    v    0   )   k-s
!     U(k) =  (   0    I    0   )   s
!             (   0    0    I   )   n-k
!                k-s   s   n-k
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
!  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
!  and A(k,k), and v overwrites A(1:k-2,k-1:k).
!
!  If UPLO = 'L', then A = L*D*L', where
!     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
!  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
!  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
!  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
!  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
!  that if the diagonal block D(k) is of order s (s = 1 or 2), then
!
!             (   I    0     0   )  k-1
!     L(k) =  (   0    I     0   )  s
!             (   0    v     I   )  n-k-s+1
!                k-1   s  n-k-s+1
!
!  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
!  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
!  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
!
!  =====================================================================
!
!     .. Parameters ..
double precision   zero, one
parameter          ( zero = 0.0d+0, one = 1.0d+0 )
double precision   eight, sevten
parameter          ( eight = 8.0d+0, sevten = 17.0d+0 )
!     ..
!     .. Local Scalars ..
logical            upper
integer            imax, jmax, k, kk, kp, kstep
double precision   absakk, alpha, c, colmax, r1, r2, rowmax, s, t
!     ..
!     .. External Functions ..
logical            lsame
external           lsame
!     ..
!     .. External Subroutines ..
external           dlaev2, dsyr, xerbla
!     ..
!     .. Intrinsic Functions ..
intrinsic          abs, max, sqrt
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
info = 0
upper = lsame( uplo, 'U' )
if( .not.upper .and. .not.lsame( uplo, 'L' ) ) then
   info = -1
else if( n.lt.0 ) then
   info = -2
else if( lda.lt.max( 1, n ) ) then
   info = -4
end if
if( info.ne.0 ) then
   call xerbla( 'DSYTF2', -info )
   return
end if
!
!     Initialize ALPHA for use in choosing pivot block size.
!
alpha = ( one+sqrt( sevten ) ) / eight
!
if( upper ) then
!
!        Factorize A as U*D*U' using the upper triangle of A
!
!        K is the main loop index, decreasing from N to 1 in steps of
!        1 or 2
!
   k = n
10    continue
!
!        If K < 1, exit from loop
!
   if( k.lt.1 ) &
      go to 30
   kstep = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
   absakk = abs( a( k, k ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
   if( k.gt.1 ) then
      imax = idamax( k-1, a( 1, k ), 1 )
      colmax = abs( a( imax, k ) )
   else
      colmax = zero
   end if
!
   if( max( absakk, colmax ).eq.zero ) then
!
!           Column K is zero: set INFO and continue
!
      if( info.eq.0 ) &
         info = k
      kp = k
   else
      if( absakk.ge.alpha*colmax ) then
!
!              no interchange, use 1-by-1 pivot block
!
         kp = k
      else
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
         jmax = imax + idamax( k-imax, a( imax, imax+1 ), lda )
         rowmax = abs( a( imax, jmax ) )
         if( imax.gt.1 ) then
            jmax = idamax( imax-1, a( 1, imax ), 1 )
            rowmax = max( rowmax, abs( a( jmax, imax ) ) )
         end if
!
         if( absakk.ge.alpha*colmax*( colmax / rowmax ) ) then
!
!                 no interchange, use 1-by-1 pivot block
!
            kp = k
         else if( abs( a( imax, imax ) ).ge.alpha*rowmax ) then
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
            kp = imax
         else
!
!                 interchange rows and columns K-1 and IMAX, use 2-by-2
!                 pivot block
!
            kp = imax
            kstep = 2
         end if
      end if
!
      kk = k - kstep + 1
      if( kp.ne.kk ) then
!
!              Interchange rows and columns KK and KP in the leading
!              submatrix A(1:k,1:k)
!
         call dswap( kp-1, a( 1, kk ), 1, a( 1, kp ), 1 )
         call dswap( kk-kp-1, a( kp+1, kk ), 1, a( kp, kp+1 ), &
                     lda )
         t = a( kk, kk )
         a( kk, kk ) = a( kp, kp )
         a( kp, kp ) = t
         if( kstep.eq.2 ) then
            t = a( k-1, k )
            a( k-1, k ) = a( kp, k )
            a( kp, k ) = t
         end if
      end if
!
!           Update the leading submatrix
!
      if( kstep.eq.1 ) then
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = U(k)*D(k)
!
!              where U(k) is the k-th column of U
!
!              Perform a rank-1 update of A(1:k-1,1:k-1) as
!
!              A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
!
         r1 = one / a( k, k )
         call dsyr( uplo, k-1, -r1, a( 1, k ), 1, a, lda )
!
!              Store U(k) in column k
!
         call dscal( k-1, r1, a( 1, k ), 1 )
      else
!
!              2-by-2 pivot block D(k): columns k and k-1 now hold
!
!              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
!
!              where U(k) and U(k-1) are the k-th and (k-1)-th columns
!              of U
!
!              Perform a rank-2 update of A(1:k-2,1:k-2) as
!
!              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
!                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
!
!              Convert this to two rank-1 updates by using the eigen-
!              decomposition of D(k)
!
         call dlaev2( a( k-1, k-1 ), a( k-1, k ), a( k, k ), r1, &
                      r2, c, s )
         r1 = one / r1
         r2 = one / r2
         call drot( k-2, a( 1, k-1 ), 1, a( 1, k ), 1, c, s )
         call dsyr( uplo, k-2, -r1, a( 1, k-1 ), 1, a, lda )
         call dsyr( uplo, k-2, -r2, a( 1, k ), 1, a, lda )
!
!              Store U(k) and U(k-1) in columns k and k-1
!
         call dscal( k-2, r1, a( 1, k-1 ), 1 )
         call dscal( k-2, r2, a( 1, k ), 1 )
         call drot( k-2, a( 1, k-1 ), 1, a( 1, k ), 1, c, -s )
      end if
   end if
!
!        Store details of the interchanges in IPIV
!
   if( kstep.eq.1 ) then
      ipiv( k ) = kp
   else
      ipiv( k ) = -kp
      ipiv( k-1 ) = -kp
   end if
!
!        Decrease K and return to the start of the main loop
!
   k = k - kstep
   go to 10
!
else
!
!        Factorize A as L*D*L' using the lower triangle of A
!
!        K is the main loop index, increasing from 1 to N in steps of
!        1 or 2
!
   k = 1
20    continue
!
!        If K > N, exit from loop
!
   if( k.gt.n ) &
      go to 30
   kstep = 1
!
!        Determine rows and columns to be interchanged and whether
!        a 1-by-1 or 2-by-2 pivot block will be used
!
   absakk = abs( a( k, k ) )
!
!        IMAX is the row-index of the largest off-diagonal element in
!        column K, and COLMAX is its absolute value
!
   if( k.lt.n ) then
      imax = k + idamax( n-k, a( k+1, k ), 1 )
      colmax = abs( a( imax, k ) )
   else
      colmax = zero
   end if
!
   if( max( absakk, colmax ).eq.zero ) then
!
!           Column K is zero: set INFO and continue
!
      if( info.eq.0 ) &
         info = k
      kp = k
   else
      if( absakk.ge.alpha*colmax ) then
!
!              no interchange, use 1-by-1 pivot block
!
         kp = k
      else
!
!              JMAX is the column-index of the largest off-diagonal
!              element in row IMAX, and ROWMAX is its absolute value
!
         jmax = k - 1 + idamax( imax-k, a( imax, k ), lda )
         rowmax = abs( a( imax, jmax ) )
         if( imax.lt.n ) then
            jmax = imax + idamax( n-imax, a( imax+1, imax ), 1 )
            rowmax = max( rowmax, abs( a( jmax, imax ) ) )
         end if
!
         if( absakk.ge.alpha*colmax*( colmax / rowmax ) ) then
!
!                 no interchange, use 1-by-1 pivot block
!
            kp = k
         else if( abs( a( imax, imax ) ).ge.alpha*rowmax ) then
!
!                 interchange rows and columns K and IMAX, use 1-by-1
!                 pivot block
!
            kp = imax
         else
!
!                 interchange rows and columns K+1 and IMAX, use 2-by-2
!                 pivot block
!
            kp = imax
            kstep = 2
         end if
      end if
!
      kk = k + kstep - 1
      if( kp.ne.kk ) then
!
!              Interchange rows and columns KK and KP in the trailing
!              submatrix A(k:n,k:n)
!
         if( kp.lt.n ) &
            call dswap( n-kp, a( kp+1, kk ), 1, a( kp+1, kp ), 1 )
         call dswap( kp-kk-1, a( kk+1, kk ), 1, a( kp, kk+1 ), &
                     lda )
         t = a( kk, kk )
         a( kk, kk ) = a( kp, kp )
         a( kp, kp ) = t
         if( kstep.eq.2 ) then
            t = a( k+1, k )
            a( k+1, k ) = a( kp, k )
            a( kp, k ) = t
         end if
      end if
!
!           Update the trailing submatrix
!
      if( kstep.eq.1 ) then
!
!              1-by-1 pivot block D(k): column k now holds
!
!              W(k) = L(k)*D(k)
!
!              where L(k) is the k-th column of L
!
         if( k.lt.n ) then
!
!                 Perform a rank-1 update of A(k+1:n,k+1:n) as
!
!                 A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
!
            r1 = one / a( k, k )
            call dsyr( uplo, n-k, -r1, a( k+1, k ), 1, &
                       a( k+1, k+1 ), lda )
!
!                 Store L(k) in column K
!
            call dscal( n-k, r1, a( k+1, k ), 1 )
         end if
      else
!
!              2-by-2 pivot block D(k): columns K and K+1 now hold
!
!              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
!
!              where L(k) and L(k+1) are the k-th and (k+1)-th columns
!              of L
!
         if( k.lt.n-1 ) then
!
!                 Perform a rank-2 update of A(k+2:n,k+2:n) as
!
!                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )'
!                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )'
!
!                 Convert this to two rank-1 updates by using the eigen-
!                 decomposition of D(k)
!
            call dlaev2( a( k, k ), a( k+1, k ), a( k+1, k+1 ), &
                         r1, r2, c, s )
            r1 = one / r1
            r2 = one / r2
            call drot( n-k-1, a( k+2, k ), 1, a( k+2, k+1 ), 1, c, &
                       s )
            call dsyr( uplo, n-k-1, -r1, a( k+2, k ), 1, &
                       a( k+2, k+2 ), lda )
            call dsyr( uplo, n-k-1, -r2, a( k+2, k+1 ), 1, &
                       a( k+2, k+2 ), lda )
!
!                 Store L(k) and L(k+1) in columns k and k+1
!
            call dscal( n-k-1, r1, a( k+2, k ), 1 )
            call dscal( n-k-1, r2, a( k+2, k+1 ), 1 )
            call drot( n-k-1, a( k+2, k ), 1, a( k+2, k+1 ), 1, c, &
                       -s )
         end if
      end if
   end if
!
!        Store details of the interchanges in IPIV
!
   if( kstep.eq.1 ) then
      ipiv( k ) = kp
   else
      ipiv( k ) = -kp
      ipiv( k+1 ) = -kp
   end if
!
!        Increase K and return to the start of the main loop
!
   k = k + kstep
   go to 20
!
end if
!
30 continue
return
!
!     End of DSYTF2
!
end
#endif

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
