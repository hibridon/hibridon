module mod_hivector
contains
!*************************************************************************
!                                                                        *
!                    vector routines library supplement             *
!                                                                        *
!*************************************************************************
!                        routines included                               *
!                                                                        *
!  1. matcopy   puts a nr x nc matrix a into nr x nc matrix b             *
!  2. maxmgv   finds largest value in a vector                           *
!  3. vsmul    multiplies the elements of a vector by a scalar           *
!  4. vmul     multiplies the elements of two vectors                    *
!  4a.dsum     sum of elements of a vector
!  5. blas     basis linear algebra routines from lapack
!     lsame,ilaenv,xerbla
!  6. blas extensions from lapack
!     dlaev2, dlasyf
!     isamax, saxpy, scopy, sdot, sscal, sswap, srot (grot)
!  7. fzero (vector zero) (no longer part of code)
!  8. blas extensions from ibm essl
!     idamin, idmin
!                                                                        *
!*************************************************************************
#if defined(HIB_UNIX) || defined(HIB_CRAY)
subroutine matcopy (a, b, nr, nc, na, nb)
use mod_hiblas, only: dcopy
!  subroutine to put nr x nc matrix a into nr x nc matrix b
!  author:  millard alexander
!  current revision date: 24-sept-87
! ------------------------------------------------------------------
!  variables in call list
!    a,b:     input matrices, stored as one-dimensional arrays
!    nr:      actual row dimension of matrix a
!    nc:      actual column dimension of matrix a
!    na:      maximum row dimension of matrix a
!    nb:      maximum row dimension of matrix b
! ------------------------------------------------------------------
!  the two matrices are treated as vectors here, with column-by-column
!  ordering assumed
!  the coding uses the blas routine scopy
! ------------------------------------------------------------------
implicit double precision (a-h,o-z)
integer ia, ib, j, na, nb, nr, nc
dimension a(1), b(1)
ia = 0
ib = 0
!  ia and ib point to one position before the top of the jth
!  column of each matrix
do 20  j = 1, nc
#endif
#if defined(HIB_UNIX_CONVEX)
call scopy (nr, a(ia+1), 1, b(ib+1), 1)
#endif
#if (defined(HIB_UNIX) || defined(HIB_CRAY)) && !defined(HIB_UNIX_CONVEX)
  call dcopy (nr, a(ia+1), 1, b(ib+1), 1)
#endif
#if defined(HIB_UNIX) || defined(HIB_CRAY)
  ia = ia + na
  ib = ib + nb
20 continue
return
end
!  -----------------------------------------------------------------------
subroutine maxmgv (a, na, c, nc, n)
!  subroutine to scan a vector for its maximum magnitude (absolute value)
!  element
!  current revision date: 24-sept-87
!  -----------------------------------------------------------------------
!  variables in call list:
!  a:   floating point input vector
!  na:  integer element step for a
!  c:   floating point output scalar: on return contains value of
!       maximum magnitude (absolute value) element
!  nc:  integer index of maximum magnitude element
!  n:   integer element count
!  subroutines called:
!  isamax: blas routine to find index of maximum magnitude (absolute value)
!          element
!  -----------------------------------------------------------------------
#if (defined(HIB_UNIX) || defined(HIB_CRAY)) && !defined(HIB_UNIX_CONVEX)
use mod_hiblas, only: idamax
#endif
implicit double precision (a-h,o-z)
dimension a(1)
#endif
#if defined(HIB_UNIX_CONVEX)
nc = ( isamax (n, a, na) - 1) * na + 1
#endif
#if (defined(HIB_UNIX) || defined(HIB_CRAY)) && !defined(HIB_UNIX_CONVEX)
nc = ( idamax (n, a, na) - 1) * na + 1
#endif
#if defined(HIB_UNIX) || defined(HIB_CRAY)
c = abs( a(nc) )
return
end
!  -----------------------------------------------------------------------
subroutine vsmul (a, na, b, c, nc, n)
use mod_hiblas, only: dscal, dcopy
!  subroutine to multiply the elements of a vector by a scalar
!  current revision date: 24-sept-87
!  -----------------------------------------------------------------------
!  variables in call list:
!  a:   floating point input vector
!  na:  integer element step for a
!  b:   floating point input scalar
!  c:   floating point output vector
!  nc:  integer element step for c
!  n:   integer element count
!  c(m) = a(m) * b for m=1 to n
!  -----------------------------------------------------------------------
implicit double precision (a-h,o-z)
integer n, na, nc
dimension a(1), c(1)
!  first copy vector a into vector c
!  then multiply by scalar constant
#endif
#if defined(HIB_UNIX_CONVEX)
call scopy (n, a, na, c, nc)
call sscal (n, b, c, nc)
#endif
#if defined(HIB_UNIX) || defined(HIB_CRAY) && !defined(HIB_UNIX_CONVEX)
call dcopy (n, a, na, c, nc)
call dscal (n, b, c, nc)
#endif
#if defined(HIB_UNIX) || defined(HIB_CRAY)
return
end
!  -----------------------------------------------------------------------
subroutine vmul (a, na, b, nb, c, nc, n)
!  subroutine to multiply the elements of two vectors
!  current revision date: 23-sept-87
!  -----------------------------------------------------------------------
!  variables in call list:
!  a:   floating point input vector
!  na:  integer element step for a
!  b:   floating point input vector
!  nb:  integer element step for b
!  c:   floating point output vector
!  nc:  integer element step for c
!  n:   integer element count
!  c(m) = a(m) * b(m) for m=1 to n
!  -----------------------------------------------------------------------
implicit double precision (a-h,o-z)
real(8), intent(in) :: a(:)
integer, intent(in) :: na
real(8), intent(in) :: b(:)
integer, intent(in) :: nb
real(8), intent(out) :: c(:)
integer, intent(in) :: nc
integer, intent(in) :: n
integer i, inda, indb, indc
inda = 1
indb = 1
indc = 1
do 4  i = 1, n
  c(indc) = b(indb) * a(inda)
  inda = inda + na
  indb = indb + nb
  indc = indc + nc
4 continue
return
end
!  -----------------------------------------------------------------------
subroutine vadd (ic,a, na, b, nb, n)
!  subroutine to add or subtract the elements of two vectors
!  current revision date: 6-dec-1991
!  -----------------------------------------------------------------------
!  variables in call list:
!  ic:  factor
!  a:   floating point input vector
!  na:  integer element step for a
!  b:   floating point input vector
!  nb:  integer element step for b
!  n:   integer element count
!  a(m) = a(m) + ic*b(m) for m=1 to n
!  -----------------------------------------------------------------------
implicit double precision (a-h,o-z)
integer i,ic, inda, indb, n, na, nb
dimension a(1), b(1)
inda = 1
indb = 1
if (ic .gt. 0) then
  do 4  i = 1, n
    a(inda) = a(inda) + b(indb)
    inda = inda + na
    indb = indb + nb
4   continue
else
  do 5  i = 1, n
    a(inda) = a(inda) - b(indb)
    inda = inda + na
    indb = indb + nb
5   continue
endif
return
end
double precision function dsum(n,dx,incx)
!
!     returns sum of double precision dx
!     dasum = sum from 0 to n-1 of dx(1+i*incx)
!     adapted from blas dasum by mha  4-apr-1996
!
double precision dx(1)
dsum = 0.d0
if(n.le.0)return
if(incx.eq.1)goto 20
!
!        code for increments not equal to 1
!
ns = n*incx
    do 10 i=1,ns,incx
    dsum = dsum + dx(i)
10     continue
return
!
!        code for increments equal to 1.
!
!
!        clean-up loop so remaining vector length is a multiple of 6.
!
20 m = mod(n,6)
if( m .eq. 0 ) go to 40
do 30 i = 1,m
   dsum = dsum + dx(i)
30 continue
if( n .lt. 6 ) return
40 mp1 = m + 1
do 50 i = mp1,n,6
   dsum = dsum + dx(i) + dx(i+1) + dx(i+2) &
   + dx(i+3) + dx(i+4) + dx(i+5)
50 continue
return
end
#endif
#if defined(HIB_UNIX_SUN)
logical          function lsame( ca, cb )
!
!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
character          ca, cb
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
!     .. Intrinsic Functions ..
intrinsic          ichar
!     ..
!     .. Local Scalars ..
integer            inta, intb, zcode
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
lsame = ca.eq.cb
if( lsame ) &
   return
!
!     Now test for equivalence if both characters are alphabetic.
!
zcode = ichar( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
inta = ichar( ca )
intb = ichar( cb )
!
if( zcode.eq.90 .or. zcode.eq.122 ) then
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
   if( inta.ge.97 .and. inta.le.122 ) inta = inta - 32
   if( intb.ge.97 .and. intb.le.122 ) intb = intb - 32
!
else if( zcode.eq.233 .or. zcode.eq.169 ) then
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
   if( inta.ge.129 .and. inta.le.137 .or. &
       inta.ge.145 .and. inta.le.153 .or. &
       inta.ge.162 .and. inta.le.169 ) inta = inta + 64
   if( intb.ge.129 .and. intb.le.137 .or. &
       intb.ge.145 .and. intb.le.153 .or. &
       intb.ge.162 .and. intb.le.169 ) intb = intb + 64
!
else if( zcode.eq.218 .or. zcode.eq.250 ) then
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
   if( inta.ge.225 .and. inta.le.250 ) inta = inta - 32
   if( intb.ge.225 .and. intb.le.250 ) intb = intb - 32
end if
lsame = inta.eq.intb
!
!     RETURN
!
!     End of LSAME
!
end
subroutine xerbla( srname, info )
!
!  -- LAPACK auxiliary routine (preliminary version) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
character*6        srname
integer            info
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
!
write( *, fmt = 9999 )srname, info
!
stop
!
9999 format( ' ** On entry to ', a6, ' parameter number ', i2, ' had ', &
      'an illegal value' )
!
!     End of XERBLA
!
end
#endif
#if !defined(HIB_UNIX_DEC)
subroutine dset(n,da,dx,incx)
!
!     sets a vector equal to a constant
!     n : number of elements to set
!     da : the constant value a to set
!     dx : the start of the vector x 
!     incx : the number of doubles to add to reach the next element of the vector
!     uses unrolled loops for increment equal to one.
!     modified by mha from linpack dscal, written originally by
!     jack dongarra, linpack, 3/11/78.
!     modified to correct problem with negative increment, 8/21/90.
!
double precision dx(1), da
integer i,incx,ix,m,mp1,n

!
if(n.le.0)return
if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
ix = 1
if(incx.lt.0)ix = (-n+1)*incx + 1
do 10 i = 1,n
  dx(ix) = da
  ix = ix + incx
10 continue
return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
20 m = mod(n,5)
if( m .eq. 0 ) go to 40
do 30 i = 1,m
  dx(i) = da
30 continue
if( n .lt. 5 ) return
40 mp1 = m + 1
do 50 i = mp1,n,5
  dx(i) = da
  dx(i + 1) = da
  dx(i + 2) = da
  dx(i + 3) = da
  dx(i + 4) = da
50 continue
return
end
#endif
#if defined(HIB_UNIX) && (!defined(HIB_UNIX_IBM) || !defined(HIB_UNIX_DEC))
integer function idamin(n, sx, incx)
!
!     find smallest index of minimum magnitude of double precision s
!     isamax =  first i,  i = 1 to n,  to minimize  abs(sx(1-incx+i*incx)
!
implicit none
integer, intent(in) :: n
real(8), intent(in) :: sx(n)
integer, intent(in) :: incx


real(8) :: smin, xmag
integer :: i, ii, ns

idamin = 0
if (n .le. 0) return
idamin = 1
if (n .le. 1) return
if (incx .eq. 1) go to 20
!
!        code for increments not equal to 1.
!
smin = abs(sx(1))
ns = n * incx
ii = 1
    do 10 i = 1, ns, incx
    xmag = abs(sx(i))
    if (xmag .ge. smin) go to 5
    idamin = ii
    smin = xmag
5     ii = ii + 1
10     continue
return
!
!        code for increments equal to 1.
!
20 smin = abs(sx(1))
do 30 i = 2, n
   xmag = abs(sx(i))
   if (xmag .ge. smin) go to 30
   idamin = i
   smin = xmag
30 continue
return
end
integer function idmin(n, sx, incx)
!
!     find smallest index of minimum element in double precision s
!     isamax =  first i,  i = 1 to n,  to minimize  sx(1-incx+i*incx)
!
implicit none
integer, intent(in) :: n
real(8), intent(in) :: sx(n)
integer, intent(in) :: incx


real(8) :: smin, xmag
integer :: i, ii, ns

idmin = 0
if (n .le. 0) return
idmin = 1
if (n .le. 1) return

if (incx .eq. 1) go to 20
!
!        code for increments not equal to 1.
!
smin = sx(1)
ns = n * incx
ii = 1
    do 10 i = 1, ns, incx
    xmag = sx(i)
    if (xmag .ge. smin) go to 5
    idmin = ii
    smin = xmag
5     ii = ii + 1
10     continue
return
!
!        code for increments equal to 1.
!
20 smin = abs(sx(1))
do 30 i = 2, n
   xmag = sx(i)
   if (xmag .ge. smin) go to 30
   idmin = i
   smin = xmag
30 continue
return
end
#endif
end module mod_hivector
