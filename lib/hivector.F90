module mod_hivector
contains
!*************************************************************************
!                                                                        *
!                    vector routines library supplement             *
!                                                                        *
!*************************************************************************
!                        routines included                               *
!                                                                        *
!  1. matmov   puts a nr x nc matrix a into nr x nc matrix b             *
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
subroutine matmov (a, b, nr, nc, na, nb)
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
use mod_hiblas, only: dscal, dcopy
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
integer            idamax
external           lsame, idamax
!     ..
!     .. External Subroutines ..
external           dgemm, dgemv, dswap
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
