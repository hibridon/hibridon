**************************************************************************
*                                                                        *
*                    vector routines library supplement             *
*                                                                        *
**************************************************************************
*                        routines included                               *
*                                                                        *
*  1. matmov   puts a nr x nc matrix a into nr x nc matrix b             *
*  2. maxmgv   finds largest value in a vector                           *
*  3. vsmul    multiplies the elements of a vector by a scalar           *
*  4. vmul     multiplies the elements of two vectors                    *
*  4a.dsum     sum of elements of a vector
*  5. blas     basis linear algebra routines from lapack
*     lsame,ilaenv,xerbla
*  6. blas extensions from lapack
*     dlaev2, dlasyf
*     isamax, saxpy, scopy, sdot, sscal, sswap, srot (grot)
*  7. fzero (vector zero) (no longer part of code)
*  8. blas extensions from ibm essl
*     idamin, idmin
*                                                                        *
**************************************************************************
cstart unix cray
      subroutine matmov (a, b, nr, nc, na, nb)
*  subroutine to put nr x nc matrix a into nr x nc matrix b
*  author:  millard alexander
*  current revision date: 24-sept-87
c ------------------------------------------------------------------
*  variables in call list
*    a,b:     input matrices, stored as one-dimensional arrays
*    nr:      actual row dimension of matrix a
*    nc:      actual column dimension of matrix a
*    na:      maximum row dimension of matrix a
*    nb:      maximum row dimension of matrix b
c ------------------------------------------------------------------
*  the two matrices are treated as vectors here, with column-by-column
*  ordering assumed
*  the coding uses the blas routine scopy
c ------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer ia, ib, j, na, nb, nr, nc
      dimension a(1), b(1)
      ia = 0
      ib = 0
*  ia and ib point to one position before the top of the jth
*  column of each matrix
      do 20  j = 1, nc
cend
cstart unix-convex
c;      call scopy (nr, a(ia+1), 1, b(ib+1), 1)
cend
cstart unix cray .and. .not. unix-convex
        call dcopy (nr, a(ia+1), 1, b(ib+1), 1)
cend
cstart unix cray
        ia = ia + na
        ib = ib + nb
20    continue
      return
      end
*  -----------------------------------------------------------------------
      subroutine maxmgv (a, na, c, nc, n)
*  subroutine to scan a vector for its maximum magnitude (absolute value)
*  element
*  current revision date: 24-sept-87
*  -----------------------------------------------------------------------
*  variables in call list:
*  a:   floating point input vector
*  na:  integer element step for a
*  c:   floating point output scalar: on return contains value of
*       maximum magnitude (absolute value) element
*  nc:  integer index of maximum magnitude element
*  n:   integer element count
*  subroutines called:
*  isamax: blas routine to find index of maximum magnitude (absolute value)
*          element
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension a(1)
cend
cstart unix-convex
c;      nc = ( isamax (n, a, na) - 1) * na + 1
cend
cstart unix cray .and. .not. unix-convex
      nc = ( idamax (n, a, na) - 1) * na + 1
cend
cstart unix cray
      c = abs( a(nc) )
      return
      end
*  -----------------------------------------------------------------------
      subroutine vsmul (a, na, b, c, nc, n)
*  subroutine to multiply the elements of a vector by a scalar
*  current revision date: 24-sept-87
*  -----------------------------------------------------------------------
*  variables in call list:
*  a:   floating point input vector
*  na:  integer element step for a
*  b:   floating point input scalar
*  c:   floating point output vector
*  nc:  integer element step for c
*  n:   integer element count
*  c(m) = a(m) * b for m=1 to n
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer n, na, nc
      dimension a(1), c(1)
*  first copy vector a into vector c
*  then multiply by scalar constant
cend
cstart unix-convex
c;      call scopy (n, a, na, c, nc)
c;      call sscal (n, b, c, nc)
cend
cstart unix  cray .and. .not. unix-convex
      call dcopy (n, a, na, c, nc)
      call dscal (n, b, c, nc)
cend
cstart unix cray
      return
      end
*  -----------------------------------------------------------------------
      subroutine vmul (a, na, b, nb, c, nc, n)
*  subroutine to multiply the elements of two vectors
*  current revision date: 23-sept-87
*  -----------------------------------------------------------------------
*  variables in call list:
*  a:   floating point input vector
*  na:  integer element step for a
*  b:   floating point input vector
*  nb:  integer element step for b
*  c:   floating point output vector
*  nc:  integer element step for c
*  n:   integer element count
*  c(m) = a(m) * b(m) for m=1 to n
*  -----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer i, inda, indb, indc, n, na, nb, nc
      dimension a(1), b(1), c(1)
      inda = 1
      indb = 1
      indc = 1
      do 4  i = 1, n
        c(indc) = b(indb) * a(inda)
        inda = inda + na
        indb = indb + nb
        indc = indc + nc
4     continue
      return
      end
*  -----------------------------------------------------------------------
      subroutine vadd (ic,a, na, b, nb, n)
*  subroutine to add or subtract the elements of two vectors
*  current revision date: 6-dec-1991
*  -----------------------------------------------------------------------
*  variables in call list:
*  ic:  factor
*  a:   floating point input vector
*  na:  integer element step for a
*  b:   floating point input vector
*  nb:  integer element step for b
*  n:   integer element count
*  a(m) = a(m) + ic*b(m) for m=1 to n
*  -----------------------------------------------------------------------
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
4       continue
      else
        do 5  i = 1, n
          a(inda) = a(inda) - b(indb)
          inda = inda + na
          indb = indb + nb
5       continue
      endif
      return
      end
      double precision function dsum(n,dx,incx)
c
c     returns sum of double precision dx
c     dasum = sum from 0 to n-1 of dx(1+i*incx)
c     adapted from blas dasum by mha  4-apr-1996
c
      double precision dx(1)
      dsum = 0.d0
      if(n.le.0)return
      if(incx.eq.1)goto 20
c
c        code for increments not equal to 1
c
      ns = n*incx
          do 10 i=1,ns,incx
          dsum = dsum + dx(i)
   10     continue
      return
c
c        code for increments equal to 1.
c
c
c        clean-up loop so remaining vector length is a multiple of 6.
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         dsum = dsum + dx(i)
   30 continue
      if( n .lt. 6 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,6
         dsum = dsum + dx(i) + dx(i+1) + dx(i+2)
     1   + dx(i+3) + dx(i+4) + dx(i+5)
   50 continue
      return
      end
cend
cstart  unix-sun
c;      logical          function lsame( ca, cb )
c;*
c;*  -- LAPACK auxiliary routine (version 1.0) --
c;*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c;*     Courant Institute, Argonne National Lab, and Rice University
c;*     February 29, 1992
c;*
c;*     .. Scalar Arguments ..
c;      character          ca, cb
c;*     ..
c;*
c;*  Purpose
c;*  =======
c;*
c;*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
c;*  case.
c;*
c;*  Arguments
c;*  =========
c;*
c;*  CA      (input) CHARACTER*1
c;*  CB      (input) CHARACTER*1
c;*          CA and CB specify the single characters to be compared.
c;*
c;*     .. Intrinsic Functions ..
c;      intrinsic          ichar
c;*     ..
c;*     .. Local Scalars ..
c;      integer            inta, intb, zcode
c;*     ..
c;*     .. Executable Statements ..
c;*
c;*     Test if the characters are equal
c;*
c;      lsame = ca.eq.cb
c;      if( lsame )
c;     $   return
c;*
c;*     Now test for equivalence if both characters are alphabetic.
c;*
c;      zcode = ichar( 'Z' )
c;*
c;*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
c;*     machines, on which ICHAR returns a value with bit 8 set.
c;*     ICHAR('A') on Prime machines returns 193 which is the same as
c;*     ICHAR('A') on an EBCDIC machine.
c;*
c;      inta = ichar( ca )
c;      intb = ichar( cb )
c;*
c;      if( zcode.eq.90 .or. zcode.eq.122 ) then
c;*
c;*        ASCII is assumed - ZCODE is the ASCII code of either lower or
c;*        upper case 'Z'.
c;*
c;         if( inta.ge.97 .and. inta.le.122 ) inta = inta - 32
c;         if( intb.ge.97 .and. intb.le.122 ) intb = intb - 32
c;*
c;      else if( zcode.eq.233 .or. zcode.eq.169 ) then
c;*
c;*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
c;*        upper case 'Z'.
c;*
c;         if( inta.ge.129 .and. inta.le.137 .or.
c;     $       inta.ge.145 .and. inta.le.153 .or.
c;     $       inta.ge.162 .and. inta.le.169 ) inta = inta + 64
c;         if( intb.ge.129 .and. intb.le.137 .or.
c;     $       intb.ge.145 .and. intb.le.153 .or.
c;     $       intb.ge.162 .and. intb.le.169 ) intb = intb + 64
c;*
c;      else if( zcode.eq.218 .or. zcode.eq.250 ) then
c;*
c;*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
c;*        plus 128 of either lower or upper case 'Z'.
c;*
c;         if( inta.ge.225 .and. inta.le.250 ) inta = inta - 32
c;         if( intb.ge.225 .and. intb.le.250 ) intb = intb - 32
c;      end if
c;      lsame = inta.eq.intb
c;*
c;*     RETURN
c;*
c;*     End of LSAME
c;*
c;      end
c;      subroutine xerbla( srname, info )
c;*
c;*  -- LAPACK auxiliary routine (preliminary version) --
c;*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c;*     Courant Institute, Argonne National Lab, and Rice University
c;*     February 29, 1992
c;*
c;*     .. Scalar Arguments ..
c;      character*6        srname
c;      integer            info
c;*     ..
c;*
c;*  Purpose
c;*  =======
c;*
c;*  XERBLA  is an error handler for the LAPACK routines.
c;*  It is called by an LAPACK routine if an input parameter has an
c;*  invalid value.  A message is printed and execution stops.
c;*
c;*  Installers may consider modifying the STOP statement in order to
c;*  call system-specific exception-handling facilities.
c;*
c;*  Arguments
c;*  =========
c;*
c;*  SRNAME  (input) CHARACTER*6
c;*          The name of the routine which called XERBLA.
c;*
c;*  INFO    (input) INTEGER
c;*          The position of the invalid parameter in the parameter list
c;*          of the calling routine.
c;*
c;*
c;      write( *, fmt = 9999 )srname, info
c;*
c;      stop
c;*
c; 9999 format( ' ** On entry to ', a6, ' parameter number ', i2, ' had ',
c;     $      'an illegal value' )
c;*
c;*     End of XERBLA
c;*
c;      end
cend
cstart unix-hp unix-ibm unix-aix unix-iris unix-sun
c;      integer          function ilaenv( ispec, name, opts, n1, n2, n3,
c;     $                 n4 )
c;*
c;*  -- LAPACK auxiliary routine (preliminary version) --
c;*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c;*     Courant Institute, Argonne National Lab, and Rice University
c;*     February 20, 1992
c;*
c;*     .. Scalar Arguments ..
c;      character*( * )    name, opts
c;      integer            ispec, n1, n2, n3, n4
c;*     ..
c;*
c;*  Purpose
c;*  =======
c;*
c;*  ILAENV is called from the LAPACK routines to choose problem-dependent
c;*  parameters for the local environment.  See ISPEC for a description of
c;*  the parameters.
c;*
c;*  This version provides a set of parameters which should give good,
c;*  but not optimal, performance on many of the currently available
c;*  computers.  Users are encouraged to modify this subroutine to set
c;*  the tuning parameters for their particular machine using the option
c;*  and problem size information in the arguments.
c;*
c;*  This routine will not function correctly if it is converted to all
c;*  lower case.  Converting it to all upper case is allowed.
c;*
c;*  Arguments
c;*  =========
c;*
c;*  ISPEC   (input) INTEGER
c;*          Specifies the parameter to be returned as the value of
c;*          ILAENV.
c;*          = 1: the optimal blocksize; if this value is 1, an unblocked
c;*               algorithm will give the best performance.
c;*          = 2: the minimum block size for which the block routine
c;*               should be used; if the usable block size is less than
c;*               this value, an unblocked routine should be used.
c;*          = 3: the crossover point (in a block routine, for N less
c;*               than this value, an unblocked routine should be used)
c;*          = 4: the number of shifts, used in the nonsymmetric
c;*               eigenvalue routines
c;*          = 5: the minimum column dimension for blocking to be used;
c;*               rectangular blocks must have dimension at least k by m,
c;*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
c;*          = 6: the crossover point for the SVD (when reducing an m by n
c;*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
c;*               this value, a QR factorization is used first to reduce
c;*               the matrix to a triangular form.)
c;*          = 7: the number of processors
c;*          = 8: the crossover point for the multishift QR and QZ methods
c;*               for nonsymmetric eigenvalue problems.
c;*
c;*  NAME    (input) CHARACTER*(*)
c;*          The name of the calling subroutine, in either upper case or
c;*          lower case.
c;*
c;*  OPTS    (input) CHARACTER*(*)
c;*          The character options to the subroutine NAME, concatenated
c;*          into a single character string.  For example, UPLO = 'U',
c;*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
c;*          be specified as OPTS = 'UTN'.
c;*
c;*  N1      (input) INTEGER
c;*  N2      (input) INTEGER
c;*  N3      (input) INTEGER
c;*  N4      (input) INTEGER
c;*          Problem dimensions for the subroutine NAME; these may not all
c;*          be required.
c;*
c;* (ILAENV) (output) INTEGER
c;*          >= 0: the value of the parameter specified by ISPEC
c;*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
c;*
c;*  Further Details
c;*  ===============
c;*
c;*  The following conventions have been used when calling ILAENV from the
c;*  LAPACK routines:
c;*  1)  OPTS is a concatenation of all of the character options to
c;*      subroutine NAME, in the same order that they appear in the
c;*      argument list for NAME, even if they are not used in determining
c;*      the value of the parameter specified by ISPEC.
c;*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
c;*      that they appear in the argument list for NAME.  N1 is used
c;*      first, N2 second, and so on, and unused problem dimensions are
c;*      passed a value of -1.
c;*  3)  The parameter value returned by ILAENV is checked for validity in
c;*      the calling subroutine.  For example, ILAENV is used to retrieve
c;*      the optimal blocksize for STRTRI as follows:
c;*
c;*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
c;*      IF( NB.LE.1 ) NB = MAX( 1, N )
c;*
c;*  =====================================================================
c;*
c;*     .. Local Scalars ..
c;      logical            cname, sname
c;      character*1        c1
c;      character*2        c2, c4
c;      character*3        c3
c;      character*6        subnam
c;      integer            i, ic, iz, nb, nbmin, nx
c;*     ..
c;*     .. Intrinsic Functions ..
c;      intrinsic          char, ichar, int, min, real
c;*     ..
c;*     .. Executable Statements ..
c;*
c;      go to ( 100, 100, 100, 400, 500, 600, 700, 800 ) ispec
c;*
c;*     Invalid value for ISPEC
c;*
c;      ilaenv = -1
c;      return
c;*
c;  100 continue
c;*
c;*     Convert NAME to upper case if the first character is lower case.
c;*
c;      ilaenv = 1
c;      subnam = name
c;      ic = ichar( subnam( 1:1 ) )
c;      iz = ichar( 'Z' )
c;      if( iz.eq.90 .or. iz.eq.122 ) then
c;*
c;*        ASCII character set
c;*
c;         if( ic.ge.97 .and. ic.le.122 ) then
c;            subnam( 1:1 ) = char( ic-32 )
c;            do 10 i = 2, 6
c;               ic = ichar( subnam( i:i ) )
c;               if( ic.ge.97 .and. ic.le.122 )
c;     $            subnam( i:i ) = char( ic-32 )
c;   10       continue
c;         end if
c;*
c;      else if( iz.eq.233 .or. iz.eq.169 ) then
c;*
c;*        EBCDIC character set
c;*
c;         if( ( ic.ge.129 .and. ic.le.137 ) .or.
c;     $       ( ic.ge.145 .and. ic.le.153 ) .or.
c;     $       ( ic.ge.162 .and. ic.le.169 ) ) then
c;            subnam( 1:1 ) = char( ic+64 )
c;            do 20 i = 2, 6
c;               ic = ichar( subnam( i:i ) )
c;               if( ( ic.ge.129 .and. ic.le.137 ) .or.
c;     $             ( ic.ge.145 .and. ic.le.153 ) .or.
c;     $             ( ic.ge.162 .and. ic.le.169 ) )
c;     $            subnam( i:i ) = char( ic+64 )
c;   20       continue
c;         end if
c;*
c;      else if( iz.eq.218 .or. iz.eq.250 ) then
c;*
c;*        Prime machines:  ASCII+128
c;*
c;         if( ic.ge.225 .and. ic.le.250 ) then
c;            subnam( 1:1 ) = char( ic-32 )
c;            do 30 i = 2, 6
c;               ic = ichar( subnam( i:i ) )
c;               if( ic.ge.225 .and. ic.le.250 )
c;     $            subnam( i:i ) = char( ic-32 )
c;   30       continue
c;         end if
c;      end if
c;*
c;      c1 = subnam( 1:1 )
c;      sname = c1.eq.'S' .or. c1.eq.'D'
c;      cname = c1.eq.'C' .or. c1.eq.'Z'
c;      if( .not.( cname .or. sname ) )
c;     $   return
c;      c2 = subnam( 2:3 )
c;      c3 = subnam( 4:6 )
c;      c4 = c3( 2:3 )
c;*
c;      go to ( 110, 200, 300 ) ispec
c;*
c;  110 continue
c;*
c;*     ISPEC = 1:  block size
c;*
c;*     In these examples, separate code is provided for setting NB for
c;*     real and complex.  We assume that NB will take the same value in
c;*     single or double precision.
c;*
c;      nb = 1
c;*
c;      if( c2.eq.'GE' ) then
c;         if( c3.eq.'TRF' ) then
c;            if( sname ) then
c;               nb = 64
c;            else
c;               nb = 64
c;            end if
c;         else if( c3.eq.'QRF' .or. c3.eq.'RQF' .or. c3.eq.'LQF' .or.
c;     $            c3.eq.'QLF' ) then
c;            if( sname ) then
c;               nb = 32
c;            else
c;               nb = 32
c;            end if
c;         else if( c3.eq.'HRD' ) then
c;            if( sname ) then
c;               nb = 32
c;            else
c;               nb = 32
c;            end if
c;         else if( c3.eq.'BRD' ) then
c;            if( sname ) then
c;               nb = 32
c;            else
c;               nb = 32
c;            end if
c;         else if( c3.eq.'TRI' ) then
c;            if( sname ) then
c;               nb = 64
c;            else
c;               nb = 64
c;            end if
c;         end if
c;      else if( c2.eq.'PO' ) then
c;         if( c3.eq.'TRF' ) then
c;            if( sname ) then
c;               nb = 64
c;            else
c;               nb = 64
c;            end if
c;         end if
c;      else if( c2.eq.'SY' ) then
c;         if( c3.eq.'TRF' ) then
c;            if( sname ) then
c;               nb = 64
c;            else
c;               nb = 64
c;            end if
c;         else if( sname .and. c3.eq.'TRD' ) then
c;            nb = 1
c;         else if( sname .and. c3.eq.'GST' ) then
c;            nb = 64
c;         end if
c;      else if( cname .and. c2.eq.'HE' ) then
c;         if( c3.eq.'TRF' ) then
c;            nb = 64
c;         else if( c3.eq.'TRD' ) then
c;            nb = 1
c;         else if( c3.eq.'GST' ) then
c;            nb = 64
c;         end if
c;      else if( sname .and. c2.eq.'OR' ) then
c;         if( c3( 1:1 ).eq.'G' ) then
c;            if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or.
c;     $          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or.
c;     $          c4.eq.'BR' ) then
c;               nb = 32
c;            end if
c;         else if( c3( 1:1 ).eq.'M' ) then
c;            if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or.
c;     $          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or.
c;     $          c4.eq.'BR' ) then
c;               nb = 32
c;            end if
c;         end if
c;      else if( cname .and. c2.eq.'UN' ) then
c;         if( c3( 1:1 ).eq.'G' ) then
c;            if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or.
c;     $          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or.
c;     $          c4.eq.'BR' ) then
c;               nb = 32
c;            end if
c;         else if( c3( 1:1 ).eq.'M' ) then
c;            if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or.
c;     $          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or.
c;     $          c4.eq.'BR' ) then
c;               nb = 32
c;            end if
c;         end if
c;      else if( c2.eq.'GB' ) then
c;         if( c3.eq.'TRF' ) then
c;            if( sname ) then
c;               if( n4.le.64 ) then
c;                  nb = 1
c;               else
c;                  nb = 32
c;               end if
c;            else
c;               if( n4.le.64 ) then
c;                  nb = 1
c;               else
c;                  nb = 32
c;               end if
c;            end if
c;         end if
c;      else if( c2.eq.'PB' ) then
c;         if( c3.eq.'TRF' ) then
c;            if( sname ) then
c;               if( n2.le.64 ) then
c;                  nb = 1
c;               else
c;                  nb = 32
c;               end if
c;            else
c;               if( n2.le.64 ) then
c;                  nb = 1
c;               else
c;                  nb = 32
c;               end if
c;            end if
c;         end if
c;      else if( c2.eq.'TR' ) then
c;         if( c3.eq.'TRI' ) then
c;            if( sname ) then
c;               nb = 64
c;            else
c;               nb = 64
c;            end if
c;         end if
c;      else if( c2.eq.'LA' ) then
c;         if( c3.eq.'UUM' ) then
c;            if( sname ) then
c;               nb = 64
c;            else
c;               nb = 64
c;            end if
c;         end if
c;      else if( sname .and. c2.eq.'ST' ) then
c;         if( c3.eq.'EBZ' ) then
c;            nb = 1
c;         end if
c;      end if
c;      ilaenv = nb
c;      return
c;*
c;  200 continue
c;*
c;*     ISPEC = 2:  minimum block size
c;*
c;      nbmin = 2
c;      if( c2.eq.'GE' ) then
c;         if( c3.eq.'QRF' .or. c3.eq.'RQF' .or. c3.eq.'LQF' .or.
c;     $       c3.eq.'QLF' ) then
c;            if( sname ) then
c;               nbmin = 2
c;            else
c;               nbmin = 2
c;            end if
c;         else if( c3.eq.'HRD' ) then
c;            if( sname ) then
c;               nbmin = 2
c;            else
c;               nbmin = 2
c;            end if
c;         else if( c3.eq.'BRD' ) then
c;            if( sname ) then
c;               nbmin = 2
c;            else
c;               nbmin = 2
c;            end if
c;         else if( c3.eq.'TRI' ) then
c;            if( sname ) then
c;               nbmin = 2
c;            else
c;               nbmin = 2
c;            end if
c;         end if
c;      else if( c2.eq.'SY' ) then
c;         if( c3.eq.'TRF' ) then
c;            if( sname ) then
c;               nbmin = 2
c;            else
c;               nbmin = 2
c;            end if
c;         else if( sname .and. c3.eq.'TRD' ) then
c;            nbmin = 2
c;         end if
c;      else if( cname .and. c2.eq.'HE' ) then
c;         if( c3.eq.'TRD' ) then
c;            nbmin = 2
c;         end if
c;      else if( sname .and. c2.eq.'OR' ) then
c;         if( c3( 1:1 ).eq.'G' ) then
c;            if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or.
c;     $          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or.
c;     $          c4.eq.'BR' ) then
c;               nbmin = 2
c;            end if
c;         else if( c3( 1:1 ).eq.'M' ) then
c;            if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or.
c;     $          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or.
c;     $          c4.eq.'BR' ) then
c;               nbmin = 2
c;            end if
c;         end if
c;      else if( cname .and. c2.eq.'UN' ) then
c;         if( c3( 1:1 ).eq.'G' ) then
c;            if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or.
c;     $          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or.
c;     $          c4.eq.'BR' ) then
c;               nbmin = 2
c;            end if
c;         else if( c3( 1:1 ).eq.'M' ) then
c;            if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or.
c;     $          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or.
c;     $          c4.eq.'BR' ) then
c;               nbmin = 2
c;            end if
c;         end if
c;      end if
c;      ilaenv = nbmin
c;      return
c;*
c;  300 continue
c;*
c;*     ISPEC = 3:  crossover point
c;*
c;      nx = 0
c;      if( c2.eq.'GE' ) then
c;         if( c3.eq.'QRF' .or. c3.eq.'RQF' .or. c3.eq.'LQF' .or.
c;     $       c3.eq.'QLF' ) then
c;            if( sname ) then
c;               nx = 128
c;            else
c;               nx = 128
c;            end if
c;         else if( c3.eq.'HRD' ) then
c;            if( sname ) then
c;               nx = 128
c;            else
c;               nx = 128
c;            end if
c;         else if( c3.eq.'BRD' ) then
c;            if( sname ) then
c;               nx = 128
c;            else
c;               nx = 128
c;            end if
c;         end if
c;      else if( c2.eq.'SY' ) then
c;         if( sname .and. c3.eq.'TRD' ) then
c;            nx = 1
c;         end if
c;      else if( cname .and. c2.eq.'HE' ) then
c;         if( c3.eq.'TRD' ) then
c;            nx = 1
c;         end if
c;      else if( sname .and. c2.eq.'OR' ) then
c;         if( c3( 1:1 ).eq.'G' ) then
c;            if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or.
c;     $          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or.
c;     $          c4.eq.'BR' ) then
c;               nx = 128
c;            end if
c;         end if
c;      else if( cname .and. c2.eq.'UN' ) then
c;         if( c3( 1:1 ).eq.'G' ) then
c;            if( c4.eq.'QR' .or. c4.eq.'RQ' .or. c4.eq.'LQ' .or.
c;     $          c4.eq.'QL' .or. c4.eq.'HR' .or. c4.eq.'TR' .or.
c;     $          c4.eq.'BR' ) then
c;               nx = 128
c;            end if
c;         end if
c;      end if
c;      ilaenv = nx
c;      return
c;*
c;  400 continue
c;*
c;*     ISPEC = 4:  number of shifts (used by xHSEQR)
c;*
c;      ilaenv = 6
c;      return
c;*
c;  500 continue
c;*
c;*     ISPEC = 5:  minimum column dimension (not used)
c;*
c;      ilaenv = 2
c;      return
c;*
c;  600 continue
c;*
c;*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
c;*
c;      ilaenv = int( real( min( n1, n2 ) )*1.6e0 )
c;      return
c;*
c;  700 continue
c;*
c;*     ISPEC = 7:  number of processors (not used)
c;*
c;      ilaenv = 1
c;      return
c;*
c;  800 continue
c;*
c;*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
c;*
c;      ilaenv = 50
c;      return
c;*
c;*     End of ILAENV
c;*
c;      end
c;      subroutine dlaev2( a, b, c, rt1, rt2, cs1, sn1 )
c;*
c;*  -- LAPACK auxiliary routine (version 1.0b) --
c;*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c;*     Courant Institute, Argonne National Lab, and Rice University
c;*     October 31, 1992
c;*
c;*     .. Scalar Arguments ..
c;      double precision   a, b, c, cs1, rt1, rt2, sn1
c;*     ..
c;*
c;*  Purpose
c;*  =======
c;*
c;*  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
c;*     [  A   B  ]
c;*     [  B   C  ].
c;*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
c;*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
c;*  eigenvector for RT1, giving the decomposition
c;*
c;*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
c;*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
c;*
c;*  Arguments
c;*  =========
c;*
c;*  A       (input) DOUBLE PRECISION
c;*          The (1,1) entry of the 2-by-2 matrix.
c;*
c;*  B       (input) DOUBLE PRECISION
c;*          The (1,2) entry and the conjugate of the (2,1) entry of the
c;*          2-by-2 matrix.
c;*
c;*  C       (input) DOUBLE PRECISION
c;*          The (2,2) entry of the 2-by-2 matrix.
c;*
c;*  RT1     (output) DOUBLE PRECISION
c;*          The eigenvalue of larger absolute value.
c;*
c;*  RT2     (output) DOUBLE PRECISION
c;*          The eigenvalue of smaller absolute value.
c;*
c;*  CS1     (output) DOUBLE PRECISION
c;*  SN1     (output) DOUBLE PRECISION
c;*          The vector (CS1, SN1) is a unit right eigenvector for RT1.
c;*
c;*  Further Details
c;*  ===============
c;*
c;*  RT1 is accurate to a few ulps barring over/underflow.
c;*
c;*  RT2 may be inaccurate if there is massive cancellation in the
c;*  determinant A*C-B*B; higher precision or correctly rounded or
c;*  correctly truncated arithmetic would be needed to compute RT2
c;*  accurately in all cases.
c;*
c;*  CS1 and SN1 are accurate to a few ulps barring over/underflow.
c;*
c;*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
c;*  Underflow is harmless if the input data is 0 or exceeds
c;*     underflow_threshold / macheps.
c;*
c;* =====================================================================
c;*
c;*     .. Parameters ..
c;      double precision   one
c;      parameter          ( one = 1.0d0 )
c;      double precision   two
c;      parameter          ( two = 2.0d0 )
c;      double precision   zero
c;      parameter          ( zero = 0.0d0 )
c;      double precision   half
c;      parameter          ( half = 0.5d0 )
c;*     ..
c;*     .. Local Scalars ..
c;      integer            sgn1, sgn2
c;      double precision   ab, acmn, acmx, acs, adf, cs, ct, df, rt, sm,
c;     $                   tb, tn
c;*     ..
c;*     .. Intrinsic Functions ..
c;      intrinsic          abs, sqrt
c;*     ..
c;*     .. Executable Statements ..
c;*
c;*     Compute the eigenvalues
c;*
c;      sm = a + c
c;      df = a - c
c;      adf = abs( df )
c;      tb = b + b
c;      ab = abs( tb )
c;      if( abs( a ).gt.abs( c ) ) then
c;         acmx = a
c;         acmn = c
c;      else
c;         acmx = c
c;         acmn = a
c;      end if
c;      if( adf.gt.ab ) then
c;         rt = adf*sqrt( one+( ab / adf )**2 )
c;      else if( adf.lt.ab ) then
c;         rt = ab*sqrt( one+( adf / ab )**2 )
c;      else
c;*
c;*        Includes case AB=ADF=0
c;*
c;         rt = ab*sqrt( two )
c;      end if
c;      if( sm.lt.zero ) then
c;         rt1 = half*( sm-rt )
c;         sgn1 = -1
c;*
c;*        Order of execution important.
c;*        To get fully accurate smaller eigenvalue,
c;*        next line needs to be executed in higher precision.
c;*
c;         rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
c;      else if( sm.gt.zero ) then
c;         rt1 = half*( sm+rt )
c;         sgn1 = 1
c;*
c;*        Order of execution important.
c;*        To get fully accurate smaller eigenvalue,
c;*        next line needs to be executed in higher precision.
c;*
c;         rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
c;      else
c;*
c;*        Includes case RT1 = RT2 = 0
c;*
c;         rt1 = half*rt
c;         rt2 = -half*rt
c;         sgn1 = 1
c;      end if
c;*
c;*     Compute the eigenvector
c;*
c;      if( df.ge.zero ) then
c;         cs = df + rt
c;         sgn2 = 1
c;      else
c;         cs = df - rt
c;         sgn2 = -1
c;      end if
c;      acs = abs( cs )
c;      if( acs.gt.ab ) then
c;         ct = -tb / cs
c;         sn1 = one / sqrt( one+ct*ct )
c;         cs1 = ct*sn1
c;      else
c;         if( ab.eq.zero ) then
c;            cs1 = one
c;            sn1 = zero
c;         else
c;            tn = -cs / tb
c;            cs1 = one / sqrt( one+tn*tn )
c;            sn1 = tn*cs1
c;         end if
c;      end if
c;      if( sgn1.eq.sgn2 ) then
c;         tn = cs1
c;         cs1 = -sn1
c;         sn1 = tn
c;      end if
c;      return
c;*
c;*     End of DLAEV2
c;*
c;      end
c;      subroutine dlasyf( uplo, n, nb, kb, a, lda, ipiv, w, ldw, info )
c;*
c;*  -- LAPACK routine (version 1.0) --
c;*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c;*     Courant Institute, Argonne National Lab, and Rice University
c;*     February 29, 1992
c;*
c;*     .. Scalar Arguments ..
c;      character          uplo
c;      integer            info, kb, lda, ldw, n, nb
c;*     ..
c;*     .. Array Arguments ..
c;      integer            ipiv( * )
c;      double precision   a( lda, * ), w( ldw, * )
c;*     ..
c;*
c;*  Purpose
c;*  =======
c;*
c;*  DLASYF computes a partial factorization of a real symmetric matrix A
c;*  using the Bunch-Kaufman diagonal pivoting method. The partial
c;*  factorization has the form:
c;*
c;*  A  =  ( I  U12 ) ( A11  0  ) (  I    0   )  if UPLO = 'U', or:
c;*        ( 0  U22 ) (  0   D  ) ( U12' U22' )
c;*
c;*  A  =  ( L11  0 ) (  D   0  ) ( L11' L21' )  if UPLO = 'L'
c;*        ( L21  I ) (  0  A22 ) (  0    I   )
c;*
c;*  where the order of D is at most NB. The actual order is returned in
c;*  the argument KB, and is either NB or NB-1, or N if N <= NB.
c;*
c;*  DLASYF is an auxiliary routine called by DSYTRF. It uses blocked code
c;*  (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
c;*  A22 (if UPLO = 'L').
c;*
c;*  Arguments
c;*  =========
c;*
c;*  UPLO    (input) CHARACTER*1
c;*          Specifies whether the upper or lower triangular part of the
c;*          symmetric matrix A is stored:
c;*          = 'U':  Upper triangular
c;*          = 'L':  Lower triangular
c;*
c;*  N       (input) INTEGER
c;*          The order of the matrix A.  N >= 0.
c;*
c;*  NB      (input) INTEGER
c;*          The maximum number of columns of the matrix A that should be
c;*          factored.  NB should be at least 2 to allow for 2-by-2 pivot
c;*          blocks.
c;*
c;*  KB      (output) INTEGER
c;*          The number of columns of A that were actually factored.
c;*          KB is either NB-1 or NB, or N if N <= NB.
c;*
c;*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
c;*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
c;*          n-by-n upper triangular part of A contains the upper
c;*          triangular part of the matrix A, and the strictly lower
c;*          triangular part of A is not referenced.  If UPLO = 'L', the
c;*          leading n-by-n lower triangular part of A contains the lower
c;*          triangular part of the matrix A, and the strictly upper
c;*          triangular part of A is not referenced.
c;*          On exit, A contains details of the partial factorization.
c;*
c;*  LDA     (input) INTEGER
c;*          The leading dimension of the array A.  LDA >= max(1,N).
c;*
c;*  IPIV    (output) INTEGER array, dimension (N)
c;*          Details of the interchanges and the block structure of D.
c;*          If UPLO = 'U', only the last KB elements of IPIV are set;
c;*          if UPLO = 'L', only the first KB elements are set.
c;*
c;*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
c;*          interchanged and D(k,k) is a 1-by-1 diagonal block.
c;*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
c;*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
c;*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
c;*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
c;*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
c;*
c;*  W       (workspace) DOUBLE PRECISION array, dimension (LDW,NB)
c;*
c;*  LDW     (input) INTEGER
c;*          The leading dimension of the array W.  LDW >= max(1,N).
c;*
c;*  INFO    (output) INTEGER
c;*          = 0: successful exit
c;*          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
c;*               has been completed, but the block diagonal matrix D is
c;*               exactly singular.
c;*
c;*  =====================================================================
c;*
c;*     .. Parameters ..
c;      double precision   zero, one
c;      parameter          ( zero = 0.0d+0, one = 1.0d+0 )
c;      double precision   eight, sevten
c;      parameter          ( eight = 8.0d+0, sevten = 17.0d+0 )
c;*     ..
c;*     .. Local Scalars ..
c;      integer            imax, j, jb, jj, jmax, jp, k, kk, kkw, kp,
c;     $                   kstep, kw
c;      double precision   absakk, alpha, colmax, d11, d21, d22, r1,
c;     $                   rowmax, t
c;*     ..
c;*     .. External Functions ..
c;      logical            lsame
c;      integer            idamax
c;      external           lsame, idamax
c;*     ..
c;*     .. External Subroutines ..
c;      external           dcopy, dgemm, dgemv, dscal, dswap
c;*     ..
c;*     .. Intrinsic Functions ..
c;      intrinsic          abs, max, min, sqrt
c;*     ..
c;*     .. Executable Statements ..
c;*
c;      info = 0
c;*
c;*     Initialize ALPHA for use in choosing pivot block size.
c;*
c;      alpha = ( one+sqrt( sevten ) ) / eight
c;*
c;      if( lsame( uplo, 'U' ) ) then
c;*
c;*        Factorize the trailing columns of A using the upper triangle
c;*        of A and working backwards, and compute the matrix W = U12*D
c;*        for use in updating A11
c;*
c;*        K is the main loop index, decreasing from N in steps of 1 or 2
c;*
c;*        KW is the column of W which corresponds to column K of A
c;*
c;         k = n
c;   10    continue
c;         kw = nb + k - n
c;*
c;*        Exit from loop
c;*
c;         if( ( k.le.n-nb+1 .and. nb.lt.n ) .or. k.lt.1 )
c;     $      go to 30
c;*
c;*        Copy column K of A to column KW of W and update it
c;*
c;         call dcopy( k, a( 1, k ), 1, w( 1, kw ), 1 )
c;         if( k.lt.n )
c;     $      call dgemv( 'No transpose', k, n-k, -one, a( 1, k+1 ), lda,
c;     $                  w( k, kw+1 ), ldw, one, w( 1, kw ), 1 )
c;*
c;         kstep = 1
c;*
c;*        Determine rows and columns to be interchanged and whether
c;*        a 1-by-1 or 2-by-2 pivot block will be used
c;*
c;         absakk = abs( w( k, kw ) )
c;*
c;*        IMAX is the row-index of the largest off-diagonal element in
c;*        column K, and COLMAX is its absolute value
c;*
c;         if( k.gt.1 ) then
c;            imax = idamax( k-1, w( 1, kw ), 1 )
c;            colmax = abs( w( imax, kw ) )
c;         else
c;            colmax = zero
c;         end if
c;*
c;         if( max( absakk, colmax ).eq.zero ) then
c;*
c;*           Column K is zero: set INFO and continue
c;*
c;            if( info.eq.0 )
c;     $         info = k
c;            kp = k
c;         else
c;            if( absakk.ge.alpha*colmax ) then
c;*
c;*              no interchange, use 1-by-1 pivot block
c;*
c;               kp = k
c;            else
c;*
c;*              Copy column IMAX to column KW-1 of W and update it
c;*
c;               call dcopy( imax, a( 1, imax ), 1, w( 1, kw-1 ), 1 )
c;               call dcopy( k-imax, a( imax, imax+1 ), lda,
c;     $                     w( imax+1, kw-1 ), 1 )
c;               if( k.lt.n )
c;     $            call dgemv( 'No transpose', k, n-k, -one, a( 1, k+1 ),
c;     $                        lda, w( imax, kw+1 ), ldw, one,
c;     $                        w( 1, kw-1 ), 1 )
c;*
c;*              JMAX is the column-index of the largest off-diagonal
c;*              element in row IMAX, and ROWMAX is its absolute value
c;*
c;               jmax = imax + idamax( k-imax, w( imax+1, kw-1 ), 1 )
c;               rowmax = abs( w( jmax, kw-1 ) )
c;               if( imax.gt.1 ) then
c;                  jmax = idamax( imax-1, w( 1, kw-1 ), 1 )
c;                  rowmax = max( rowmax, abs( w( jmax, kw-1 ) ) )
c;               end if
c;*
c;               if( absakk.ge.alpha*colmax*( colmax / rowmax ) ) then
c;*
c;*                 no interchange, use 1-by-1 pivot block
c;*
c;                  kp = k
c;               else if( abs( w( imax, kw-1 ) ).ge.alpha*rowmax ) then
c;*
c;*                 interchange rows and columns K and IMAX, use 1-by-1
c;*                 pivot block
c;*
c;                  kp = imax
c;*
c;*                 copy column KW-1 of W to column KW
c;*
c;                  call dcopy( k, w( 1, kw-1 ), 1, w( 1, kw ), 1 )
c;               else
c;*
c;*                 interchange rows and columns K-1 and IMAX, use 2-by-2
c;*                 pivot block
c;*
c;                  kp = imax
c;                  kstep = 2
c;               end if
c;            end if
c;*
c;            kk = k - kstep + 1
c;            kkw = nb + kk - n
c;*
c;*           Updated column KP is already stored in column KKW of W
c;*
c;            if( kp.ne.kk ) then
c;*
c;*              Copy non-updated column KK to column KP
c;*
c;               a( kp, k ) = a( kk, k )
c;               call dcopy( k-1-kp, a( kp+1, kk ), 1, a( kp, kp+1 ),
c;     $                     lda )
c;               call dcopy( kp, a( 1, kk ), 1, a( 1, kp ), 1 )
c;*
c;*              Interchange rows KK and KP in last KK columns of A and W
c;*
c;               call dswap( n-kk+1, a( kk, kk ), lda, a( kp, kk ), lda )
c;               call dswap( n-kk+1, w( kk, kkw ), ldw, w( kp, kkw ),
c;     $                     ldw )
c;            end if
c;*
c;            if( kstep.eq.1 ) then
c;*
c;*              1-by-1 pivot block D(k): column KW of W now holds
c;*
c;*              W(k) = U(k)*D(k)
c;*
c;*              where U(k) is the k-th column of U
c;*
c;*              Store U(k) in column k of A
c;*
c;               call dcopy( k, w( 1, kw ), 1, a( 1, k ), 1 )
c;               r1 = one / a( k, k )
c;               call dscal( k-1, r1, a( 1, k ), 1 )
c;            else
c;*
c;*              2-by-2 pivot block D(k): columns KW and KW-1 of W now
c;*              hold
c;*
c;*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
c;*
c;*              where U(k) and U(k-1) are the k-th and (k-1)-th columns
c;*              of U
c;*
c;               if( k.gt.2 ) then
c;*
c;*                 Store U(k) and U(k-1) in columns k and k-1 of A
c;*
c;                  d21 = w( k-1, kw )
c;                  d11 = w( k, kw ) / d21
c;                  d22 = w( k-1, kw-1 ) / d21
c;                  t = one / ( d11*d22-one )
c;                  d21 = t / d21
c;                  do 20 j = 1, k - 2
c;                     a( j, k-1 ) = d21*( d11*w( j, kw-1 )-w( j, kw ) )
c;                     a( j, k ) = d21*( d22*w( j, kw )-w( j, kw-1 ) )
c;   20             continue
c;               end if
c;*
c;*              Copy D(k) to A
c;*
c;               a( k-1, k-1 ) = w( k-1, kw-1 )
c;               a( k-1, k ) = w( k-1, kw )
c;               a( k, k ) = w( k, kw )
c;            end if
c;         end if
c;*
c;*        Store details of the interchanges in IPIV
c;*
c;         if( kstep.eq.1 ) then
c;            ipiv( k ) = kp
c;         else
c;            ipiv( k ) = -kp
c;            ipiv( k-1 ) = -kp
c;         end if
c;*
c;*        Decrease K and return to the start of the main loop
c;*
c;         k = k - kstep
c;         go to 10
c;*
c;   30    continue
c;*
c;*        Update the upper triangle of A11 (= A(1:k,1:k)) as
c;*
c;*        A11 := A11 - U12*D*U12' = A11 - U12*W'
c;*
c;*        computing blocks of NB columns at a time
c;*
c;         do 50 j = ( ( k-1 ) / nb )*nb + 1, 1, -nb
c;            jb = min( nb, k-j+1 )
c;*
c;*           Update the upper triangle of the diagonal block
c;*
c;            do 40 jj = j, j + jb - 1
c;               call dgemv( 'No transpose', jj-j+1, n-k, -one,
c;     $                     a( j, k+1 ), lda, w( jj, kw+1 ), ldw, one,
c;     $                     a( j, jj ), 1 )
c;   40       continue
c;*
c;*           Update the rectangular superdiagonal block
c;*
c;            call dgemm( 'No transpose', 'Transpose', j-1, jb, n-k, -one,
c;     $                  a( 1, k+1 ), lda, w( j, kw+1 ), ldw, one,
c;     $                  a( 1, j ), lda )
c;   50    continue
c;*
c;*        Put U12 in standard form by partially undoing the interchanges
c;*        in columns k+1:n
c;*
c;         j = k + 1
c;   60    continue
c;         jj = j
c;         jp = ipiv( j )
c;         if( jp.lt.0 ) then
c;            jp = -jp
c;            j = j + 1
c;         end if
c;         j = j + 1
c;         if( jp.ne.jj .and. j.le.n )
c;     $      call dswap( n-j+1, a( jp, j ), lda, a( jj, j ), lda )
c;         if( j.le.n )
c;     $      go to 60
c;*
c;*        Set KB to the number of columns factorized
c;*
c;         kb = n - k
c;*
c;      else
c;*
c;*        Factorize the leading columns of A using the lower triangle
c;*        of A and working forwards, and compute the matrix W = L21*D
c;*        for use in updating A22
c;*
c;*        K is the main loop index, increasing from 1 in steps of 1 or 2
c;*
c;         k = 1
c;   70    continue
c;*
c;*        Exit from loop
c;*
c;         if( ( k.ge.nb .and. nb.lt.n ) .or. k.gt.n )
c;     $      go to 90
c;*
c;*        Copy column K of A to column K of W and update it
c;*
c;         call dcopy( n-k+1, a( k, k ), 1, w( k, k ), 1 )
c;         call dgemv( 'No transpose', n-k+1, k-1, -one, a( k, 1 ), lda,
c;     $               w( k, 1 ), ldw, one, w( k, k ), 1 )
c;*
c;         kstep = 1
c;*
c;*        Determine rows and columns to be interchanged and whether
c;*        a 1-by-1 or 2-by-2 pivot block will be used
c;*
c;         absakk = abs( w( k, k ) )
c;*
c;*        IMAX is the row-index of the largest off-diagonal element in
c;*        column K, and COLMAX is its absolute value
c;*
c;         if( k.lt.n ) then
c;            imax = k + idamax( n-k, w( k+1, k ), 1 )
c;            colmax = abs( w( imax, k ) )
c;         else
c;            colmax = zero
c;         end if
c;*
c;         if( max( absakk, colmax ).eq.zero ) then
c;*
c;*           Column K is zero: set INFO and continue
c;*
c;            if( info.eq.0 )
c;     $         info = k
c;            kp = k
c;         else
c;            if( absakk.ge.alpha*colmax ) then
c;*
c;*              no interchange, use 1-by-1 pivot block
c;*
c;               kp = k
c;            else
c;*
c;*              Copy column IMAX to column K+1 of W and update it
c;*
c;               call dcopy( imax-k, a( imax, k ), lda, w( k, k+1 ), 1 )
c;               call dcopy( n-imax+1, a( imax, imax ), 1, w( imax, k+1 ),
c;     $                     1 )
c;               call dgemv( 'No transpose', n-k+1, k-1, -one, a( k, 1 ),
c;     $                     lda, w( imax, 1 ), ldw, one, w( k, k+1 ), 1 )
c;*
c;*              JMAX is the column-index of the largest off-diagonal
c;*              element in row IMAX, and ROWMAX is its absolute value
c;*
c;               jmax = k - 1 + idamax( imax-k, w( k, k+1 ), 1 )
c;               rowmax = abs( w( jmax, k+1 ) )
c;               if( imax.lt.n ) then
c;                  jmax = imax + idamax( n-imax, w( imax+1, k+1 ), 1 )
c;                  rowmax = max( rowmax, abs( w( jmax, k+1 ) ) )
c;               end if
c;*
c;               if( absakk.ge.alpha*colmax*( colmax / rowmax ) ) then
c;*
c;*                 no interchange, use 1-by-1 pivot block
c;*
c;                  kp = k
c;               else if( abs( w( imax, k+1 ) ).ge.alpha*rowmax ) then
c;*
c;*                 interchange rows and columns K and IMAX, use 1-by-1
c;*                 pivot block
c;*
c;                  kp = imax
c;*
c;*                 copy column K+1 of W to column K
c;*
c;                  call dcopy( n-k+1, w( k, k+1 ), 1, w( k, k ), 1 )
c;               else
c;*
c;*                 interchange rows and columns K+1 and IMAX, use 2-by-2
c;*                 pivot block
c;*
c;                  kp = imax
c;                  kstep = 2
c;               end if
c;            end if
c;*
c;            kk = k + kstep - 1
c;*
c;*           Updated column KP is already stored in column KK of W
c;*
c;            if( kp.ne.kk ) then
c;*
c;*              Copy non-updated column KK to column KP
c;*
c;               a( kp, k ) = a( kk, k )
c;               call dcopy( kp-k-1, a( k+1, kk ), 1, a( kp, k+1 ), lda )
c;               call dcopy( n-kp+1, a( kp, kk ), 1, a( kp, kp ), 1 )
c;*
c;*              Interchange rows KK and KP in first KK columns of A and W
c;*
c;               call dswap( kk, a( kk, 1 ), lda, a( kp, 1 ), lda )
c;               call dswap( kk, w( kk, 1 ), ldw, w( kp, 1 ), ldw )
c;            end if
c;*
c;            if( kstep.eq.1 ) then
c;*
c;*              1-by-1 pivot block D(k): column k of W now holds
c;*
c;*              W(k) = L(k)*D(k)
c;*
c;*              where L(k) is the k-th column of L
c;*
c;*              Store L(k) in column k of A
c;*
c;               call dcopy( n-k+1, w( k, k ), 1, a( k, k ), 1 )
c;               if( k.lt.n ) then
c;                  r1 = one / a( k, k )
c;                  call dscal( n-k, r1, a( k+1, k ), 1 )
c;               end if
c;            else
c;*
c;*              2-by-2 pivot block D(k): columns k and k+1 of W now hold
c;*
c;*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
c;*
c;*              where L(k) and L(k+1) are the k-th and (k+1)-th columns
c;*              of L
c;*
c;               if( k.lt.n-1 ) then
c;*
c;*                 Store L(k) and L(k+1) in columns k and k+1 of A
c;*
c;                  d21 = w( k+1, k )
c;                  d11 = w( k+1, k+1 ) / d21
c;                  d22 = w( k, k ) / d21
c;                  t = one / ( d11*d22-one )
c;                  d21 = t / d21
c;                  do 80 j = k + 2, n
c;                     a( j, k ) = d21*( d11*w( j, k )-w( j, k+1 ) )
c;                     a( j, k+1 ) = d21*( d22*w( j, k+1 )-w( j, k ) )
c;   80             continue
c;               end if
c;*
c;*              Copy D(k) to A
c;*
c;               a( k, k ) = w( k, k )
c;               a( k+1, k ) = w( k+1, k )
c;               a( k+1, k+1 ) = w( k+1, k+1 )
c;            end if
c;         end if
c;*
c;*        Store details of the interchanges in IPIV
c;*
c;         if( kstep.eq.1 ) then
c;            ipiv( k ) = kp
c;         else
c;            ipiv( k ) = -kp
c;            ipiv( k+1 ) = -kp
c;         end if
c;*
c;*        Increase K and return to the start of the main loop
c;*
c;         k = k + kstep
c;         go to 70
c;*
c;   90    continue
c;*
c;*        Update the lower triangle of A22 (= A(k:n,k:n)) as
c;*
c;*        A22 := A22 - L21*D*L21' = A22 - L21*W'
c;*
c;*        computing blocks of NB columns at a time
c;*
c;         do 110 j = k, n, nb
c;            jb = min( nb, n-j+1 )
c;*
c;*           Update the lower triangle of the diagonal block
c;*
c;            do 100 jj = j, j + jb - 1
c;               call dgemv( 'No transpose', j+jb-jj, k-1, -one,
c;     $                     a( jj, 1 ), lda, w( jj, 1 ), ldw, one,
c;     $                     a( jj, jj ), 1 )
c;  100       continue
c;*
c;*           Update the rectangular subdiagonal block
c;*
c;            if( j+jb.le.n )
c;     $         call dgemm( 'No transpose', 'Transpose', n-j-jb+1, jb,
c;     $                     k-1, -one, a( j+jb, 1 ), lda, w( j, 1 ), ldw,
c;     $                     one, a( j+jb, j ), lda )
c;  110    continue
c;*
c;*        Put L21 in standard form by partially undoing the interchanges
c;*        in columns 1:k-1
c;*
c;         j = k - 1
c;  120    continue
c;         jj = j
c;         jp = ipiv( j )
c;         if( jp.lt.0 ) then
c;            jp = -jp
c;            j = j - 1
c;         end if
c;         j = j - 1
c;         if( jp.ne.jj .and. j.ge.1 )
c;     $      call dswap( j, a( jp, 1 ), lda, a( jj, 1 ), lda )
c;         if( j.ge.1 )
c;     $      go to 120
c;*
c;*        Set KB to the number of columns factorized
c;*
c;         kb = k - 1
c;*
c;      end if
c;      return
c;*
c;*     End of DLASYF
c;*
c;      end
cend
cstart .not. unix-dec
      subroutine  dset(n,da,dx,incx)
c
c     sets a vector equal to a constant
c     uses unrolled loops for increment equal to one.
c     modified by mha from linpack dscal, written originally by
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      double precision dx(1), da
      integer i,incx,ix,m,mp1,n

c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      do 10 i = 1,n
        dx(ix) = da
        ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
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
cend
cstart unix .and. (.not. unix-ibm .or. .not. unix-dec)
      integer function idamin(n, sx, incx)
*
*     find smallest index of minimum magnitude of double precision s
*     isamax =  first i,  i = 1 to n,  to minimize  abs(sx(1-incx+i*incx)
*
      double precision sx(1), smin, xmag
      idamin = 0
      if (n .le. 0) return
      idamin = 1
      if (n .le. 1) return
      if (incx .eq. 1) go to 20
*
*        code for increments not equal to 1.
*
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
*
*        code for increments equal to 1.
*
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
*
*     find smallest index of minimum element in double precision s
*     isamax =  first i,  i = 1 to n,  to minimize  sx(1-incx+i*incx)
*
      double precision sx(1), smin, xmag
      idmin = 0
      if (n .le. 0) return
      idmin = 1
      if (n .le. 1) return
      if (incx .eq. 1) go to 20
*
*        code for increments not equal to 1.
*
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
*
*        code for increments equal to 1.
*
   20 smin = abs(sx(1))
      do 30 i = 2, n
         xmag = sx(i)
         if (xmag .ge. smin) go to 30
         idmin = i
         smin = xmag
   30 continue
      return
      end
cend
