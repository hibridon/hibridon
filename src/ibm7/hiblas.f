**************************************************************************
*                                                                        *
*                    blas routines library supplement             *
*                                                                        *
**************************************************************************
*                        routines included                               *
*                                                                        *
*  1. blas     basis linear algebra routines from lapack
*     daxpy,dcopy,ddot,drot,dscal,dswap,dsyr,idamax,lsame,ilaenv,xerbla
*     isamax, saxpy, scopy, sdot, sscal, sswap, srot (grot)
*                                                                        *
**************************************************************************
* basic linear algebra (blas) routines
* --------------------------------------------------
* double->single compatability for convex
cstart unix-convex
c;      subroutine daxpy(n,da,dx,incx,dy,incy)
c;c
c;c     constant times a vector plus a vector.
c;c     uses unrolled loops for increments equal to one.
c;c     jack dongarra, linpack, 3/11/78.
c;c
c;      double precision dx(1),dy(1),da
c;      integer i,incx,incy,ix,iy,m,mp1,n
c;      if(n.le.0)return
c;      call saxpy(n,da,dx,incx,dy,incy)
c;      return
c;      end
c;      subroutine  dcopy(n,dx,incx,dy,incy)
c;c
c;c     copies a vector, x, to a vector, y.
c;c     uses unrolled loops for increments equal to one.
c;c     jack dongarra, linpack, 3/11/78.
c;c
c;      double precision dx(1),dy(1)
c;      integer i,incx,incy,ix,iy,m,mp1,n
c;      if(n.le.0)return
c;      call  scopy(n,dx,incx,dy,incy)
c;      return
c;      end
c;      double precision function ddot(n,dx,incx,dy,incy)
c;c
c;c     forms the dot product of two vectors.
c;c     uses unrolled loops for increments equal to one.
c;c     jack dongarra, linpack, 3/11/78.
c;c
c;      double precision dx(1),dy(1),dtemp
c;      integer i,incx,incy,ix,iy,m,mp1,n
c;c
c;      ddot=0.d0
c;      if(n.le.0)return
c;      ddot = sdot(n,dx,incx,dy,incy)
c;      return
c;      end
c;      subroutine  drot (n,dx,incx,dy,incy,c,s)
c;c
c;c     applies a plane rotation.
c;c     jack dongarra, linpack, 3/11/78.
c;c
c;      double precision dx(1),dy(1),dtemp,c,s
c;      integer i,incx,incy,ix,iy,n
c;      if(n.le.0)return
c;      call  srot (n,dx,incx,dy,incy,c,s)
c;      return
c;      end
c;      subroutine  dscal(n,da,dx,incx)
c;c
c;c     scales a vector by a constant.
c;c     uses unrolled loops for increment equal to one.
c;c     jack dongarra, linpack, 3/11/78.
c;c     modified to correct problem with negative increment, 8/21/90.
c;c
c;      double precision da,dx(1)
c;      integer i,incx,ix,m,mp1,n
c;      if(n.le.0)return
c;      call  sscal(n,da,dx,incx)
c;      return
c;      end
c;      subroutine  dswap (n,dx,incx,dy,incy)
c;c
c;c     interchanges two vectors.
c;c     uses unrolled loops for increments equal one.
c;c     jack dongarra, linpack, 3/11/78.
c;c
c;      double precision dx(1),dy(1),dtemp
c;      integer i,incx,incy,ix,iy,m,mp1,n
c;      if(n.le.0)return
c;      call  sswap (n,dx,incx,dy,incy)
c;      return
c;      end
cend
cstart unix-noblas
c;      double precision function dasum(n,dx,incx)
c;c
c;c     returns sum of magnitudes of double precision dx
c;c     dasum = sum from 0 to n-1 of dabs(dx(1+i*incx))
c;c
c;      double precision dx(1)
c;      dasum = 0.d0
c;      if(n.le.0)return
c;      if(incx.eq.1)goto 20
c;c
c;c        code for increments not equal to 1
c;c
c;      ns = n*incx
c;          do 10 i=1,ns,incx
c;          dasum = dasum + dabs(dx(i))
c;   10     continue
c;      return
c;c
c;c        code for increments equal to 1.
c;c
c;c
c;c        clean-up loop so remaining vector length is a multiple of 6.
c;c
c;   20 m = mod(n,6)
c;      if( m .eq. 0 ) go to 40
c;      do 30 i = 1,m
c;         dasum = dasum + dabs(dx(i))
c;   30 continue
c;      if( n .lt. 6 ) return
c;   40 mp1 = m + 1
c;      do 50 i = mp1,n,6
c;         dasum = dasum + dabs(dx(i)) + dabs(dx(i+1)) + dabs(dx(i+2))
c;     1   + dabs(dx(i+3)) + dabs(dx(i+4)) + dabs(dx(i+5))
c;   50 continue
c;      return
c;      end
c;      subroutine daxpy(n,da,dx,incx,dy,incy)
c;c
c;c     constant times a vector plus a vector.
c;c     uses unrolled loops for increments equal to one.
c;c     jack dongarra, linpack, 3/11/78.
c;c
c;      double precision dx(1),dy(1),da
c;      integer i,incx,incy,ix,iy,m,mp1,n
c;c
c;      if(n.le.0)return
c;      if (da .eq. 0.0d0) return
c;      if(incx.eq.1.and.incy.eq.1)go to 20
c;c
c;c        code for unequal increments or equal increments
c;c          not equal to 1
c;c
c;      ix = 1
c;      iy = 1
c;      if(incx.lt.0)ix = (-n+1)*incx + 1
c;      if(incy.lt.0)iy = (-n+1)*incy + 1
c;      do 10 i = 1,n
c;        dy(iy) = dy(iy) + da*dx(ix)
c;        ix = ix + incx
c;        iy = iy + incy
c;   10 continue
c;      return
c;c
c;c        code for both increments equal to 1
c;c
c;c
c;c        clean-up loop
c;c
c;   20 m = mod(n,4)
c;      if( m .eq. 0 ) go to 40
c;      do 30 i = 1,m
c;        dy(i) = dy(i) + da*dx(i)
c;   30 continue
c;      if( n .lt. 4 ) return
c;   40 mp1 = m + 1
c;      do 50 i = mp1,n,4
c;        dy(i) = dy(i) + da*dx(i)
c;        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
c;        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
c;        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
c;   50 continue
c;      return
c;      end
c;      subroutine  dcopy(n,dx,incx,dy,incy)
c;c
c;c     copies a vector, x, to a vector, y.
c;c     uses unrolled loops for increments equal to one.
c;c     jack dongarra, linpack, 3/11/78.
c;c
c;      double precision dx(1),dy(1)
c;      integer i,incx,incy,ix,iy,m,mp1,n
c;c
c;      if(n.le.0)return
c;      if(incx.eq.1.and.incy.eq.1)go to 20
c;c
c;c        code for unequal increments or equal increments
c;c          not equal to 1
c;c
c;      ix = 1
c;      iy = 1
c;      if(incx.lt.0)ix = (-n+1)*incx + 1
c;      if(incy.lt.0)iy = (-n+1)*incy + 1
c;      do 10 i = 1,n
c;        dy(iy) = dx(ix)
c;        ix = ix + incx
c;        iy = iy + incy
c;   10 continue
c;      return
c;c
c;c        code for both increments equal to 1
c;c
c;c
c;c        clean-up loop
c;c
c;   20 m = mod(n,7)
c;      if( m .eq. 0 ) go to 40
c;      do 30 i = 1,m
c;        dy(i) = dx(i)
c;   30 continue
c;      if( n .lt. 7 ) return
c;   40 mp1 = m + 1
c;      do 50 i = mp1,n,7
c;        dy(i) = dx(i)
c;        dy(i + 1) = dx(i + 1)
c;        dy(i + 2) = dx(i + 2)
c;        dy(i + 3) = dx(i + 3)
c;        dy(i + 4) = dx(i + 4)
c;        dy(i + 5) = dx(i + 5)
c;        dy(i + 6) = dx(i + 6)
c;   50 continue
c;      return
c;      end
c;      double precision function ddot(n,dx,incx,dy,incy)
c;c
c;c     forms the dot product of two vectors.
c;c     uses unrolled loops for increments equal to one.
c;c     jack dongarra, linpack, 3/11/78.
c;c
c;      double precision dx(1),dy(1),dtemp
c;      integer i,incx,incy,ix,iy,m,mp1,n
c;c
c;      ddot = 0.0d0
c;      dtemp = 0.0d0
c;      if(n.le.0)return
c;      if(incx.eq.1.and.incy.eq.1)go to 20
c;c
c;c        code for unequal increments or equal increments
c;c          not equal to 1
c;c
c;      ix = 1
c;      iy = 1
c;      if(incx.lt.0)ix = (-n+1)*incx + 1
c;      if(incy.lt.0)iy = (-n+1)*incy + 1
c;      do 10 i = 1,n
c;        dtemp = dtemp + dx(ix)*dy(iy)
c;        ix = ix + incx
c;        iy = iy + incy
c;   10 continue
c;      ddot = dtemp
c;      return
c;c
c;c        code for both increments equal to 1
c;c
c;c
c;c        clean-up loop
c;c
c;   20 m = mod(n,5)
c;      if( m .eq. 0 ) go to 40
c;      do 30 i = 1,m
c;        dtemp = dtemp + dx(i)*dy(i)
c;   30 continue
c;      if( n .lt. 5 ) go to 60
c;   40 mp1 = m + 1
c;      do 50 i = mp1,n,5
c;        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
c;     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
c;   50 continue
c;   60 ddot = dtemp
c;      return
c;      end
c;      subroutine  drot (n,dx,incx,dy,incy,c,s)
c;c
c;c     applies a plane rotation.
c;c     jack dongarra, linpack, 3/11/78.
c;c
c;      double precision dx(1),dy(1),dtemp,c,s
c;      integer i,incx,incy,ix,iy,n
c;c
c;      if(n.le.0)return
c;      if(incx.eq.1.and.incy.eq.1)go to 20
c;c
c;c       code for unequal increments or equal increments not equal
c;c         to 1
c;c
c;      ix = 1
c;      iy = 1
c;      if(incx.lt.0)ix = (-n+1)*incx + 1
c;      if(incy.lt.0)iy = (-n+1)*incy + 1
c;      do 10 i = 1,n
c;        dtemp = c*dx(ix) + s*dy(iy)
c;        dy(iy) = c*dy(iy) - s*dx(ix)
c;        dx(ix) = dtemp
c;        ix = ix + incx
c;        iy = iy + incy
c;   10 continue
c;      return
c;c
c;c       code for both increments equal to 1
c;c
c;   20 do 30 i = 1,n
c;        dtemp = c*dx(i) + s*dy(i)
c;        dy(i) = c*dy(i) - s*dx(i)
c;        dx(i) = dtemp
c;   30 continue
c;      return
c;      end
c;      subroutine  dscal(n,da,dx,incx)
c;c
c;c     scales a vector by a constant.
c;c     uses unrolled loops for increment equal to one.
c;c     jack dongarra, linpack, 3/11/78.
c;c     modified to correct problem with negative increment, 8/21/90.
c;c
c;      double precision da,dx(1)
c;      integer i,incx,ix,m,mp1,n
c;c
c;      if(n.le.0)return
c;      if(incx.eq.1)go to 20
c;c
c;c        code for increment not equal to 1
c;c
c;      ix = 1
c;      if(incx.lt.0)ix = (-n+1)*incx + 1
c;      do 10 i = 1,n
c;        dx(ix) = da*dx(ix)
c;        ix = ix + incx
c;   10 continue
c;      return
c;c
c;c        code for increment equal to 1
c;c
c;c
c;c        clean-up loop
c;c
c;   20 m = mod(n,5)
c;      if( m .eq. 0 ) go to 40
c;      do 30 i = 1,m
c;        dx(i) = da*dx(i)
c;   30 continue
c;      if( n .lt. 5 ) return
c;   40 mp1 = m + 1
c;      do 50 i = mp1,n,5
c;        dx(i) = da*dx(i)
c;        dx(i + 1) = da*dx(i + 1)
c;        dx(i + 2) = da*dx(i + 2)
c;        dx(i + 3) = da*dx(i + 3)
c;        dx(i + 4) = da*dx(i + 4)
c;   50 continue
c;      return
c;      end
c;      subroutine  dswap (n,dx,incx,dy,incy)
c;c
c;c     interchanges two vectors.
c;c     uses unrolled loops for increments equal one.
c;c     jack dongarra, linpack, 3/11/78.
c;c
c;      double precision dx(1),dy(1),dtemp
c;      integer i,incx,incy,ix,iy,m,mp1,n
c;c
c;      if(n.le.0)return
c;      if(incx.eq.1.and.incy.eq.1)go to 20
c;c
c;c       code for unequal increments or equal increments not equal
c;c         to 1
c;c
c;      ix = 1
c;      iy = 1
c;      if(incx.lt.0)ix = (-n+1)*incx + 1
c;      if(incy.lt.0)iy = (-n+1)*incy + 1
c;      do 10 i = 1,n
c;        dtemp = dx(ix)
c;        dx(ix) = dy(iy)
c;        dy(iy) = dtemp
c;        ix = ix + incx
c;        iy = iy + incy
c;   10 continue
c;      return
c;c
c;c       code for both increments equal to 1
c;c
c;c
c;c       clean-up loop
c;c
c;   20 m = mod(n,3)
c;      if( m .eq. 0 ) go to 40
c;      do 30 i = 1,m
c;        dtemp = dx(i)
c;        dx(i) = dy(i)
c;        dy(i) = dtemp
c;   30 continue
c;      if( n .lt. 3 ) return
c;   40 mp1 = m + 1
c;      do 50 i = mp1,n,3
c;        dtemp = dx(i)
c;        dx(i) = dy(i)
c;        dy(i) = dtemp
c;        dtemp = dx(i + 1)
c;        dx(i + 1) = dy(i + 1)
c;        dy(i + 1) = dtemp
c;        dtemp = dx(i + 2)
c;        dx(i + 2) = dy(i + 2)
c;        dy(i + 2) = dtemp
c;   50 continue
c;      return
c;      end
c;      subroutine dsyr  ( uplo, n, alpha, x, incx, a, lda )
c;*     .. Scalar Arguments ..
c;      double precision   alpha
c;      integer            incx, lda, n
c;      character*1        uplo
c;*     .. Array Arguments ..
c;      double precision   a( lda, * ), x( * )
c;*     ..
c;*
c;*  Purpose
c;*  =======
c;*
c;*  DSYR   performs the symmetric rank 1 operation
c;*
c;*     A := alpha*x*x' + A,
c;*
c;*  where alpha is a real scalar, x is an n element vector and A is an
c;*  n by n symmetric matrix.
c;*
c;*  Parameters
c;*  ==========
c;*
c;*  UPLO   - CHARACTER*1.
c;*           On entry, UPLO specifies whether the upper or lower
c;*           triangular part of the array A is to be referenced as
c;*           follows:
c;*
c;*              UPLO = 'U' or 'u'   Only the upper triangular part of A
c;*                                  is to be referenced.
c;*
c;*              UPLO = 'L' or 'l'   Only the lower triangular part of A
c;*                                  is to be referenced.
c;*
c;*           Unchanged on exit.
c;*
c;*  N      - INTEGER.
c;*           On entry, N specifies the order of the matrix A.
c;*           N must be at least zero.
c;*           Unchanged on exit.
c;*
c;*  ALPHA  - DOUBLE PRECISION.
c;*           On entry, ALPHA specifies the scalar alpha.
c;*           Unchanged on exit.
c;*
c;*  X      - DOUBLE PRECISION array of dimension at least
c;*           ( 1 + ( n - 1 )*abs( INCX ) ).
c;*           Before entry, the incremented array X must contain the n
c;*           element vector x.
c;*           Unchanged on exit.
c;*
c;*  INCX   - INTEGER.
c;*           On entry, INCX specifies the increment for the elements of
c;*           X. INCX must not be zero.
c;*           Unchanged on exit.
c;*
c;*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
c;*           Before entry with  UPLO = 'U' or 'u', the leading n by n
c;*           upper triangular part of the array A must contain the upper
c;*           triangular part of the symmetric matrix and the strictly
c;*           lower triangular part of A is not referenced. On exit, the
c;*           upper triangular part of the array A is overwritten by the
c;*           upper triangular part of the updated matrix.
c;*           Before entry with UPLO = 'L' or 'l', the leading n by n
c;*           lower triangular part of the array A must contain the lower
c;*           triangular part of the symmetric matrix and the strictly
c;*           upper triangular part of A is not referenced. On exit, the
c;*           lower triangular part of the array A is overwritten by the
c;*           lower triangular part of the updated matrix.
c;*
c;*  LDA    - INTEGER.
c;*           On entry, LDA specifies the first dimension of A as declared
c;*           in the calling (sub) program. LDA must be at least
c;*           max( 1, n ).
c;*           Unchanged on exit.
c;*
c;*
c;*  Level 2 Blas routine.
c;*
c;*  -- Written on 22-October-1986.
c;*     Jack Dongarra, Argonne National Lab.
c;*     Jeremy Du Croz, Nag Central Office.
c;*     Sven Hammarling, Nag Central Office.
c;*     Richard Hanson, Sandia National Labs.
c;*
c;*
c;*     .. Parameters ..
c;      double precision   zero
c;      parameter        ( zero = 0.0d+0 )
c;*     .. Local Scalars ..
c;      double precision   temp
c;      integer            i, info, ix, j, jx, kx
c;*     .. External Functions ..
c;      logical            lsame
c;      external           lsame
c;*     .. External Subroutines ..
c;      external           xerbla
c;*     .. Intrinsic Functions ..
c;      intrinsic          max
c;*     ..
c;*     .. Executable Statements ..
c;*
c;*     Test the input parameters.
c;*
c;      info = 0
c;      if     ( .not.lsame( uplo, 'U' ).and.
c;     $         .not.lsame( uplo, 'L' )      )then
c;         info = 1
c;      else if( n.lt.0 )then
c;         info = 2
c;      else if( incx.eq.0 )then
c;         info = 5
c;      else if( lda.lt.max( 1, n ) )then
c;         info = 7
c;      end if
c;      if( info.ne.0 )then
c;         call xerbla( 'DSYR  ', info )
c;         return
c;      end if
c;*
c;*     Quick return if possible.
c;*
c;      if( ( n.eq.0 ).or.( alpha.eq.zero ) )
c;     $   return
c;*
c;*     Set the start point in X if the increment is not unity.
c;*
c;      if( incx.le.0 )then
c;         kx = 1 - ( n - 1 )*incx
c;      else if( incx.ne.1 )then
c;         kx = 1
c;      end if
c;*
c;*     Start the operations. In this version the elements of A are
c;*     accessed sequentially with one pass through the triangular part
c;*     of A.
c;*
c;      if( lsame( uplo, 'U' ) )then
c;*
c;*        Form  A  when A is stored in upper triangle.
c;*
c;         if( incx.eq.1 )then
c;            do 20, j = 1, n
c;               if( x( j ).ne.zero )then
c;                  temp = alpha*x( j )
c;                  do 10, i = 1, j
c;                     a( i, j ) = a( i, j ) + x( i )*temp
c;   10             continue
c;               end if
c;   20       continue
c;         else
c;            jx = kx
c;            do 40, j = 1, n
c;               if( x( jx ).ne.zero )then
c;                  temp = alpha*x( jx )
c;                  ix   = kx
c;                  do 30, i = 1, j
c;                     a( i, j ) = a( i, j ) + x( ix )*temp
c;                     ix        = ix        + incx
c;   30             continue
c;               end if
c;               jx = jx + incx
c;   40       continue
c;         end if
c;      else
c;*
c;*        Form  A  when A is stored in lower triangle.
c;*
c;         if( incx.eq.1 )then
c;            do 60, j = 1, n
c;               if( x( j ).ne.zero )then
c;                  temp = alpha*x( j )
c;                  do 50, i = j, n
c;                     a( i, j ) = a( i, j ) + x( i )*temp
c;   50             continue
c;               end if
c;   60       continue
c;         else
c;            jx = kx
c;            do 80, j = 1, n
c;               if( x( jx ).ne.zero )then
c;                  temp = alpha*x( jx )
c;                  ix   = jx
c;                  do 70, i = j, n
c;                     a( i, j ) = a( i, j ) + x( ix )*temp
c;                     ix        = ix        + incx
c;   70             continue
c;               end if
c;               jx = jx + incx
c;   80       continue
c;         end if
c;      end if
c;*
c;      return
c;*
c;*     End of DSYR  .
c;*
c;      end
c;      integer function idamax(n,sx,incx)
c;c
c;c     finds the index of element having max. absolute value.
c;c     jack dongarra, linpack, 3/11/78.
c;c     modified to correct problem with negative increment, 8/21/90.
c;c
c;      double precision sx(1),smax
c;      integer i,incx,ix,n
c;c
c;      idamax = 0
c;      if( n .lt. 1 ) return
c;      idamax = 1
c;      if(n.eq.1)return
c;      if(incx.eq.1)go to 20
c;c
c;c        code for increment not equal to 1
c;c
c;      ix = 1
c;      if(incx.lt.0)ix = (-n+1)*incx + 1
c;      smax = dabs(sx(ix))
c;      ix = ix + incx
c;      do 10 i = 2,n
c;         if(dabs(sx(ix)).le.smax) go to 5
c;         idamax = i
c;         smax = dabs(sx(ix))
c;    5    ix = ix + incx
c;   10 continue
c;      return
c;c
c;c        code for increment equal to 1
c;c
c;   20 smax = dabs(sx(1))
c;      do 30 i = 2,n
c;         if(dabs(sx(i)).le.smax) go to 30
c;         idamax = i
c;         smax = dabs(sx(i))
c;   30 continue
c;      return
c;      end
c;c$$$cmlib:blas           dnrm2
c;      double precision function dnrm2 ( n, dx, incx)
c;      integer          next
c;      double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
c;      data   zero, one /0.0d0, 1.0d0/
c;c
c;c     euclidean norm of the n-vector stored in dx() with storage
c;c     increment incx .
c;c     if    n .le. 0 return with result = 0.
c;c     if n .ge. 1 then incx must be .ge. 1
c;c
c;c           c.l.lawson, 1978 jan 08
c;c
c;c     four phase method     using two built-in constants that are
c;c     hopefully applicable to all machines.
c;c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c;c         cuthi = minimum of  dsqrt(v)      over all known machines.
c;c     where
c;c         eps = smallest no. such that eps + 1. .gt. 1.
c;c         u   = smallest positive no.   (underflow limit)
c;c         v   = largest  no.            (overflow  limit)
c;c
c;c     brief outline of algorithm..
c;c
c;c     phase 1    scans zero components.
c;c     move to phase 2 when a component is nonzero and .le. cutlo
c;c     move to phase 3 when a component is .gt. cutlo
c;c     move to phase 4 when a component is .ge. cuthi/m
c;c     where m = n for x() real and m = 2*n for complex.
c;c
c;c     values for cutlo and cuthi..
c;c     from the environmental parameters listed in the imsl converter
c;c     document the limiting values are as follows..
c;c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c;c                   univac and dec at 2**(-103)
c;c                   thus cutlo = 2**(-51) = 4.44089e-16
c;c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c;c                   thus cuthi = 2**(63.5) = 1.30438e19
c;c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c;c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c;c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c;c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c;c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
c;      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c;c
c;      if(n .gt. 0) go to 10
c;         dnrm2  = zero
c;         go to 300
c;c
c;   10 assign 30 to next
c;      sum = zero
c;      nn = n * incx
c;c                                                 begin main loop
c;      i = 1
c;   20    go to next,(30, 50, 70, 110)
c;   30 if( dabs(dx(i)) .gt. cutlo) go to 85
c;      assign 50 to next
c;      xmax = zero
c;c
c;c                        phase 1.  sum is zero
c;c
c;   50 if( dx(i) .eq. zero) go to 200
c;      if( dabs(dx(i)) .gt. cutlo) go to 85
c;c
c;c                                prepare for phase 2.
c;      assign 70 to next
c;      go to 105
c;c
c;c                                prepare for phase 4.
c;c
c;  100 i = j
c;      assign 110 to next
c;      sum = (sum / dx(i)) / dx(i)
c;  105 xmax = dabs(dx(i))
c;      go to 115
c;c
c;c                   phase 2.  sum is small.
c;c                             scale to avoid destructive underflow.
c;c
c;   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c;c
c;c                     common code for phases 2 and 4.
c;c                     in phase 4 sum is large.  scale to avoid overflow.
c;c
c;  110 if( dabs(dx(i)) .le. xmax ) go to 115
c;         sum = one + sum * (xmax / dx(i))**2
c;         xmax = dabs(dx(i))
c;         go to 200
c;c
c;  115 sum = sum + (dx(i)/xmax)**2
c;      go to 200
c;c
c;c
c;c                  prepare for phase 3.
c;c
c;   75 sum = (sum * xmax) * xmax
c;c
c;c
c;c     for real or d.p. set hitest = cuthi/n
c;c     for complex      set hitest = cuthi/(2*n)
c;c
c;   85 hitest = cuthi/float( n )
c;c
c;c                   phase 3.  sum is mid-range.  no scaling.
c;c
c;      do 95 j =i,nn,incx
c;      if(dabs(dx(j)) .ge. hitest) go to 100
c;   95    sum = sum + dx(j)**2
c;      dnrm2 = dsqrt( sum )
c;      go to 300
c;c
c;  200 continue
c;      i = i + incx
c;      if ( i .le. nn ) go to 20
c;c
c;c              end of main loop.
c;c
c;c              compute square root and adjust for scaling.
c;c
c;      dnrm2 = xmax * dsqrt(sum)
c;  300 continue
c;      return
c;      end
cend
* --------------------------------------------------
* single-precision basic linear algebra (blas) routines for compatability
* with old subroutines
* --------------------------------------------------
cstart none
c;      integer function isamax(n, sx, incx)
c;*
c;*     find smallest index of maximum magnitude of single precision s
c;*     isamax =  first i,  i = 1 to n,  to minimize  abs(sx(1-incx+i*incx)
c;*
c;      double precision sx(1), smax, xmag
c;      isamax = 0
c;      if (n .le. 0) return
c;      isamax = 1
c;      if (n .le. 1) return
c;      if (incx .eq. 1) go to 20
c;*
c;*        code for increments not equal to 1.
c;*
c;      smax = abs(sx(1))
c;      ns = n * incx
c;      ii = 1
c;          do 10 i = 1, ns, incx
c;          xmag = abs(sx(i))
c;          if (xmag .le. smax) go to 5
c;          isamax = ii
c;          smax = xmag
c;    5     ii = ii + 1
c;   10     continue
c;      return
c;*
c;*        code for increments equal to 1.
c;*
c;   20 smax = abs(sx(1))
c;      do 30 i = 2, n
c;         xmag = abs(sx(i))
c;         if (xmag .le. smax) go to 30
c;         isamax = i
c;         smax = xmag
c;   30 continue
c;      return
c;      end
cend
cstart none
c;      subroutine saxpy (n, sa, sx, incx, sy, incy)
c;*
c;*     overwrite single precision sy with single precision sa*sx +sy.
c;*     for i = 0 to n-1,  replace  sy(ly+i*incy) with sa*sx(lx+i*incx) +
c;*       sy(ly+i*incy),  where lx = 1 if incx .ge. 0,  else lx = (-incx)*n
c;*       and ly is defined in a similar way using incy.
c;*
c;      double precision sx(1), sy(1), sa
c;      if (n .le. 0 .or. sa .eq. 0.e0) return
c;      if (incx .eq. incy) if (incx - 1) 5, 20, 60
c;    5 continue
c;*
c;*        code for nonequal or nonpositive increments.
c;*
c;      ix = 1
c;      iy = 1
c;      if (incx .lt. 0) ix = (-n+1) * incx + 1
c;      if (incy .lt. 0) iy = (-n+1) * incy + 1
c;      do 10 i = 1, n
c;        sy(iy) = sy(iy) + sa * sx(ix)
c;        ix = ix + incx
c;        iy = iy + incy
c;   10 continue
c;      return
c;*
c;*        code for both increments equal to 1
c;*
c;*
c;*        clean-up loop so remaining vector length is a multiple of 4
c;*
c;   20 m = mod(n, 4)
c;      if ( m  .eq.  0 ) go to 40
c;      do 30 i = 1, m
c;        sy(i) = sy(i) + sa * sx(i)
c;   30 continue
c;      if ( n .lt. 4 ) return
c;   40 mp1 = m + 1
c;      do 50 i = mp1, n, 4
c;        sy(i) = sy(i) + sa * sx(i)
c;        sy(i + 1) = sy(i + 1) + sa * sx(i + 1)
c;        sy(i + 2) = sy(i + 2) + sa * sx(i + 2)
c;        sy(i + 3) = sy(i + 3) + sa * sx(i + 3)
c;   50 continue
cend
ccstart unix-dec unix-iris
cstart none
c;*        clean-up loop so remaining vector length is a multiple of 6
c;*
c;   20 m = mod(n, 6)
c;      if ( m  .eq.  0 ) go to 40
c;      do 30 i = 1, m
c;        sy(i) = sy(i) + sa * sx(i)
c;   30 continue
c;      if ( n .lt. 6 ) return
c;   40 mp1 = m + 1
c;      do 50 i = mp1, n, 6
c;        sy(i) = sy(i) + sa * sx(i)
c;        sy(i + 1) = sy(i + 1) + sa * sx(i + 1)
c;        sy(i + 2) = sy(i + 2) + sa * sx(i + 2)
c;        sy(i + 3) = sy(i + 3) + sa * sx(i + 3)
c;        sy(i + 4) = sy(i + 4) + sa * sx(i + 4)
c;        sy(i + 5) = sy(i + 5) + sa * sx(i + 5)
c;   50 continue
cend
cstart none
c;      return
c;*
c;*        code for equal,  positive,  nonunit increments.
c;*
c;   60 continue
c;      ns = n * incx
c;          do 70 i = 1, ns, incx
c;          sy(i) = sa * sx(i) + sy(i)
c;   70     continue
c;      return
c;      end
c;      subroutine scopy (n, sx, incx, sy, incy)
c;*
c;*     copy single precision sx to single precision sy.
c;*     for i = 0 to n-1,  copy  sx(lx+i*incx) to sy(ly+i*incy),
c;*     where lx = 1 if incx .ge. 0,  else lx = (-incx)*n,  and ly is
c;*     defined in a similar way using incy.
c;*
c;      double precision sx(1), sy(1)
c;      if (n .le. 0) return
c;      if (incx .eq. incy) if (incx - 1) 5, 20, 60
c;    5 continue
c;*
c;*        code for unequal or nonpositive increments.
c;*
c;      ix = 1
c;      iy = 1
c;      if (incx .lt. 0) ix = (-n+1) * incx + 1
c;      if (incy .lt. 0) iy = (-n+1) * incy + 1
c;      do 10 i = 1, n
c;        sy(iy) = sx(ix)
c;        ix = ix + incx
c;        iy = iy + incy
c;   10 continue
c;      return
c;*
c;*        code for both increments equal to 1
c;*
c;*
c;*        clean-up loop so remaining vector length is a multiple of 7
c;*
c;   20 m = mod(n, 7)
c;      if ( m .eq. 0 ) go to 40
c;      do 30 i = 1, m
c;        sy(i) = sx(i)
c;   30 continue
c;      if ( n .lt. 7 ) return
c;   40 mp1 = m + 1
c;      do 50 i = mp1, n, 7
c;        sy(i) = sx(i)
c;        sy(i + 1) = sx(i + 1)
c;        sy(i + 2) = sx(i + 2)
c;        sy(i + 3) = sx(i + 3)
c;        sy(i + 4) = sx(i + 4)
c;        sy(i + 5) = sx(i + 5)
c;        sy(i + 6) = sx(i + 6)
c;   50 continue
c;      return
c;*
c;*        code for equal,  positive,  nonunit increments.
c;*
c;   60 continue
c;      ns = n * incx
c;          do 70 i = 1, ns, incx
c;          sy(i) = sx(i)
c;   70     continue
c;      return
c;      end
c;      double precision function sdot(n, sx, incx, sy, incy)
c;*     returns the dot product of single precision sx and sy.
c;*     sdot = sum for i = 0 to n - 1 of  sx(lx + i * incx) * sy(ly + i * incy
c;*     where lx = 1 if incx .ge. 0,  else lx = (- incx) * n,  and ly is
c;*     defined in a similar way using incy.
c;      double precision sx(1), sy(1)
c;      sdot = 0.0
c;      if (n .le. 0) return
c;      if (incx .eq. incy) if (incx - 1) 5, 20, 60
c;    5 continue
c;*        code for unequal increments or nonpositive increments.
c;      ix = 1
c;      iy = 1
c;      if (incx .lt. 0)ix = ( - n + 1) * incx + 1
c;      if (incy .lt. 0)iy = ( - n + 1) * incy + 1
c;      do 10 i = 1, n
c;        sdot = sdot + sx(ix) * sy(iy)
c;        ix = ix + incx
c;        iy = iy + incy
c;   10 continue
c;      return
c;*         code for both increments equal to 1
cend
ccstart unix-hp mac
cstart none
c;*         clean - up loop so remaining vector length is a multiple of 5
c;   20 m = mod(n, 5)
c;      if ( m  .eq.  0 ) go to 40
c;      do 30 i = 1, m
c;        sdot = sdot + sx(i) * sy(i)
c;   30 continue
c;      if ( n  .lt.  5 ) return
c;   40 mp1 = m + 1
c;      do 50 i = mp1, n, 5
c;        sdot = sdot + sx(i) * sy(i) + sx(i + 1) * sy(i + 1) +
c;     :                sx(i + 2) * sy(i + 2) + sx(i + 3) * sy(i + 3) +
c;     :                sx(i + 4) * sy(i + 4)
c;   50 continue
c;      return
cend
ccstart unix-dec unix-iris
cstart none
c;   20 do 30 i = 1, n
c;        sdot = sdot + sx(i) * sy(i)
c;*        ncount=ncount+2
c;   30 continue
c;      return
cend
ccstart unix-hp unix-dec mac unix-iris
cstart none
c;*         code for positive equal increments .ne. 1
c;   60 continue
c;      ns = n * incx
c;      do 70 i = 1, ns, incx
c;        sdot = sdot + sx(i) * sy(i)
c;   70   continue
c;      return
c;      end
c;      subroutine sscal (n, sa, sx, incx)
c;*
c;*     replace single precision sx by single precision sa * sx.
c;*     for i = 0 to n-1,  replace sx(1+i*incx) with  sa * sx(1+i*incx)
c;*
c;      double precision sa, sx(1)
c;      if (n .le. 0) return
c;      if (incx .eq. 1) go to 20
c;*
c;*        code for increments not equal to 1.
c;*
c;      ns = n * incx
c;          do 10 i = 1, ns, incx
c;          sx(i) = sa * sx(i)
c;   10     continue
c;      return
c;*
c;*        code for increments equal to 1.
c;*
c;*
c;*        clean-up loop so remaining vector length is a multiple of 5
c;*
c;   20 m = mod(n, 5)
c;      if ( m .eq. 0 ) go to 40
c;      do 30 i = 1, m
c;        sx(i) = sa * sx(i)
c;   30 continue
c;      if ( n .lt. 5 ) return
c;   40 mp1 = m + 1
c;      do 50 i = mp1, n, 5
c;        sx(i) = sa * sx(i)
c;        sx(i + 1) = sa * sx(i + 1)
c;        sx(i + 2) = sa * sx(i + 2)
c;        sx(i + 3) = sa * sx(i + 3)
c;        sx(i + 4) = sa * sx(i + 4)
c;   50 continue
c;      return
c;      end
c;      subroutine sswap (n, sx, incx, sy, incy)
c;*
c;*     interchange single precision sx and single precision sy.
c;*     for i = 0 to n-1,  interchange  sx(lx+i*incx) and sy(ly+i*incy),
c;*     where lx = 1 if incx .ge. 0,  else lx = (-incx)*n,  and ly is
c;*     defined in a similar way using incy.
c;*
c;      double precision sx(1), sy(1), stemp1, stemp2, stemp3
c;      if (n .le. 0) return
c;      if (incx .eq. incy) if (incx - 1) 5, 20, 60
c;    5 continue
c;*
c;*       code for unequal or nonpositive increments.
c;*
c;      ix = 1
c;      iy = 1
c;      if (incx .lt. 0) ix = (-n+1) * incx + 1
c;      if (incy .lt. 0) iy = (-n+1) * incy + 1
c;      do 10 i = 1, n
c;        stemp1 = sx(ix)
c;        sx(ix) = sy(iy)
c;        sy(iy) = stemp1
c;        ix = ix + incx
c;        iy = iy + incy
c;   10 continue
c;      return
c;*
c;*       code for both increments equal to 1
c;*
c;*
c;*       clean-up loop so remaining vector length is a multiple of 3.
c;*
c;   20 m = mod(n, 3)
c;      if ( m  .eq. 0 ) go to 40
c;      do 30 i = 1, m
c;        stemp1 = sx(i)
c;        sx(i) = sy(i)
c;        sy(i) = stemp1
c;   30 continue
c;      if ( n .lt. 3 ) return
c;   40 mp1 = m + 1
c;      do 50 i = mp1, n, 3
c;        stemp1 = sx(i)
c;        stemp2 = sx(i+1)
c;        stemp3 = sx(i+2)
c;        sx(i) = sy(i)
c;        sx(i+1) = sy(i+1)
c;        sx(i+2) = sy(i+2)
c;        sy(i) = stemp1
c;        sy(i+1) = stemp2
c;        sy(i+2) = stemp3
c;   50 continue
c;      return
c;   60 continue
c;*
c;*     code for equal,  positive,  nonunit increments.
c;*
c;      ns = n * incx
c;        do 70 i = 1, ns, incx
c;        stemp1 = sx(i)
c;        sx(i) = sy(i)
c;        sy(i) = stemp1
c;   70   continue
c;      return
c;      end
cend
cstart none
c;      subroutine srot(n,dx,incx,dy,incy,dc,ds)
c;c
c;c     multiply the 2 x 2 matrix  ( dc ds) times the 2 x n matrix (dx**t)
c;c                                (-ds dc)                        (dy**t)
c;c     where **t indicates transpose.    the elements of dx are in
c;c     dx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
c;c     lx = (-incx)*n, and similarly for dy using ly and incy.
c;      double precision dx,dy,dc,ds,zero,one,w,z
c;      dimension dx(1),dy(1)
c;c
c;      data zero,one/0.d0,1.d0/
c;      if(n .le. 0 .or. (ds .eq. zero .and. dc .eq. one)) go to 40
c;      if(.not. (incx .eq. incy .and. incx .gt. 0)) go to 20
c;c
c;           nsteps=incx*n
c;           do 10 i=1,nsteps,incx
c;                w=dx(i)
c;                z=dy(i)
c;                dx(i)=dc*w+ds*z
c;                dy(i)=-ds*w+dc*z
c;   10           continue
c;           go to 40
c;c
c;   20 continue
c;           kx=1
c;           ky=1
c;c
c;           if(incx .lt. 0) kx=1-(n-1)*incx
c;           if(incy .lt. 0) ky=1-(n-1)*incy
c;c
c;           do 30 i=1,n
c;                w=dx(kx)
c;                z=dy(ky)
c;                dx(kx)=dc*w+ds*z
c;                dy(ky)=-ds*w+dc*z
c;                kx=kx+incx
c;                ky=ky+incy
c;   30           continue
c;   40 continue
c;c
c;      return
c;      end
cend
