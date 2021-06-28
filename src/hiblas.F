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
#if defined(HIB_UNIX_CONVEX)
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
      if(n.le.0)return
      call saxpy(n,da,dx,incx,dy,incy)
      return
      end
      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
      if(n.le.0)return
      call  scopy(n,dx,incx,dy,incy)
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot=0.d0
      if(n.le.0)return
      ddot = sdot(n,dx,incx,dy,incy)
      return
      end
      subroutine drot (n,dx,incx,dy,incy,c,s)
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp,c,s
      integer i,incx,incy,ix,iy,n
      if(n.le.0)return
      call  srot (n,dx,incx,dy,incy,c,s)
      return
      end
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      double precision da,dx(1)
      integer i,incx,ix,m,mp1,n
      if(n.le.0)return
      call  sscal(n,da,dx,incx)
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
      if(n.le.0)return
      call  sswap (n,dx,incx,dy,incy)
      return
      end
#endif
#if defined(HIB_UNIX_NOBLAS)
      double precision function dasum(n,dx,incx)
c
c     returns sum of magnitudes of double precision dx
c     dasum = sum from 0 to n-1 of dabs(dx(1+i*incx))
c
      double precision dx(1)
      dasum = 0.d0
      if(n.le.0)return
      if(incx.eq.1)goto 20
c
c        code for increments not equal to 1
c
      ns = n*incx
          do 10 i=1,ns,incx
          dasum = dasum + dabs(dx(i))
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
         dasum = dasum + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,6
         dasum = dasum + dabs(dx(i)) + dabs(dx(i+1)) + dabs(dx(i+2))
     1   + dabs(dx(i+3)) + dabs(dx(i+4)) + dabs(dx(i+5))
   50 continue
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      subroutine drot (n,dx,incx,dy,incy,c,s)
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
      integer, intent(in), parameter :: n ! the number of displacements to rotate
      real(8), dimension(:), intent(in, out) :: dx, dy  ! the displacement vectors that are about to be rotated 
      integer, intent(in), parameter :: incx, incy ! steps between consecutive elements
      real(8), intent(in), parameter :: c, s  ! cosine and sine of rotation to apply
      double precision dtemp
      integer i,ix,iy
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = c*dx(ix) + s*dy(iy)
        dy(iy) = c*dy(iy) - s*dx(ix)
        dx(ix) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
   20 do 30 i = 1,n
        dtemp = c*dx(i) + s*dy(i)
        dy(i) = c*dy(i) - s*dx(i)
        dx(i) = dtemp
   30 continue
      return
      end
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      double precision da,dx(1)
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
        dx(ix) = da*dx(ix)
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
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      subroutine dsyr  ( uplo, n, alpha, x, incx, a, lda )
*     .. Scalar Arguments ..
      double precision   alpha
      integer            incx, lda, n
      character*1        uplo
*     .. Array Arguments ..
      double precision   a( lda, * ), x( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYR   performs the symmetric rank 1 operation
*
*     A := alpha*x*x' + A,
*
*  where alpha is a real scalar, x is an n element vector and A is an
*  n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced. On exit, the
*           upper triangular part of the array A is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced. On exit, the
*           lower triangular part of the array A is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      double precision   zero
      parameter        ( zero = 0.0d+0 )
*     .. Local Scalars ..
      double precision   temp
      integer            i, info, ix, j, jx, kx
*     .. External Functions ..
      logical            lsame
      external           lsame
*     .. External Subroutines ..
      external           xerbla
*     .. Intrinsic Functions ..
      intrinsic          max
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      info = 0
      if     ( .not.lsame( uplo, 'U' ).and.
     $         .not.lsame( uplo, 'L' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( lda.lt.max( 1, n ) )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'DSYR  ', info )
         return
      end if
*
*     Quick return if possible.
*
      if( ( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
*
*     Set the start point in X if the increment is not unity.
*
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through the triangular part
*     of A.
*
      if( lsame( uplo, 'U' ) )then
*
*        Form  A  when A is stored in upper triangle.
*
         if( incx.eq.1 )then
            do 20, j = 1, n
               if( x( j ).ne.zero )then
                  temp = alpha*x( j )
                  do 10, i = 1, j
                     a( i, j ) = a( i, j ) + x( i )*temp
   10             continue
               end if
   20       continue
         else
            jx = kx
            do 40, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  ix   = kx
                  do 30, i = 1, j
                     a( i, j ) = a( i, j ) + x( ix )*temp
                     ix        = ix        + incx
   30             continue
               end if
               jx = jx + incx
   40       continue
         end if
      else
*
*        Form  A  when A is stored in lower triangle.
*
         if( incx.eq.1 )then
            do 60, j = 1, n
               if( x( j ).ne.zero )then
                  temp = alpha*x( j )
                  do 50, i = j, n
                     a( i, j ) = a( i, j ) + x( i )*temp
   50             continue
               end if
   60       continue
         else
            jx = kx
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  ix   = jx
                  do 70, i = j, n
                     a( i, j ) = a( i, j ) + x( ix )*temp
                     ix        = ix        + incx
   70             continue
               end if
               jx = jx + incx
   80       continue
         end if
      end if
*
      return
*
*     End of DSYR  .
*
      end
      integer function idamax(n,sx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      double precision sx(1),smax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      smax = dabs(sx(ix))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(sx(ix)).le.smax) go to 5
         idamax = i
         smax = dabs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = dabs(sx(1))
      do 30 i = 2,n
         if(dabs(sx(i)).le.smax) go to 30
         idamax = i
         smax = dabs(sx(i))
   30 continue
      return
      end
c$$$cmlib:blas           dnrm2
      double precision function dnrm2 ( n, dx, incx)
      integer          next
      double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/
c
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  dsqrt(u/eps)  over all known machines.
c         cuthi = minimum of  dsqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
   20    go to next,(30, 50, 70, 110)
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 i = j
      assign 110 to next
      sum = (sum / dx(i)) / dx(i)
  105 xmax = dabs(dx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
c
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300
c
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end
#endif
* --------------------------------------------------
* single-precision basic linear algebra (blas) routines for compatability
* with old subroutines
* --------------------------------------------------
#if defined(HIB_NONE)
      integer function isamax(n, sx, incx)
*
*     find smallest index of maximum magnitude of single precision s
*     isamax =  first i,  i = 1 to n,  to minimize  abs(sx(1-incx+i*incx)
*
      double precision sx(1), smax, xmag
      isamax = 0
      if (n .le. 0) return
      isamax = 1
      if (n .le. 1) return
      if (incx .eq. 1) go to 20
*
*        code for increments not equal to 1.
*
      smax = abs(sx(1))
      ns = n * incx
      ii = 1
          do 10 i = 1, ns, incx
          xmag = abs(sx(i))
          if (xmag .le. smax) go to 5
          isamax = ii
          smax = xmag
    5     ii = ii + 1
   10     continue
      return
*
*        code for increments equal to 1.
*
   20 smax = abs(sx(1))
      do 30 i = 2, n
         xmag = abs(sx(i))
         if (xmag .le. smax) go to 30
         isamax = i
         smax = xmag
   30 continue
      return
      end
#endif
#if defined(HIB_NONE)
      subroutine saxpy (n, sa, sx, incx, sy, incy)
*
*     overwrite single precision sy with single precision sa*sx +sy.
*     for i = 0 to n-1,  replace  sy(ly+i*incy) with sa*sx(lx+i*incx) +
*       sy(ly+i*incy),  where lx = 1 if incx .ge. 0,  else lx = (-incx)*n
*       and ly is defined in a similar way using incy.
*
      double precision sx(1), sy(1), sa
      if (n .le. 0 .or. sa .eq. 0.e0) return
      if (incx .eq. incy) if (incx - 1) 5, 20, 60
    5 continue
*
*        code for nonequal or nonpositive increments.
*
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1) * incx + 1
      if (incy .lt. 0) iy = (-n+1) * incy + 1
      do 10 i = 1, n
        sy(iy) = sy(iy) + sa * sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
*
*        code for both increments equal to 1
*
*
*        clean-up loop so remaining vector length is a multiple of 4
*
   20 m = mod(n, 4)
      if ( m  .eq.  0 ) go to 40
      do 30 i = 1, m
        sy(i) = sy(i) + sa * sx(i)
   30 continue
      if ( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1, n, 4
        sy(i) = sy(i) + sa * sx(i)
        sy(i + 1) = sy(i + 1) + sa * sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa * sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa * sx(i + 3)
   50 continue
#endif
c#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS)
#if defined(HIB_NONE)
*        clean-up loop so remaining vector length is a multiple of 6
*
   20 m = mod(n, 6)
      if ( m  .eq.  0 ) go to 40
      do 30 i = 1, m
        sy(i) = sy(i) + sa * sx(i)
   30 continue
      if ( n .lt. 6 ) return
   40 mp1 = m + 1
      do 50 i = mp1, n, 6
        sy(i) = sy(i) + sa * sx(i)
        sy(i + 1) = sy(i + 1) + sa * sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa * sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa * sx(i + 3)
        sy(i + 4) = sy(i + 4) + sa * sx(i + 4)
        sy(i + 5) = sy(i + 5) + sa * sx(i + 5)
   50 continue
#endif
#if defined(HIB_NONE)
      return
*
*        code for equal,  positive,  nonunit increments.
*
   60 continue
      ns = n * incx
          do 70 i = 1, ns, incx
          sy(i) = sa * sx(i) + sy(i)
   70     continue
      return
      end
      subroutine scopy (n, sx, incx, sy, incy)
*
*     copy single precision sx to single precision sy.
*     for i = 0 to n-1,  copy  sx(lx+i*incx) to sy(ly+i*incy),
*     where lx = 1 if incx .ge. 0,  else lx = (-incx)*n,  and ly is
*     defined in a similar way using incy.
*
      double precision sx(1), sy(1)
      if (n .le. 0) return
      if (incx .eq. incy) if (incx - 1) 5, 20, 60
    5 continue
*
*        code for unequal or nonpositive increments.
*
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1) * incx + 1
      if (incy .lt. 0) iy = (-n+1) * incy + 1
      do 10 i = 1, n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
*
*        code for both increments equal to 1
*
*
*        clean-up loop so remaining vector length is a multiple of 7
*
   20 m = mod(n, 7)
      if ( m .eq. 0 ) go to 40
      do 30 i = 1, m
        sy(i) = sx(i)
   30 continue
      if ( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1, n, 7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
   50 continue
      return
*
*        code for equal,  positive,  nonunit increments.
*
   60 continue
      ns = n * incx
          do 70 i = 1, ns, incx
          sy(i) = sx(i)
   70     continue
      return
      end
      double precision function sdot(n, sx, incx, sy, incy)
*     returns the dot product of single precision sx and sy.
*     sdot = sum for i = 0 to n - 1 of  sx(lx + i * incx) * sy(ly + i * incy
*     where lx = 1 if incx .ge. 0,  else lx = (- incx) * n,  and ly is
*     defined in a similar way using incy.
      double precision sx(1), sy(1)
      sdot = 0.0
      if (n .le. 0) return
      if (incx .eq. incy) if (incx - 1) 5, 20, 60
    5 continue
*        code for unequal increments or nonpositive increments.
      ix = 1
      iy = 1
      if (incx .lt. 0)ix = ( - n + 1) * incx + 1
      if (incy .lt. 0)iy = ( - n + 1) * incy + 1
      do 10 i = 1, n
        sdot = sdot + sx(ix) * sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
*         code for both increments equal to 1
#endif
c#if defined(HIB_UNIX_HP) || defined(HIB_MAC)
#if defined(HIB_NONE)
*         clean - up loop so remaining vector length is a multiple of 5
   20 m = mod(n, 5)
      if ( m  .eq.  0 ) go to 40
      do 30 i = 1, m
        sdot = sdot + sx(i) * sy(i)
   30 continue
      if ( n  .lt.  5 ) return
   40 mp1 = m + 1
      do 50 i = mp1, n, 5
        sdot = sdot + sx(i) * sy(i) + sx(i + 1) * sy(i + 1) +
     :                sx(i + 2) * sy(i + 2) + sx(i + 3) * sy(i + 3) +
     :                sx(i + 4) * sy(i + 4)
   50 continue
      return
#endif
c#if defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS)
#if defined(HIB_NONE)
   20 do 30 i = 1, n
        sdot = sdot + sx(i) * sy(i)
*        ncount=ncount+2
   30 continue
      return
#endif
c#if defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_MAC) || defined(HIB_UNIX_IRIS)
#if defined(HIB_NONE)
*         code for positive equal increments .ne. 1
   60 continue
      ns = n * incx
      do 70 i = 1, ns, incx
        sdot = sdot + sx(i) * sy(i)
   70   continue
      return
      end
      subroutine sscal (n, sa, sx, incx)
*
*     replace single precision sx by single precision sa * sx.
*     for i = 0 to n-1,  replace sx(1+i*incx) with  sa * sx(1+i*incx)
*
      double precision sa, sx(1)
      if (n .le. 0) return
      if (incx .eq. 1) go to 20
*
*        code for increments not equal to 1.
*
      ns = n * incx
          do 10 i = 1, ns, incx
          sx(i) = sa * sx(i)
   10     continue
      return
*
*        code for increments equal to 1.
*
*
*        clean-up loop so remaining vector length is a multiple of 5
*
   20 m = mod(n, 5)
      if ( m .eq. 0 ) go to 40
      do 30 i = 1, m
        sx(i) = sa * sx(i)
   30 continue
      if ( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1, n, 5
        sx(i) = sa * sx(i)
        sx(i + 1) = sa * sx(i + 1)
        sx(i + 2) = sa * sx(i + 2)
        sx(i + 3) = sa * sx(i + 3)
        sx(i + 4) = sa * sx(i + 4)
   50 continue
      return
      end
      subroutine sswap (n, sx, incx, sy, incy)
*
*     interchange single precision sx and single precision sy.
*     for i = 0 to n-1,  interchange  sx(lx+i*incx) and sy(ly+i*incy),
*     where lx = 1 if incx .ge. 0,  else lx = (-incx)*n,  and ly is
*     defined in a similar way using incy.
*
      double precision sx(1), sy(1), stemp1, stemp2, stemp3
      if (n .le. 0) return
      if (incx .eq. incy) if (incx - 1) 5, 20, 60
    5 continue
*
*       code for unequal or nonpositive increments.
*
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1) * incx + 1
      if (incy .lt. 0) iy = (-n+1) * incy + 1
      do 10 i = 1, n
        stemp1 = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp1
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
*
*       code for both increments equal to 1
*
*
*       clean-up loop so remaining vector length is a multiple of 3.
*
   20 m = mod(n, 3)
      if ( m  .eq. 0 ) go to 40
      do 30 i = 1, m
        stemp1 = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp1
   30 continue
      if ( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1, n, 3
        stemp1 = sx(i)
        stemp2 = sx(i+1)
        stemp3 = sx(i+2)
        sx(i) = sy(i)
        sx(i+1) = sy(i+1)
        sx(i+2) = sy(i+2)
        sy(i) = stemp1
        sy(i+1) = stemp2
        sy(i+2) = stemp3
   50 continue
      return
   60 continue
*
*     code for equal,  positive,  nonunit increments.
*
      ns = n * incx
        do 70 i = 1, ns, incx
        stemp1 = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp1
   70   continue
      return
      end
#endif
#if defined(HIB_NONE)
      subroutine srot(n,dx,incx,dy,incy,dc,ds)
c
c     multiply the 2 x 2 matrix  ( dc ds) times the 2 x n matrix (dx**t)
c                                (-ds dc)                        (dy**t)
c     where **t indicates transpose.    the elements of dx are in
c     dx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
c     lx = (-incx)*n, and similarly for dy using ly and incy.
      double precision dx,dy,dc,ds,zero,one,w,z
      dimension dx(1),dy(1)
c
      data zero,one/0.d0,1.d0/
      if(n .le. 0 .or. (ds .eq. zero .and. dc .eq. one)) go to 40
      if(.not. (incx .eq. incy .and. incx .gt. 0)) go to 20
c
           nsteps=incx*n
           do 10 i=1,nsteps,incx
                w=dx(i)
                z=dy(i)
                dx(i)=dc*w+ds*z
                dy(i)=-ds*w+dc*z
   10           continue
           go to 40
c
   20 continue
           kx=1
           ky=1
c
           if(incx .lt. 0) kx=1-(n-1)*incx
           if(incy .lt. 0) ky=1-(n-1)*incy
c
           do 30 i=1,n
                w=dx(kx)
                z=dy(ky)
                dx(kx)=dc*w+ds*z
                dy(ky)=-ds*w+dc*z
                kx=kx+incx
                ky=ky+incy
   30           continue
   40 continue
c
      return
      end
#endif
