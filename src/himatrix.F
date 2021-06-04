* NB #if defined(HIB_UNIX_AIX) rather than unix-ibm for fortran not essl routines
#if defined(HIB_UNIX_IBM)
c;      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
c;c
c;      integer n,nm,ierr,matz
c;      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c;c
c;c     this subroutine calls the essl routine dspev to get the
c;c     eigenvalues and eigenvectors of a real symmetric matrix.
c;c
c;c     on input
c;c
c;c        nm  must be set to the row dimension of the two-dimensional
c;c        array parameters as declared in the calling program
c;c        dimension statement.
c;c
c;c        n  is the order of the matrix  a.
c;c
c;c        a  contains the real symmetric matrix.
c;c
c;c        matz  is an integer variable set equal to zero if
c;c        only eigenvalues are desired.  otherwise it is set to
c;c        any non-zero integer for both eigenvalues and eigenvectors.
c;c
c;c     on output
c;c
c;c        w  contains the eigenvalues in ascending order.
c;c
c;c        z  contains the eigenvectors if matz is not zero.


c;c        ierr  is an integer output variable set equal to an error
c;c           completion code described in the documentation for tqlrat
c;c           and tql2.  the normal completion code is zero.
c;c
c;c        fv1  and  fv2  are temporary storage arrays.
c;c
c;c
c;c     this version dated august 1983.
c;c
c;c     ------------------------------------------------------------------
c;c
c;      ione = 1
c;      naux=2*nm
c;* compress matrix into lower triangle
c;      call triang(a,nm,n)
c;* call essl routine
c;      call dspev(ione,a,w,z,nm,n,fv1,naux)
c;      return
c;      end
c;      subroutine triang(a,nrow,n)
c;* -------------------
c;* subroutine to pack lower triangle of symmetric matrix, stored as
c;* column form
c;* input
c;*   a:         on input: full matrix stored column by column
c;*              on return: packed lower triangle
c;*   nrow:      maximum row dimension of matrix
c;*   n:         order of matrix
c;* written by:  millard alexander
c;* current revision date;  1-aug-91
c;* -------------------
c;      implicit double precision (a-h,o-z)
c;      dimension a(25)
c;      if (n .gt. 1) then
c;        indnew=n+1
c;        indold=nrow+2
c;        do  50 i = 2, n
c;          nelt=n-i+1
c;          call dcopy(nelt,a(indold),1,a(indnew),1)
c;          indold=indold+nrow+1
c;          indnew=indnew+nelt
c;50      continue
c;      endif
c;      return
c;      end
#endif
c#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
#if defined(HIB_UNIX) && !defined(HIB_UNIX_IBM) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;      subroutine rsg(nm,n,a,b,w,matz,z,fv1,fv2,ierr)
c;c
c;      integer n,nm,ierr,matz
c;      double precision a(nm,n),b(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c;c
c;c     this subroutine calls the recommended sequence of
c;c     subroutines from the eigensystem subroutine package (eispack)
c;c     to find the eigenvalues and eigenvectors (if desired)
c;c     for the real symmetric generalized eigenproblem  ax = (lambda)bx.
c;c
c;c     on input
c;c
c;c        nm  must be set to the row dimension of the two-dimensional
c;c        array parameters as declared in the calling program
c;c        dimension statement.
c;c
c;c        n  is the order of the matrices  a  and  b.
c;c
c;c        a  contains a real symmetric matrix.
c;c
c;c        b  contains a positive definite real symmetric matrix.
c;c
c;c        matz  is an integer variable set equal to zero if
c;c        only eigenvalues are desired.  otherwise it is set to
c;c        any non-zero integer for both eigenvalues and eigenvectors.
c;c
c;c     on output
c;c
c;c        w  contains the eigenvalues in ascending order.
c;c
c;c        z  contains the eigenvectors if matz is not zero.
c;c
c;c        ierr  is an integer output variable set equal to an error
c;c           completion code described in the documentation for tqlrat
c;c           and tql2.  the normal completion code is zero.
c;c
c;c        fv1  and  fv2  are temporary storage arrays.
c;c
c;c     questions and comments should be directed to burton s. garbow,
c;c     mathematics and computer science div, argonne national laboratory
c;c
c;c     this version dated august 1983.
c;c
c;c     ------------------------------------------------------------------
c;c
c;      if (n .le. nm) go to 10
c;      ierr = 10 * n
c;      go to 50
c;c
c;   10 call  reduc(nm,n,a,b,fv2,ierr)
c;      if (ierr .ne. 0) go to 50
c;      if (matz .ne. 0) go to 20
c;c     .......... find eigenvalues only ..........
c;      call  tred1(nm,n,a,w,fv1,fv2)
c;      call  tqlrat(n,w,fv2,ierr)
c;      go to 50
c;c     .......... find both eigenvalues and eigenvectors ..........
c;   20 call  tred2(nm,n,a,w,fv1,z)
c;      call  tql2(nm,n,w,fv1,z,ierr)
c;      if (ierr .ne. 0) go to 50
c;      call  rebak(nm,n,b,fv2,n,z)
c;   50 return
c;      end
c;      subroutine rebak(nm,n,b,dl,m,z)
c;c
c;      integer i,j,k,m,n,i1,ii,nm
c;      double precision b(nm,n),dl(n),z(nm,m)
c;      double precision x
c;c
c;c     this subroutine is a translation of the algol procedure rebaka,
c;c     num. math. 11, 99-110(1968) by martin and wilkinson.
c;c     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
c;c
c;c     this subroutine forms the eigenvectors of a generalized
c;c     symmetric eigensystem by back transforming those of the
c;c     derived symmetric matrix determined by  reduc.
c;c
c;c     on input
c;c
c;c        nm must be set to the row dimension of two-dimensional
c;c          array parameters as declared in the calling program
c;c          dimension statement.
c;c
c;c        n is the order of the matrix system.
c;c
c;c        b contains information about the similarity transformation
c;c          (cholesky decomposition) used in the reduction by  reduc
c;c          in its strict lower triangle.
c;c
c;c        dl contains further information about the transformation.
c;c
c;c        m is the number of eigenvectors to be back transformed.
c;c
c;c        z contains the eigenvectors to be back transformed
c;c          in its first m columns.
c;c
c;c     on output
c;c
c;c        z contains the transformed eigenvectors
c;c          in its first m columns.
c;c
c;c     questions and comments should be directed to burton s. garbow,
c;c     mathematics and computer science div, argonne national laboratory
c;c
c;c     this version dated august 1983.
c;c
c;c     ------------------------------------------------------------------
c;c
c;      if (m .eq. 0) go to 200
c;c
c;      do 100 j = 1, m
c;c     .......... for i=n step -1 until 1 do -- ..........
c;         do 100 ii = 1, n
c;            i = n + 1 - ii
c;            i1 = i + 1
c;            x = z(i,j)
c;            if (i .eq. n) go to 80
c;c
c;            do 60 k = i1, n
c;   60       x = x - b(k,i) * z(k,j)
c;c
c;   80       z(i,j) = x / dl(i)
c;  100 continue
c;c
c;  200 return
c;      end
c;      subroutine reduc(nm,n,a,b,dl,ierr)
c;c
c;      integer i,j,k,n,i1,j1,nm,nn,ierr
c;      double precision a(nm,n),b(nm,n),dl(n)
c;      double precision x,y
c;c
c;c     this subroutine is a translation of the algol procedure reduc1,
c;c     num. math. 11, 99-110(1968) by martin and wilkinson.
c;c     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
c;c
c;c     this subroutine reduces the generalized symmetric eigenproblem
c;c     ax=(lambda)bx, where b is positive definite, to the standard
c;c     symmetric eigenproblem using the cholesky factorization of b.
c;c
c;c     on input
c;c
c;c        nm must be set to the row dimension of two-dimensional
c;c          array parameters as declared in the calling program
c;c          dimension statement.
c;c
c;c        n is the order of the matrices a and b.  if the cholesky
c;c          factor l of b is already available, n should be prefixed
c;c          with a minus sign.
c;c
c;c        a and b contain the real symmetric input matrices.  only the
c;c          full upper triangles of the matrices need be supplied.  if
c;c          n is negative, the strict lower triangle of b contains,
c;c          instead, the strict lower triangle of its cholesky factor l.
c;c
c;c        dl contains, if n is negative, the diagonal elements of l.
c;c
c;c     on output
c;c
c;c        a contains in its full lower triangle the full lower triangle
c;c          of the symmetric matrix derived from the reduction to the
c;c          standard form.  the strict upper triangle of a is unaltered.
c;c
c;c        b contains in its strict lower triangle the strict lower
c;c          triangle of its cholesky factor l.  the full upper
c;c          triangle of b is unaltered.
c;c
c;c        dl contains the diagonal elements of l.
c;c
c;c        ierr is set to
c;c          zero       for normal return,
c;c          7*n+1      if b is not positive definite.
c;c
c;c     questions and comments should be directed to burton s. garbow,
c;c     mathematics and computer science div, argonne national laboratory
c;c
c;c     this version dated august 1983.
c;c
c;c     ------------------------------------------------------------------
c;c
c;      ierr = 0
c;      nn = iabs(n)
c;      if (n .lt. 0) go to 100
c;c     .......... form l in the arrays b and dl ..........
c;      do 80 i = 1, n
c;         i1 = i - 1
c;c
c;         do 80 j = i, n
c;            x = b(i,j)
c;            if (i .eq. 1) go to 40
c;c
c;            do 20 k = 1, i1
c;   20       x = x - b(i,k) * b(j,k)
c;c
c;   40       if (j .ne. i) go to 60
c;            if (x .le. 0.0d0) go to 1000
c;            y = dsqrt(x)
c;            dl(i) = y
c;            go to 80
c;   60       b(j,i) = x / y
c;   80 continue
c;c     .......... form the transpose of the upper triangle of inv(l)*a
c;c                in the lower triangle of the array a ..........
c;  100 do 200 i = 1, nn
c;         i1 = i - 1
c;         y = dl(i)
c;c
c;         do 200 j = i, nn
c;            x = a(i,j)
c;            if (i .eq. 1) go to 180
c;c
c;            do 160 k = 1, i1
c;  160       x = x - b(i,k) * a(j,k)
c;c
c;  180       a(j,i) = x / y
c;  200 continue
c;c     .......... pre-multiply by inv(l) and overwrite ..........
c;      do 300 j = 1, nn
c;         j1 = j - 1
c;c
c;         do 300 i = j, nn
c;            x = a(i,j)
c;            if (i .eq. j) go to 240
c;            i1 = i - 1
c;c
c;            do 220 k = j, i1
c;  220       x = x - a(k,j) * b(i,k)
c;c
c;  240       if (j .eq. 1) go to 280
c;c
c;            do 260 k = 1, j1
c;  260       x = x - a(j,k) * b(i,k)
c;c
c;  280       a(i,j) = x / dl(i)
c;  300 continue
c;c
c;      go to 1001
c;c     .......... set error -- b is not positive definite ..........
c; 1000 ierr = 7 * n + 1
c; 1001 return
c;      end
c;*** from netlib, tue may 24 12:02:59 cdt 1988 ***
c;*  eispack eigenvalue package
c;      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
c;c
c;      integer n,nm,ierr,matz
c;      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
c;c
c;c     this subroutine calls the recommended sequence of
c;c     subroutines from the eigensystem subroutine package (eispack)
c;c     to find the eigenvalues and eigenvectors (if desired)
c;c     of a real symmetric matrix.
c;c
c;c     on input
c;c
c;c        nm  must be set to the row dimension of the two-dimensional
c;c        array parameters as declared in the calling program
c;c        dimension statement.
c;c
c;c        n  is the order of the matrix  a.
c;c
c;c        a  contains the real symmetric matrix.
c;c
c;c        matz  is an integer variable set equal to zero if
c;c        only eigenvalues are desired.  otherwise it is set to
c;c        any non-zero integer for both eigenvalues and eigenvectors.
c;c
c;c     on output
c;c
c;c        w  contains the eigenvalues in ascending order.
c;c
c;c        z  contains the eigenvectors if matz is not zero.
c;c
c;c        ierr  is an integer output variable set equal to an error
c;c           completion code described in the documentation for tqlrat
c;c           and tql2.  the normal completion code is zero.
c;c
c;c        fv1  and  fv2  are temporary storage arrays.
c;c
c;c     questions and comments should be directed to burton s. garbow,
c;c     mathematics and computer science div, argonne national laboratory
c;c
c;c     this version dated august 1983.
c;c
c;c     ------------------------------------------------------------------
c;c
c;      if (n .le. nm) go to 10
c;      ierr = 10 * n
c;      go to 50
c;c
c;   10 if (matz .ne. 0) go to 20
c;c     .......... find eigenvalues only ..........
c;      call  tred1(nm,n,a,w,fv1,fv2)
c;*  tqlrat encounters catastrophic underflow on the vax
c;      call  tqlrat(n,w,fv2,ierr)
c;*     call  tql1(n,w,fv1,ierr)
c;      go to 50
c;c     .......... find both eigenvalues and eigenvectors ..........
c;   20 call  tred2(nm,n,a,w,fv1,z)
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_IBM) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;      call  tql2(nm,n,w,fv1,z,ierr)
c;   50 return
c;      end
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
c;      subroutine tql2(nm,n,d,e,z,ierr)
c;c
c;      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
c;      double precision d(n),e(n),z(nm,n)
c;      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
c;c
c;c     this subroutine is a translation of the algol procedure tql2,
c;c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c;c     wilkinson.
c;c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c;c
c;c     this subroutine finds the eigenvalues and eigenvectors
c;c     of a symmetric tridiagonal matrix by the ql method.
c;c     the eigenvectors of a full symmetric matrix can also
c;c     be found if  tred2  has been used to reduce this
c;c     full matrix to tridiagonal form.
c;c
c;c     on input
c;c
c;c        nm must be set to the row dimension of two-dimensional
c;c          array parameters as declared in the calling program
c;c          dimension statement.
c;c
c;c        n is the order of the matrix.
c;c
c;c        d contains the diagonal elements of the input matrix.
c;c
c;c        e contains the subdiagonal elements of the input matrix
c;c          in its last n-1 positions.  e(1) is arbitrary.
c;c
c;c        z contains the transformation matrix produced in the
c;c          reduction by  tred2, if performed.  if the eigenvectors
c;c          of the tridiagonal matrix are desired, z must contain
c;c          the identity matrix.
c;c
c;c      on output
c;c
c;c        d contains the eigenvalues in ascending order.  if an
c;c          error exit is made, the eigenvalues are correct but
c;c          unordered for indices 1,2,...,ierr-1.
c;c
c;c        e has been destroyed.
c;c
c;c        z contains orthonormal eigenvectors of the symmetric
c;c          tridiagonal (or full) matrix.  if an error exit is made,
c;c          z contains the eigenvectors associated with the stored
c;c          eigenvalues.
c;c
c;c        ierr is set to
c;c          zero       for normal return,
c;c          j          if the j-th eigenvalue has not been
c;c                     determined after 30 iterations.
c;c
c;c     calls pythag for  dsqrt(a*a + b*b) .
c;c
c;c     questions and comments should be directed to burton s. garbow,
c;c     mathematics and computer science div, argonne national laboratory
c;c
c;c     this version dated august 1983.
c;c
c;c     ------------------------------------------------------------------
c;c
c;      ierr = 0
c;      if (n .eq. 1) go to 1001
c;c
c;      do 100 i = 2, n
c;  100 e(i-1) = e(i)
c;c
c;      f = 0.0d0
c;      tst1 = 0.0d0
c;      e(n) = 0.0d0
c;c
c;      do 240 l = 1, n
c;         j = 0
c;         h = dabs(d(l)) + dabs(e(l))
c;         if (tst1 .lt. h) tst1 = h
c;c     .......... look for small sub-diagonal element ..........
c;         do 110 m = l, n
c;            tst2 = tst1 + dabs(e(m))
c;            if (tst2 .eq. tst1) go to 120
c;c     .......... e(n) is always zero, so there is no exit
c;c                through the bottom of the loop ..........
c;  110    continue
c;c
c;  120    if (m .eq. l) go to 220
c;  130    if (j .eq. 30) go to 1000
c;         j = j + 1
c;c     .......... form shift ..........
c;         l1 = l + 1
c;         l2 = l1 + 1
c;         g = d(l)
c;         p = (d(l1) - g) / (2.0d0 * e(l))
c;         r = pythag(p,1.0d0)
c;         d(l) = e(l) / (p + dsign(r,p))
c;         d(l1) = e(l) * (p + dsign(r,p))
c;         dl1 = d(l1)
c;         h = g - d(l)
c;         if (l2 .gt. n) go to 145
c;c
c;         do 140 i = l2, n
c;  140    d(i) = d(i) - h
c;c
c;  145    f = f + h
c;c     .......... ql transformation ..........
c;         p = d(m)
c;         c = 1.0d0
c;         c2 = c
c;         el1 = e(l1)
c;         s = 0.0d0
c;         mml = m - l
c;c     .......... for i=m-1 step -1 until l do -- ..........
c;c        do 200 ii = 1, mml
c;         do 200 i=m-1,l,-1
c;            c3 = c2
c;            c2 = c
c;            s2 = s
c;c           i = m - ii
c;            g = c * e(i)
c;            h = c * p
c;            r = pythag(p,e(i))
c;            e(i+1) = s * r
c;            s = e(i) / r
c;            c = p / r
c;            p = c * d(i) - s * g
c;            d(i+1) = h + s * (c * g + s * d(i))
c;c  replace this inner loop with blas call (millard alexander 12/24/90)
c;c     .......... form vector ..........
c;c            do 180 k = 1, n
c;c               h = z(k,i+1)
c;c               hh = z(k, i)
c;c               z(k,i+1) = s * hh + c * h
c;c               z(k,i) = c * hh - s * h
c;c  180       continue
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
c;           call drot (n, z(1,i+1), 1, z(1,i), 1, c, s)
#endif
#if defined(HIB_UNIX_CONVEX)
c;           call srot (n, z(1,i+1), 1, z(1,i), 1, c, s)
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
c;  200    continue
c;c
c;         p = -s * s2 * c3 * el1 * e(l) / dl1
c;         e(l) = s * p
c;         d(l) = c * p
c;         tst2 = tst1 + dabs(e(l))
c;         if (tst2 .gt. tst1) go to 130
c;  220    d(l) = d(l) + f
c;  240 continue
c;c     .......... order eigenvalues and eigenvectors ..........
c;      do 300 ii = 2, n
c;         i = ii - 1
c;         k = i
c;         p = d(i)
c;c
c;         do 260 j = ii, n
c;            if (d(j) .ge. p) go to 260
c;            k = j
c;            p = d(j)
c;  260    continue
c;c
c;         if (k .eq. i) go to 300
c;         d(k) = d(i)
c;         d(i) = p
c;c
c;         do 280 j = 1, n
c;            p = z(j,i)
c;            z(j,i) = z(j,k)
c;            z(j,k) = p
c;  280    continue
c;c
c;  300 continue
c;c
c;      go to 1001
c;c     .......... set error -- no convergence to an
c;c                eigenvalue after 30 iterations ..........
c; 1000 ierr = l
c; 1001 return
c;      end
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;      subroutine tqlrat(n,d,e2,ierr)
c;c
c;      integer i,j,l,m,n,l1,mml,ierr
c;      double precision d(n),e2(n)
c;      double precision b,c,f,g,h,p,r,s,t,epslon,pythag
c;c
c;c     this subroutine is a translation of the algol procedure tqlrat,
c;c     algorithm 464, comm. acm 16, 689(1973) by reinsch.
c;c
c;c     this subroutine finds the eigenvalues of a symmetric
c;c     tridiagonal matrix by the rational ql method.
c;c
c;c     on input
c;c
c;c        n is the order of the matrix.
c;c
c;c        d contains the diagonal elements of the input matrix.
c;c
c;c        e2 contains the squares of the subdiagonal elements of the
c;c          input matrix in its last n-1 positions.  e2(1) is arbitrary.
c;c
c;c      on output
c;c
c;c        d contains the eigenvalues in ascending order.  if an
c;c          error exit is made, the eigenvalues are correct and
c;c          ordered for indices 1,2,...ierr-1, but may not be
c;c          the smallest eigenvalues.
c;c
c;c        e2 has been destroyed.
c;c
c;c        ierr is set to
c;c          zero       for normal return,
c;c          j          if the j-th eigenvalue has not been
c;c                     determined after 30 iterations.
c;c
c;c     calls pythag for  dsqrt(a*a + b*b) .
c;c
c;c     questions and comments should be directed to burton s. garbow,
c;c     mathematics and computer science div, argonne national laboratory
c;c
c;c     this version dated august 1987.
c;c     modified by c. moler to fix underflow/overflow difficulties,
c;c     especially on the vax and other machines where epslon(1.0d0)**2
c;c     nearly underflows.  see the loop involving statement 102 and
c;c     the two statements just before statement 200.
c;c
c;c     ------------------------------------------------------------------
c;c
c;      ierr = 0
c;      if (n .eq. 1) go to 1001
c;c
c;      do 100 i = 2, n
c;  100 e2(i-1) = e2(i)
c;c
c;      f = 0.0d0
c;      t = 0.0d0
c;      e2(n) = 0.0d0
c;c
c;      do 290 l = 1, n
c;         j = 0
c;         h = dabs(d(l)) + dsqrt(e2(l))
c;         if (t .gt. h) go to 105
c;         t = h
c;         b = epslon(t)
c;         c = b * b
c;         if (c .ne. 0.0d0) go to 105
c;c        spliting tolerance underflowed.  look for larger value.
c;         do 102 i = l, n
c;            h = dabs(d(i)) + dsqrt(e2(i))
c;            if (h .gt. t) t = h
c;  102    continue
c;         b = epslon(t)
c;         c = b * b
c;c     .......... look for small squared sub-diagonal element ..........
c;  105    do 110 m = l, n
c;            if (e2(m) .le. c) go to 120
c;c     .......... e2(n) is always zero, so there is no exit
c;c                through the bottom of the loop ..........
c;  110    continue
c;c
c;  120    if (m .eq. l) go to 210
c;  130    if (j .eq. 30) go to 1000
c;         j = j + 1
c;c     .......... form shift ..........
c;         l1 = l + 1
c;         s = dsqrt(e2(l))
c;         g = d(l)
c;         p = (d(l1) - g) / (2.0d0 * s)
c;         r = pythag(p,1.0d0)
c;         d(l) = s / (p + dsign(r,p))
c;         h = g - d(l)
c;c
c;         do 140 i = l1, n
c;  140    d(i) = d(i) - h
c;c
c;         f = f + h
c;c     .......... rational ql transformation ..........
c;         g = d(m)
c;         if (g .eq. 0.0d0) g = b
c;         h = g
c;         s = 0.0d0
c;         mml = m - l
c;c     .......... for i=m-1 step -1 until l do -- ..........
c;c           do 200 ii = 1, mml
c;c           i = m - ii
c;            do 200 i=m-1,l,-1
c;            p = g * h
c;            r = p + e2(i)
c;            e2(i+1) = s * r
c;            s = e2(i) / r
c;            d(i+1) = h + s * (h + d(i))
c;            g = d(i) - e2(i) / g
c;c           avoid division by zero on next pass
c;            if (g .eq. 0.0d0) g = epslon(d(i))
c;            h = g * (p / r)
c;  200    continue
c;c
c;         e2(l) = s * g
c;         d(l) = h
c;c     .......... guard against underflow in convergence test ..........
c;         if (h .eq. 0.0d0) go to 210
c;         if (dabs(e2(l)) .le. dabs(c/h)) go to 210
c;         e2(l) = h * e2(l)
c;         if (e2(l) .ne. 0.0d0) go to 130
c;  210    p = d(l) + f
c;c     .......... order eigenvalues ..........
c;         if (l .eq. 1) go to 250
c;c     .......... for i=l step -1 until 2 do -- ..........
c;c        do 230 ii = 2, l
c;c           i = l + 2 - ii
c;            do 230 i=l,2,-1
c;            if (p .ge. d(i-1)) go to 270
c;            d(i) = d(i-1)
c;  230    continue
c;c
c;  250    i = 1
c;  270    d(i) = p
c;  290 continue
c;c
c;      go to 1001
c;c     .......... set error -- no convergence to an
c;c                eigenvalue after 30 iterations ..........
c; 1000 ierr = l
c; 1001 return
c;      end
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;      subroutine tred1(nm,n,a,d,e,e2)
c;c
c;      integer i,j,k,l,n,nm,jp1
c;      double precision a(nm,n),d(n),e(n),e2(n)
c;      double precision f,g,h,scale
c;c
c;c     this subroutine is a translation of the algol procedure tred1,
c;c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c;c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c;c
c;c     this subroutine reduces a real symmetric matrix
c;c     to a symmetric tridiagonal matrix using
c;c     orthogonal similarity transformations.
c;c
c;c     on input
c;c
c;c        nm must be set to the row dimension of two-dimensional
c;c          array parameters as declared in the calling program
c;c          dimension statement.
c;c
c;c        n is the order of the matrix.
c;c
c;c        a contains the real symmetric input matrix.  only the
c;c          lower triangle of the matrix need be supplied.
c;c
c;c     on output
c;c
c;c        a contains information about the orthogonal trans-
c;c          formations used in the reduction in its strict lower
c;c          triangle.  the full upper triangle of a is unaltered.
c;c
c;c        d contains the diagonal elements of the tridiagonal matrix.
c;c
c;c        e contains the subdiagonal elements of the tridiagonal
c;c          matrix in its last n-1 positions.  e(1) is set to zero.
c;c
c;c        e2 contains the squares of the corresponding elements of e.
c;c          e2 may coincide with e if the squares are not needed.
c;c
c;c     questions and comments should be directed to burton s. garbow,
c;c     mathematics and computer science div, argonne national laboratory
c;c
c;c     this version dated august 1983.
c;c
c;c     ------------------------------------------------------------------
c;c
c;      do 100 i = 1, n
c;         d(i) = a(n,i)
c;         a(n,i) = a(i,i)
c;  100 continue
c;c     .......... for i=n step -1 until 1 do -- ..........
c;c     do 300 ii = 1, n
c;c        i = n + 1 - ii
c;         do 300 i=n,1,-1
c;         l = i - 1
c;         h = 0.0d0
c;         scale = 0.0d0
c;         if (l .lt. 1) go to 130
c;c     .......... scale row (algol tol then not needed) ..........
c;         do 120 k = 1, l
c;  120    scale = scale + dabs(d(k))
c;c
c;         if (scale .ne. 0.0d0) go to 140
c;c
c;         do 125 j = 1, l
c;            d(j) = a(l,j)
c;            a(l,j) = a(i,j)
c;            a(i,j) = 0.0d0
c;  125    continue
c;c
c;  130    e(i) = 0.0d0
c;         e2(i) = 0.0d0
c;         go to 300
c;c
c;  140    do 150 k = 1, l
c;            d(k) = d(k) / scale
c;            h = h + d(k) * d(k)
c;  150    continue
c;c
c;         e2(i) = scale * scale * h
c;         f = d(l)
c;         g = -dsign(dsqrt(h),f)
c;         e(i) = scale * g
c;         h = h - f * g
c;         d(l) = f - g
c;         if (l .eq. 1) go to 285
c;c     .......... form a*u ..........
c;         do 170 j = 1, l
c;  170    e(j) = 0.0d0
c;c
c;         do 240 j = 1, l
c;            f = d(j)
c;            g = e(j) + a(j,j) * f
c;            jp1 = j + 1
c;            if (l .lt. jp1) go to 220
c;c
c;            do 200 k = jp1, l
c;               g = g + a(k,j) * d(k)
c;               e(k) = e(k) + a(k,j) * f
c;  200       continue
c;c
c;  220       e(j) = g
c;  240    continue
c;c     .......... form p ..........
c;         f = 0.0d0
c;c
c;         do 245 j = 1, l
c;            e(j) = e(j) / h
c;            f = f + e(j) * d(j)
c;  245    continue
c;c
c;         h = f / (h + h)
c;c     .......... form q ..........
c;         do 250 j = 1, l
c;  250    e(j) = e(j) - h * d(j)
c;c     .......... form reduced a ..........
c;         do 280 j = 1, l
c;            f = d(j)
c;            g = e(j)
c;c
c;            do 260 k = j, l
c;  260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
c;c
c;  280    continue
c;c
c;  285    do 290 j = 1, l
c;            f = d(j)
c;            d(j) = a(l,j)
c;            a(l,j) = a(i,j)
c;            a(i,j) = f * scale
c;  290    continue
c;c
c;  300 continue
c;c
c;      return
c;      end
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;      subroutine tred2(nm,n,a,d,e,z)
c;c
c;      integer i,j,k,l,n,nm,jp1
c;      double precision a(nm,n),d(n),e(n),z(nm,n)
c;      double precision f,g,h,hh,scale
c;c
c;c     this subroutine is a translation of the algol procedure tred2,
c;c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c;c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c;c
c;c     this subroutine reduces a real symmetric matrix to a
c;c     symmetric tridiagonal matrix using and accumulating
c;c     orthogonal similarity transformations.
c;c
c;c     on input
c;c
c;c        nm must be set to the row dimension of two-dimensional
c;c          array parameters as declared in the calling program
c;c          dimension statement.
c;c
c;c        n is the order of the matrix.
c;c
c;c        a contains the real symmetric input matrix.  only the
c;c          lower triangle of the matrix need be supplied.
c;c
c;c     on output
c;c
c;c        d contains the diagonal elements of the tridiagonal matrix.
c;c
c;c        e contains the subdiagonal elements of the tridiagonal
c;c          matrix in its last n-1 positions.  e(1) is set to zero.
c;c
c;c        z contains the orthogonal transformation matrix
c;c          produced in the reduction.
c;c
c;c        a and z may coincide.  if distinct, a is unaltered.
c;c
c;c     questions and comments should be directed to burton s. garbow,
c;c     mathematics and computer science div, argonne national laboratory
c;c
c;c     this version dated august 1983.
c;c
c;c     ------------------------------------------------------------------
c;c
c;      do 100 i = 1, n
c;c
c;         do 80 j = i, n
c;   80    z(j,i) = a(j,i)
c;c
c;         d(i) = a(n,i)
c;  100 continue
c;c
c;      if (n .eq. 1) go to 510
c;c     .......... for i=n step -1 until 2 do -- ..........
c;c     do 300 ii = 2, n
c;c        i = n + 2 - ii
c;      do 300 i=n,2,-1
c;         l = i - 1
c;         h = 0.0d0
c;         scale = 0.0d0
c;         if (l .lt. 2) go to 130
c;c     .......... scale row (algol tol then not needed) ..........
c;         do 120 k = 1, l
c;  120    scale = scale + dabs(d(k))
c;c
c;         if (scale .ne. 0.0d0) go to 140
c;  130    e(i) = d(l)
c;c
c;         do 135 j = 1, l
c;            d(j) = z(l,j)
c;            z(i,j) = 0.0d0
c;            z(j,i) = 0.0d0
c;  135    continue
c;c
c;         go to 290
c;c
c;  140    do 150 k = 1, l
c;            d(k) = d(k) / scale
c;            h = h + d(k) * d(k)
c;  150    continue
c;c
c;         f = d(l)
c;         g = -dsign(dsqrt(h),f)
c;         e(i) = scale * g
c;         h = h - f * g
c;         d(l) = f - g
c;c     .......... form a*u ..........
c;         do 170 j = 1, l
c;  170    e(j) = 0.0d0
c;c
c;         do 240 j = 1, l
c;            f = d(j)
c;            z(j,i) = f
c;            g = e(j) + z(j,j) * f
c;            jp1 = j + 1
c;            if (l .lt. jp1) go to 220
c;c
c;            do 200 k = jp1, l
c;               g = g + z(k,j) * d(k)
c;               e(k) = e(k) + z(k,j) * f
c;  200       continue
c;c
c;  220       e(j) = g
c;  240    continue
c;c     .......... form p ..........
c;         f = 0.0d0
c;c
c;         do 245 j = 1, l
c;            e(j) = e(j) / h
c;            f = f + e(j) * d(j)
c;  245    continue
c;c
c;         hh = f / (h + h)
c;c     .......... form q ..........
c;         do 250 j = 1, l
c;  250    e(j) = e(j) - hh * d(j)
c;c     .......... form reduced a ..........
c;         do 280 j = 1, l
c;            f = d(j)
c;            g = e(j)
c;c
c;            do 260 k = j, l
c;  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c;c
c;            d(j) = z(l,j)
c;            z(i,j) = 0.0d0
c;  280    continue
c;c
c;  290    d(i) = h
c;  300 continue
c;c     .......... accumulation of transformation matrices ..........
c;      do 500 i = 2, n
c;         l = i - 1
c;         z(n,l) = z(l,l)
c;         z(l,l) = 1.0d0
c;         h = d(i)
c;         if (h .eq. 0.0d0) go to 380
c;c
c;         do 330 k = 1, l
c;  330    d(k) = z(k,i) / h
c;c
c;         do 360 j = 1, l
c;            g = 0.0d0
c;c
c;            do 340 k = 1, l
c;  340       g = g + z(k,i) * z(k,j)
c;c
c;            do 360 k = 1, l
c;               z(k,j) = z(k,j) - g * d(k)
c;  360    continue
c;c
c;  380    do 400 k = 1, l
c;  400    z(k,i) = 0.0d0
c;c
c;  500 continue
c;c
c;  510 do 520 i = 1, n
c;         d(i) = z(n,i)
c;         z(n,i) = 0.0d0
c;  520 continue
c;c
c;      z(n,n) = 1.0d0
c;      e(1) = 0.0d0
c;      return
c;      end
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;      double precision function epslon (x)
c;      double precision x
c;c
c;c     estimate unit roundoff in quantities of size x.
c;c
c;      double precision a,b,c,eps
c;c
c;c     this program should function properly on all systems
c;c     satisfying the following two assumptions,
c;c        1.  the base used in representing floating point
c;c            numbers is not a power of three.
c;c        2.  the quantity  a  in statement 10 is represented to
c;c            the accuracy used in floating point variables
c;c            that are stored in memory.
c;c     the statement number 10 and the go to 10 are intended to
c;c     force optimizing compilers to generate code satisfying
c;c     assumption 2.
c;c     under these assumptions, it should be true that,
c;c            a  is not exactly equal to four-thirds,
c;c            b  has a zero for its last bit or digit,
c;c            c  is not exactly equal to one,
c;c            eps  measures the separation of 1.0 from
c;c                 the next larger floating point number.
c;c     the developers of eispack would appreciate being informed
c;c     about any systems where these assumptions do not hold.
c;c
c;c     this version dated 4/6/83.
c;c
c;      a = 4.0d0/3.0d0
c;   10 b = a - 1.0d0
c;      c = b + b + b
c;      eps = dabs(c-1.0d0)
c;      if (eps .eq. 0.0d0) go to 10
c;      epslon = eps*dabs(x)
c;      return
c;      end
c;      double precision function pythag(a,b)
c;      double precision a,b
c;c
c;c     finds dsqrt(a**2+b**2) without overflow or destructive underflow
c;c
c;      double precision p,r,s,t,u
c;      p = dmax1(dabs(a),dabs(b))
c;      if (p .eq. 0.0d0) go to 20
c;      r = (dmin1(dabs(a),dabs(b))/p)**2
c;   10 continue
c;         t = 4.0d0 + r
c;         if (t .eq. 4.0d0) go to 20
c;         s = r/t
c;         u = 1.0d0 + 2.0d0*s
c;         p = u*p
c;         r = (s/u)**2 * r
c;      go to 10
c;   20 pythag = p
c;      return
c;      end
#endif
***************************************************************************
*                                                                         *
*                          linpack library                                *
*                                                                         *
***************************************************************************
*                         routines included:                              *
*                                                                         *
*   1.  ssidi     computes determinant, inertia and inverse of a real     *
*                 symmetric matrix using the factors from ssifa           *
*                                                                         *
*   2.  ssifa     factors a real symmetric matrix by elimination          *
*                 with symmetric pivoting.                                *
*                                                                         *
*   3.  sgefa     factors a real matrix by gaussian elimination.          *
*                                                                         *
*   4.  sgesl     solves the real system a * x = b  or  trans(a) * x = b  *
*                 using the factors computed by sgeco or sgefa.           *
*                                                                         *
***************************************************************************
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;* ------------------------------------------------------------
c;      subroutine sgefa(a,lda,n,ipvt,info)
c;      implicit double precision (a-h,o-z)
c;      integer lda,n,ipvt(1),info
c;      dimension a(lda,1)
c;c      real a(lda,1)
c;c
c;c     sgefa factors a real matrix by gaussian elimination.
c;c
c;c     sgefa is usually called by sgeco, but it can be called
c;c     directly with a saving in time if  rcond  is not needed.
c;c     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
c;c
c;c     on entry
c;c
c;c        a       real(lda, n)
c;c                the matrix to be factored.
c;c
c;c        lda     integer
c;c                the leading dimension of the array  a .
c;c
c;c        n       integer
c;c                the order of the matrix  a .
c;c
c;c     on return
c;c
c;c        a       an upper triangular matrix and the multipliers
c;c                which were used to obtain it.
c;c                the factorization can be written  a = l*u  where
c;c                l  is a product of permutation and unit lower
c;c                triangular matrices and  u  is upper triangular.
c;c
c;c        ipvt    integer(n)
c;c                an integer vector of pivot indices.
c;c
c;c        info    integer
c;c                = 0  normal value.
c;c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c;c                     condition for this subroutine, but it does
c;c                     indicate that sgesl or sgedi will divide by zero
c;c                     if called.  use  rcond  in sgeco for a reliable
c;c                     indication of singularity.
c;c
c;c     linpack. this version dated 08/14/78 .
c;c     cleve moler, university of new mexico, argonne national lab.
c;c
c;c     subroutines and functions
c;c
c;c     blas saxpy,sscal,isamax
c;c
c;c     internal variables
c;c
c;c      real t
c;      integer idamax,j,k,kp1,l,nm1,idummy
c;      idummy=1
c;c
c;c
c;c     gaussian elimination with partial pivoting
c;c
c;      info = 0
c;      nm1 = n - 1
c;      if (nm1 .lt. 1) go to 70
c;      do 60 k = 1, nm1
c;         kp1 = k + 1
c;c
c;c        find l = pivot index
c;c
c;         l = idamax(n-k+1,a(k,k),1) + k - 1
c;         ipvt(k) = l
c;c
c;c        zero pivot implies this column already triangularized
c;c
c;         if (a(l,k) .eq. 0.0) go to 40
c;c
c;c           interchange if necessary
c;c
c;            if (l .eq. k) go to 10
c;               t = a(l,k)
c;               a(l,k) = a(k,k)
c;               a(k,k) = t
c;   10       continue
c;c
c;c           compute multipliers
c;c
c;            t = -1.0/a(k,k)
c;            call dscal(n-k,t,a(k+1,k),idummy)
c;c
c;c           row elimination with column indexing
c;c
c;            do 30 j = kp1, n
c;               t = a(l,j)
c;               if (l .eq. k) go to 20
c;                  a(l,j) = a(k,j)
c;                  a(k,j) = t
c;   20          continue
c;             call daxpy(n-k,t,a(k+1,k),idummy,a(k+1,j),idummy)
c;   30       continue
c;         go to 50
c;   40    continue
c;            info = k
c;   50    continue
c;   60 continue
c;   70 continue
c;      ipvt(n) = n
c;      if (a(n,n) .eq. 0.0) info = n
c;      return
c;      end
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;c ------------------------------------------------
c;      subroutine sgesl(a,lda,n,ipvt,b,job)
c;c current revision date: 29/09/87
c;      implicit double precision (a-h,o-z)
c;      integer lda,n,ipvt(1),job
c;      dimension a(lda,1),b(1)
c;c      real a(lda,1),b(1)
c;c
c;c     sgesl solves the real system
c;c     a * x = b  or  trans(a) * x = b
c;c     using the factors computed by sgeco or sgefa.
c;c
c;c     on entry
c;c
c;c        a       real(lda, n)
c;c                the output from sgeco or sgefa.
c;c
c;c        lda     integer
c;c                the leading dimension of the array  a .
c;c
c;c        n       integer
c;c                the order of the matrix  a .
c;c
c;c        ipvt    integer(n)
c;c                the pivot vector from sgeco or sgefa.
c;c
c;c        b       real(n)
c;c                the right hand side vector.
c;c
c;c        job     integer
c;c                = 0         to solve  a*x = b ,
c;c                = nonzero   to solve  trans(a)*x = b  where
c;c                            trans(a)  is the transpose.
c;c
c;c     on return
c;c
c;c        b       the solution vector  x .
c;c
c;c     error condition
c;c
c;c        a division by zero will occur if the input factor contains a
c;c        zero on the diagonal.  technically this indicates singularity
c;c        but it is often caused by improper arguments or improper
c;c        setting of lda .  it will not occur if the subroutines are
c;c        called correctly and if sgeco has set rcond .gt. 0.0
c;c        or sgefa has set info .eq. 0 .
c;c
c;c     to compute  inverse(a) * c  where  c  is a matrix
c;c     with  p  columns
c;c           call sgeco(a,lda,n,ipvt,rcond,z)
c;c           if (rcond is too small) go to ...
c;c           do 10 j = 1, p
c;c              call sgesl(a,lda,n,ipvt,c(1,j),0)
c;c        10 continue
c;c
c;c     linpack. this version dated 08/14/78 .
c;c     cleve moler, university of new mexico, argonne national lab.
c;c
c;c     subroutines and functions
c;c
c;c     blas saxpy,sdot
c;c
c;c     internal variables
c;c
c;c      real sdot,t
c;      integer k,kb,l,nm1,idummy
c;      idummy=1
c;c
c;      nm1 = n - 1
c;      if (job .ne. 0) go to 50
c;c
c;c        job = 0 , solve  a * x = b
c;c        first solve  l*y = b
c;c
c;         if (nm1 .lt. 1) go to 30
c;         do 20 k = 1, nm1
c;            l = ipvt(k)
c;            t = b(l)
c;            if (l .eq. k) go to 10
c;               b(l) = b(k)
c;               b(k) = t
c;   10       continue
c;            call daxpy(n-k,t,a(k+1,k),idummy,b(k+1),idummy)
c;   20    continue
c;   30    continue
c;c
c;c        now solve  u*x = y
c;c
c;         do 40 kb = 1, n
c;            k = n + 1 - kb
c;            b(k) = b(k)/a(k,k)
c;            t = -b(k)
c;            call daxpy(k-1,t,a(1,k),idummy,b(1),idummy)
c;   40    continue
c;      go to 100
c;   50 continue
c;c
c;c        job = nonzero, solve  trans(a) * x = b
c;c        first solve  trans(u)*y = b
c;c
c;         do 60 k = 1, n
c;            t = ddot(k-1,a(1,k),idummy,b(1),idummy)
c;            b(k) = (b(k) - t)/a(k,k)
c;   60    continue
c;c
c;c        now solve trans(l)*x = y
c;c
c;         if (nm1 .lt. 1) go to 90
c;         do 80 kb = 1, nm1
c;            k = n - kb
c;            b(k) = b(k) + ddot(n-k,a(k+1,k),idummy,b(k+1),idummy)
c;            l = ipvt(k)
c;            if (l .eq. k) go to 70
c;               t = b(l)
c;               b(l) = b(k)
c;               b(k) = t
c;   70       continue
c;   80    continue
c;   90    continue
c;  100 continue
c;      return
c;      end
#endif
* -----------------------------------------------------------------------
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
      subroutine syminv (a,lda,n,ierr)
      implicit double precision (a-h,o-z)
c
c     -----------------------------------------------------------------
c     This subroutine uses LAPACK DGETRF and DGETRI
c     to invert a real symmetric indefinite matrix.
c     latest revision:  18-feb-2008 by millard alexander (replace calls to "sy" with "ge"
c         lapack routines
c     -----------------------------------------------------------------
c
      dimension a(lda,n)
      dimension ipiv(4*n),work(256*n)

c
      lwork = 256*n
*  form full matrix
      do j = 2,n
         do i = 1,j-1
            a(i,j) = a(j,i)
         enddo
      enddo


*      write(6,444) lda,n
*444   format ('calling dgetrf:  lda=',i8,'  n=',i8)


      call dgetrf (n,n,a,lda,ipiv,ierr)
      if (ierr .ne. 0) stop 'error in dgetrf'
      call dgetri(n,a,lda,ipiv,work,lwork,ierr)
      if (ierr .ne. 0) stop 'error in dgetri'
      return
      end
#endif
#if (defined(HIB_UNIX) || defined(HIB_CRAY)) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;       subroutine smxinv (a, nmax, n, scr, kpvt, ierr)
c;*  subroutine to invert the real symmetric matrix a using
c;*  lapack routines dsytrf and dsytri
c;*  written by:  millard alexander
c;*  latest revision date:  13-nov-1995
c;* -----------------------------------------------------------------------
c;*  variables in call list:
c;*  nmax:   maximum row dimension of matrices
c;*  n:      actual order of matrices
c;*  a:      symmetric matrix of order n x n, stored in packed column form
c;*          on input:  lower triangle contains matrix to be inverted
c;*          on return: contains inverse(a), both lower and upper triangles
c;*  scr:    scratch vector of length n
c;*  kpvt:   scratch vector of length n
c;*  ierr:   on return:  set equal to 0 if normal return
c;*                      set equal to nn if singularity due to row nn of matri
c;*  subroutines used:
c;*  dsytrf, dsytri:   lapack routines to factor and to invert a symmetric, re
c;*                  matrix
c;* -----------------------------------------------------------------------
c;      implicit double  precision (a-h,o-z)
c;      integer icol, icolpt, ierr, ione, irowpt, izero, n, ncol,
c;     :        nmax, nmaxp1
c;      integer kpvt
c;      dimension a(1), scr(1), kpvt(1)
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;      common /cosc11/ sc11(1)
#endif
#if defined(HIB_UNIX_IBM)
c;      common /cokaux/ kaux
c;      dimension det(2)
#endif
#if (defined(HIB_UNIX) || defined(HIB_CRAY)) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;      data ione, izero  /1, 0/
#endif
#if defined(HIB_UNIX_IBM)
c;      naux=kaux
#endif
#if (defined(HIB_CRAY) || defined(HIB_UNIX)) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;      ierr = izero
c;*  first fill in upper half of original matrix
c;      nmaxp1 = nmax + 1
c;      icolpt = 2
c;      irowpt = nmaxp1
c;      do  50  icol = 1, n - 1
c;*  icolpt points to first sub-diagonal element in column icol
c;*  irowpt points to first super-diagonal element in row icol
c;*  ncol is number of subdiagonal elements in column icol
c;        ncol = n - icol
c;        call dcopy (ncol, a(icolpt), 1, a(irowpt), nmax)
c;        icolpt = icolpt + nmaxp1
c;        irowpt = irowpt + nmaxp1
c; 50   continue
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
c;       lwork=64*nmax
c;       lwork=-1
c;       call dsytrf('U',n,a,nmax,kpvt,sc11,lwork,info)
c;       print *, ' n, info, work(1):  ', n, info, sc11(1)
c;       stop 'after dsytrf'
c;       if (info .lt. 0) then
c;         write (9, 54) info
c;         write (6, 54) info
c;54       format (' *** ERROR IN',i4,
c;     :          'TH VARIABLE IN DSYTRF ***')
c;         stop
c;       endif
c;       if (info .gt. 0) then
c;         write (9, 55) info
c;         write (6, 55) info
c;55       format (' *** POSSIBLE SINGULARITY IN',i4,
c;     :          'TH ROW OF MATRIX TO BE INVERTED ***')
c;         stop
c;       end if
#endif
#if defined(HIB_UNIX_IBM)
c;      call dgeicd(a, nmax, n,izero,rcont,det,sc11,naux)
#endif
* old linpack calls
#if defined(HIB_CRAY)
c;      call ssifa (a, nmax, n, kpvt, ierr)
#endif
#if defined(HIB_UNIX_CONVEX)
c;      call dsifa (a, nmax, n, kpvt, ierr)
#endif
#if defined(HIB_CRAY) || defined(HIB_UNIX) || defined(HIB_UNIX_CONVEX)
c;      if (ierr .ne. 0) then
c;        write (9, 60) ierr
c;        write (6, 60) ierr
c;60      format (' *** POSSIBLE SINGULARITY IN',i4,
c;     :          'TH ROW OF MATRIX TO BE INVERTED ***')
c;        ierr = 0
c;      end if
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
c;       call dsytri('U',n,a,nmax,kpvt,sc11,info)
c;       if (info .lt. 0) then
c;         write (9, 64) info
c;         write (6, 64) info
c;64       format (' *** ERROR IN',i4,
c;     :          'TH VARIABLE IN DSYTRI ***')
c;         stop
c;       endif
c;       if (info .gt. 0) then
c;         write (9, 65) info
c;         write (6, 65) info
c;65       format (' *** POSSIBLE SINGULARITY IN',i4,
c;     :          'TH ROW OF MATRIX TO BE INVERTED ***')
c;         stop
c;       end if
#endif
#if defined(HIB_CRAY)
c;      call ssidi (a, nmax, n, kpvt, det, inert, scr, ione)
#endif
#if defined(HIB_UNIX_CONVEX)
c;      call dsidi (a, nmax, n, kpvt, det, inert, scr, ione)
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_CRAY) || defined(HIB_UNIX_CONVEX) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_IRIS)  || defined(HIB_UNIX_S)
c;*  inverse is now in upper triangle of matrix a
c;*  copy upper triangle of inverse back into lower triangle of a
c;      icolpt = 2
c;      irowpt = nmaxp1
c;      do 100 icol = 1, n - 1
c;*  icolpt points to first sub-diagonal element in column icol
c;*  irowpt points to first super-diagonal element in row icol
c;*  ncol is number of subdiagonal elements in column icol
c;        ncol = n - icol
#endif
#if defined(HIB_UNIX_CONVEX)
c;      call scopy (ncol, a(irowpt), nmax, a(icolpt), 1)
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_CRAY) || defined(HIB_UNIX_CONVEX) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_IRIS)  || defined(HIB_UNIX_S)
c;        call dcopy (ncol, a(irowpt), nmax, a(icolpt), 1)
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_CRAY) || defined(HIB_UNIX_CONVEX) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_IRIS)  || defined(HIB_UNIX_S)
c;        icolpt = icolpt + nmaxp1
c;        irowpt = irowpt + nmaxp1
c;100   continue
#endif
#if (defined(HIB_CRAY) || defined(HIB_UNIX)) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
c;        return
c;        end
#endif
*-----------------------------
#if defined(HIB_NONE)
c;      subroutine dsytri( uplo, n, a, lda, ipiv, work, info )
c;*
c;*  -- LAPACK routine (version 1.0b) --
c;*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c;*     Courant Institute, Argonne National Lab, and Rice University
c;*     October 31, 1992
c;*
c;*     .. Scalar Arguments ..
c;      character          uplo
c;      integer            info, lda, n
c;*     ..
c;*     .. Array Arguments ..
c;      integer            ipiv( * )
c;      double precision   a( lda, * ), work( * )
c;*     ..
c;*
c;*  Purpose
c;*  =======
c;*
c;*  DSYTRI computes the inverse of a real symmetric indefinite matrix
c;*  A using the factorization A = U*D*U' or A = L*D*L' computed by
c;*  DSYTRF.
c;*
c;*  Arguments
c;*  =========
c;*
c;*  UPLO    (input) CHARACTER*1
c;*          Specifies whether the details of the factorization are stored
c;*          as an upper or lower triangular matrix.
c;*          = 'U':  Upper triangular (form is A = U*D*U')
c;*          = 'L':  Lower triangular (form is A = L*D*L')
c;*
c;*  N       (input) INTEGER
c;*          The order of the matrix A.  N >= 0.
c;*
c;*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
c;*          On entry, the block diagonal matrix D and the multipliers
c;*          used to obtain the factor U or L as computed by DSYTRF.
c;*
c;*          On exit, if INFO = 0, the (symmetric) inverse of the original
c;*          matrix.  If UPLO = 'U', the upper triangular part of the
c;*          inverse is formed and the part of A below the diagonal is not
c;*          referenced; if UPLO = 'L' the lower triangular part of the
c;*          inverse is formed and the part of A above the diagonal is
c;*          not referenced.
c;*
c;*  LDA     (input) INTEGER
c;*          The leading dimension of the array A.  LDA >= max(1,N).
c;*
c;*  IPIV    (input) INTEGER array, dimension (N)
c;*          Details of the interchanges and the block structure of D
c;*          as determined by DSYTRF.
c;*
c;*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
c;*
c;*  INFO    (output) INTEGER
c;*          = 0: successful exit
c;*          < 0: if INFO = -k, the k-th argument had an illegal value
c;*          > 0: if INFO = k, D(k,k) = 0; the matrix is singular and its
c;*               inverse could not be computed.
c;*
c;*  =====================================================================
c;*
c;*     .. Parameters ..
c;      double precision   one, zero
c;      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
c;*     ..
c;*     .. Local Scalars ..
c;      logical            upper
c;      integer            k, kp, kstep
c;      double precision   ak, akkp1, akp1, d, t, temp
c;*     ..
c;*     .. External Functions ..
c;      logical            lsame
c;      double precision   ddot
c;      external           lsame, ddot
c;*     ..
c;*     .. External Subroutines ..
c;      external           dcopy, dswap, dsymv, xerbla
c;*     ..
c;*     .. Intrinsic Functions ..
c;      intrinsic          abs, max
c;*     ..
c;*     .. Executable Statements ..
c;*
c;*     Test the input parameters.
c;*
c;      info = 0
c;      upper = lsame( uplo, 'U' )
c;      if( .not.upper .and. .not.lsame( uplo, 'L' ) ) then
c;         info = -1
c;      else if( n.lt.0 ) then
c;         info = -2
c;      else if( lda.lt.max( 1, n ) ) then
c;         info = -4
c;      end if
c;      if( info.ne.0 ) then
c;         call xerbla( 'DSYTRI', -info )
c;         return
c;      end if
c;*
c;*     Quick return if possible
c;*
c;      if( n.eq.0 )
c;     $   return
c;*
c;*     Check that the diagonal matrix D is nonsingular.
c;*
c;      if( upper ) then
c;*
c;*        Upper triangular storage: examine D from bottom to top
c;*
c;         do 10 info = n, 1, -1
c;            if( ipiv( info ).gt.0 .and. a( info, info ).eq.zero )
c;     $         return
c;   10    continue
c;      else
c;*
c;*        Lower triangular storage: examine D from top to bottom.
c;*
c;         do 20 info = 1, n
c;            if( ipiv( info ).gt.0 .and. a( info, info ).eq.zero )
c;     $         return
c;   20    continue
c;      end if
c;      info = 0
c;*
c;      if( upper ) then
c;*
c;*        Compute inv(A) from the factorization A = U*D*U'.
c;*
c;*        K is the main loop index, increasing from 1 to N in steps of
c;*        1 or 2, depending on the size of the diagonal blocks.
c;*
c;         k = 1
c;   30    continue
c;*
c;*        If K > N, exit from loop.
c;*
c;         if( k.gt.n )
c;     $      go to 40
c;*
c;         if( ipiv( k ).gt.0 ) then
c;*
c;*           1 x 1 diagonal block
c;*
c;*           Invert the diagonal block.
c;*
c;            a( k, k ) = one / a( k, k )
c;*
c;*           Compute column K of the inverse.
c;*
c;            if( k.gt.1 ) then
c;               call dcopy( k-1, a( 1, k ), 1, work, 1 )
c;               call dsymv( uplo, k-1, -one, a, lda, work, 1, zero,
c;     $                     a( 1, k ), 1 )
c;               a( k, k ) = a( k, k ) - ddot( k-1, work, 1, a( 1, k ),
c;     $                     1 )
c;            end if
c;            kstep = 1
c;         else
c;*
c;*           2 x 2 diagonal block
c;*
c;*           Invert the diagonal block.
c;*
c;            t = abs( a( k, k+1 ) )
c;            ak = a( k, k ) / t
c;            akp1 = a( k+1, k+1 ) / t
c;            akkp1 = a( k, k+1 ) / t
c;            d = t*( ak*akp1-one )
c;            a( k, k ) = akp1 / d
c;            a( k+1, k+1 ) = ak / d
c;            a( k, k+1 ) = -akkp1 / d
c;*
c;*           Compute columns K and K+1 of the inverse.
c;*
c;            if( k.gt.1 ) then
c;               call dcopy( k-1, a( 1, k ), 1, work, 1 )
c;               call dsymv( uplo, k-1, -one, a, lda, work, 1, zero,
c;     $                     a( 1, k ), 1 )
c;               a( k, k ) = a( k, k ) - ddot( k-1, work, 1, a( 1, k ),
c;     $                     1 )
c;               a( k, k+1 ) = a( k, k+1 ) -
c;     $                       ddot( k-1, a( 1, k ), 1, a( 1, k+1 ), 1 )
c;               call dcopy( k-1, a( 1, k+1 ), 1, work, 1 )
c;               call dsymv( uplo, k-1, -one, a, lda, work, 1, zero,
c;     $                     a( 1, k+1 ), 1 )
c;               a( k+1, k+1 ) = a( k+1, k+1 ) -
c;     $                         ddot( k-1, work, 1, a( 1, k+1 ), 1 )
c;            end if
c;            kstep = 2
c;         end if
c;*
c;         kp = abs( ipiv( k ) )
c;         if( kp.ne.k ) then
c;*
c;*           Interchange rows and columns K and KP in the leading
c;*           submatrix A(1:k+1,1:k+1)
c;*
c;            call dswap( kp-1, a( 1, k ), 1, a( 1, kp ), 1 )
c;            call dswap( k-kp-1, a( kp+1, k ), 1, a( kp, kp+1 ), lda )
c;            temp = a( k, k )
c;            a( k, k ) = a( kp, kp )
c;            a( kp, kp ) = temp
c;            if( kstep.eq.2 ) then
c;               temp = a( k, k+1 )
c;               a( k, k+1 ) = a( kp, k+1 )
c;               a( kp, k+1 ) = temp
c;            end if
c;         end if
c;*
c;         k = k + kstep
c;         go to 30
c;   40    continue
c;*
c;      else
c;*
c;*        Compute inv(A) from the factorization A = L*D*L'.
c;*
c;*        K is the main loop index, increasing from 1 to N in steps of
c;*        1 or 2, depending on the size of the diagonal blocks.
c;*
c;         k = n
c;   50    continue
c;*
c;*        If K < 1, exit from loop.
c;*
c;         if( k.lt.1 )
c;     $      go to 60
c;*
c;         if( ipiv( k ).gt.0 ) then
c;*
c;*           1 x 1 diagonal block
c;*
c;*           Invert the diagonal block.
c;*
c;            a( k, k ) = one / a( k, k )
c;*
c;*           Compute column K of the inverse.
c;*
c;            if( k.lt.n ) then
c;               call dcopy( n-k, a( k+1, k ), 1, work, 1 )
c;               call dsymv( uplo, n-k, -one, a( k+1, k+1 ), lda, work, 1,
c;     $                     zero, a( k+1, k ), 1 )
c;               a( k, k ) = a( k, k ) - ddot( n-k, work, 1, a( k+1, k ),
c;     $                     1 )
c;            end if
c;            kstep = 1
c;         else
c;*
c;*           2 x 2 diagonal block
c;*
c;*           Invert the diagonal block.
c;*
c;            t = abs( a( k, k-1 ) )
c;            ak = a( k-1, k-1 ) / t
c;            akp1 = a( k, k ) / t
c;            akkp1 = a( k, k-1 ) / t
c;            d = t*( ak*akp1-one )
c;            a( k-1, k-1 ) = akp1 / d
c;            a( k, k ) = ak / d
c;            a( k, k-1 ) = -akkp1 / d
c;*
c;*           Compute columns K-1 and K of the inverse.
c;*
c;            if( k.lt.n ) then
c;               call dcopy( n-k, a( k+1, k ), 1, work, 1 )
c;               call dsymv( uplo, n-k, -one, a( k+1, k+1 ), lda, work, 1,
c;     $                     zero, a( k+1, k ), 1 )
c;               a( k, k ) = a( k, k ) - ddot( n-k, work, 1, a( k+1, k ),
c;     $                     1 )
c;               a( k, k-1 ) = a( k, k-1 ) -
c;     $                       ddot( n-k, a( k+1, k ), 1, a( k+1, k-1 ),
c;     $                       1 )
c;               call dcopy( n-k, a( k+1, k-1 ), 1, work, 1 )
c;               call dsymv( uplo, n-k, -one, a( k+1, k+1 ), lda, work, 1,
c;     $                     zero, a( k+1, k-1 ), 1 )
c;               a( k-1, k-1 ) = a( k-1, k-1 ) -
c;     $                         ddot( n-k, work, 1, a( k+1, k-1 ), 1 )
c;            end if
c;            kstep = 2
c;         end if
c;*
c;         kp = abs( ipiv( k ) )
c;         if( kp.ne.k ) then
c;*
c;*           Interchange rows and columns K and KP in the trailing
c;*           submatrix A(k-1:n,k-1:n)
c;*
c;            if( kp.lt.n )
c;     $         call dswap( n-kp, a( kp+1, k ), 1, a( kp+1, kp ), 1 )
c;            call dswap( kp-k-1, a( k+1, k ), 1, a( kp, k+1 ), lda )
c;            temp = a( k, k )
c;            a( k, k ) = a( kp, kp )
c;            a( kp, kp ) = temp
c;            if( kstep.eq.2 ) then
c;               temp = a( k, k-1 )
c;               a( k, k-1 ) = a( kp, k-1 )
c;               a( kp, k-1 ) = temp
c;            end if
c;         end if
c;*
c;         k = k - kstep
c;         go to 50
c;   60    continue
c;      end if
c;*
c;      return
c;*
c;*     End of DSYTRI
c;*
c;      end
c;      subroutine dsytrf( uplo, n, a, lda, ipiv, work, lwork, info )
c;*
c;*  -- LAPACK routine (version 1.0) --
c;*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c;*     Courant Institute, Argonne National Lab, and Rice University
c;*     February 29, 1992
c;*
c;*     .. Scalar Arguments ..
c;      character          uplo
c;      integer            info, lda, lwork, n
c;*     ..
c;*     .. Array Arguments ..
c;      integer            ipiv( * )
c;      double precision   a( lda, * ), work( lwork )
c;*     ..
c;*
c;*  Purpose
c;*  =======
c;*
c;*  DSYTRF computes the factorization of a real symmetric matrix A using
c;*  the Bunch-Kaufman diagonal pivoting method:
c;*
c;*     A = U*D*U'  or  A = L*D*L'
c;*
c;*  where U (or L) is a product of permutation and unit upper (lower)
c;*  triangular matrices, U' is the transpose of U, and D is symmetric and
c;*  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
c;*
c;*  This is the blocked version of the algorithm, calling Level 3 BLAS.
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
c;*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
c;*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
c;*          n-by-n upper triangular part of A contains the upper
c;*          triangular part of the matrix A, and the strictly lower
c;*          triangular part of A is not referenced.  If UPLO = 'L', the
c;*          leading n-by-n lower triangular part of A contains the lower
c;*          triangular part of the matrix A, and the strictly upper
c;*          triangular part of A is not referenced.
c;*
c;*          On exit, the block diagonal matrix D and the multipliers used
c;*          to obtain the factor U or L (see below for further details).
c;*
c;*  LDA     (input) INTEGER
c;*          The leading dimension of the array A.  LDA >= max(1,N).
c;*
c;*  IPIV    (output) INTEGER array, dimension (N)
c;*          Details of the interchanges and the block structure of D.
c;*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
c;*          interchanged and D(k,k) is a 1-by-1 diagonal block.
c;*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
c;*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
c;*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
c;*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
c;*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
c;*
c;*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
c;*          If INFO returns 0, then WORK(1) returns N*NB, the minimum
c;*          value of LWORK required to use the optimal blocksize.
c;*
c;*  LWORK   (input) INTEGER
c;*          The length of WORK.  LWORK should be >= N*NB, where NB is the
c;*          block size returned by ILAENV.
c;*
c;*  INFO    (output) INTEGER
c;*          = 0: successful exit
c;*          < 0: if INFO = -k, the k-th argument had an illegal value
c;*          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
c;*               has been completed, but the block diagonal matrix D is
c;*               exactly singular, and division by zero will occur if it
c;*               is used to solve a system of equations.
c;*
c;*  Further Details
c;*  ===============
c;*
c;*  If UPLO = 'U', then A = U*D*U', where
c;*     U = P(n)*U(n)* ... *P(k)U(k)* ...,
c;*  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
c;*  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
c;*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
c;*  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
c;*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
c;*
c;*             (   I    v    0   )   k-s
c;*     U(k) =  (   0    I    0   )   s
c;*             (   0    0    I   )   n-k
c;*                k-s   s   n-k
c;*
c;*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
c;*  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
c;*  and A(k,k), and v overwrites A(1:k-2,k-1:k).
c;*
c;*  If UPLO = 'L', then A = L*D*L', where
c;*     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
c;*  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
c;*  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
c;*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
c;*  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
c;*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
c;*
c;*             (   I    0     0   )  k-1
c;*     L(k) =  (   0    I     0   )  s
c;*             (   0    v     I   )  n-k-s+1
c;*                k-1   s  n-k-s+1
c;*
c;*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
c;*  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
c;*  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
c;*
c;*  =====================================================================
c;*
c;*     .. Local Scalars ..
c;      logical            upper
c;      integer            iinfo, iws, j, k, kb, ldwork, nb, nbmin
c;*     ..
c;*     .. External Functions ..
c;      logical            lsame
c;      integer            ilaenv
c;      external           lsame, ilaenv
c;*     ..
c;*     .. External Subroutines ..
c;      external           dlasyf, dsytf2, xerbla
c;*     ..
c;*     .. Intrinsic Functions ..
c;      intrinsic          max
c;*     ..
c;*     .. Executable Statements ..
c;*
c;*     Test the input parameters.
c;*
c;      info = 0
c;      upper = lsame( uplo, 'U' )
c;      if( .not.upper .and. .not.lsame( uplo, 'L' ) ) then
c;         info = -1
c;      else if( n.lt.0 ) then
c;         info = -2
c;      else if( lda.lt.max( 1, n ) ) then
c;         info = -4
c;      else if( lwork.lt.1 ) then
c;         info = -7
c;      end if
c;      if( info.ne.0 ) then
c;         call xerbla( 'DSYTRF', -info )
c;         return
c;      end if
c;*
c;*     Determine the block size
c;*
c;      nb = ilaenv( 1, 'DSYTRF', uplo, n, -1, -1, -1 )
c;      nbmin = 2
c;      ldwork = n
c;      if( nb.gt.1 .and. nb.lt.n ) then
c;         iws = ldwork*nb
c;         if( lwork.lt.iws ) then
c;            nb = max( lwork / ldwork, 1 )
c;            nbmin = max( 2, ilaenv( 2, 'DSYTRF', uplo, n, -1, -1, -1 ) )
c;         end if
c;      else
c;         iws = 1
c;      end if
c;      if( nb.lt.nbmin )
c;     $   nb = n
c;*
c;      if( upper ) then
c;*
c;*        Factorize A as U*D*U' using the upper triangle of A
c;*
c;*        K is the main loop index, decreasing from N to 1 in steps of
c;*        KB, where KB is the number of columns factorized by DLASYF;
c;*        KB is either NB or NB-1, or K for the last block
c;*
c;         k = n
c;   10    continue
c;*
c;*        If K < 1, exit from loop
c;*
c;         if( k.lt.1 )
c;     $      go to 40
c;*
c;         if( k.gt.nb ) then
c;*
c;*           Factorize columns k-kb+1:k of A and use blocked code to
c;*           update columns 1:k-kb
c;*
c;            call dlasyf( uplo, k, nb, kb, a, lda, ipiv, work, ldwork,
c;     $                   iinfo )
c;         else
c;*
c;*           Use unblocked code to factorize columns 1:k of A
c;*
c;            call dsytf2( uplo, k, a, lda, ipiv, iinfo )
c;            kb = k
c;         end if
c;*
c;*        Set INFO on the first occurrence of a zero pivot
c;*
c;         if( info.eq.0 .and. iinfo.gt.0 )
c;     $      info = iinfo
c;*
c;*        Decrease K and return to the start of the main loop
c;*
c;         k = k - kb
c;         go to 10
c;*
c;      else
c;*
c;*        Factorize A as L*D*L' using the lower triangle of A
c;*
c;*        K is the main loop index, increasing from 1 to N in steps of
c;*        KB, where KB is the number of columns factorized by DLASYF;
c;*        KB is either NB or NB-1, or N-K+1 for the last block
c;*
c;         k = 1
c;   20    continue
c;*
c;*        If K > N, exit from loop
c;*
c;         if( k.gt.n )
c;     $      go to 40
c;*
c;         if( k.le.n-nb ) then
c;*
c;*           Factorize columns k:k+kb-1 of A and use blocked code to
c;*           update columns k+kb:n
c;*
c;            call dlasyf( uplo, n-k+1, nb, kb, a( k, k ), lda, ipiv( k ),
c;     $                   work, ldwork, iinfo )
c;         else
c;*
c;*           Use unblocked code to factorize columns k:n of A
c;*
c;            call dsytf2( uplo, n-k+1, a( k, k ), lda, ipiv( k ), iinfo )
c;            kb = n - k + 1
c;         end if
c;*
c;*        Set INFO on the first occurrence of a zero pivot
c;*
c;         if( info.eq.0 .and. iinfo.gt.0 )
c;     $      info = iinfo + k - 1
c;*
c;*        Adjust IPIV
c;*
c;         do 30 j = k, k + kb - 1
c;            if( ipiv( j ).gt.0 ) then
c;               ipiv( j ) = ipiv( j ) + k - 1
c;            else
c;               ipiv( j ) = ipiv( j ) - k + 1
c;            end if
c;   30    continue
c;*
c;*        Increase K and return to the start of the main loop
c;*
c;         k = k + kb
c;         go to 20
c;*
c;      end if
c;*
c;   40 continue
c;      work( 1 ) = iws
c;      return
c;*
c;*     End of DSYTRF
c;*
c;      end
c;      subroutine dsytf2( uplo, n, a, lda, ipiv, info )
c;*
c;*  -- LAPACK routine (version 1.0b) --
c;*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
c;*     Courant Institute, Argonne National Lab, and Rice University
c;*     October 31, 1992
c;*
c;*     .. Scalar Arguments ..
c;      character          uplo
c;      integer            info, lda, n
c;*     ..
c;*     .. Array Arguments ..
c;      integer            ipiv( * )
c;      double precision   a( lda, * )
c;*     ..
c;*
c;*  Purpose
c;*  =======
c;*
c;*  DSYTF2 computes the factorization of a real symmetric matrix A using
c;*  the Bunch-Kaufman diagonal pivoting method:
c;*
c;*     A = U*D*U'  or  A = L*D*L'
c;*
c;*  where U (or L) is a product of permutation and unit upper (lower)
c;*  triangular matrices, U' is the transpose of U, and D is symmetric and
c;*  block diagonal with 1-by-1 and 2-by-2 diagonal blocks.
c;*
c;*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
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
c;*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
c;*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
c;*          n-by-n upper triangular part of A contains the upper
c;*          triangular part of the matrix A, and the strictly lower
c;*          triangular part of A is not referenced.  If UPLO = 'L', the
c;*          leading n-by-n lower triangular part of A contains the lower
c;*          triangular part of the matrix A, and the strictly upper
c;*          triangular part of A is not referenced.
c;*
c;*          On exit, the block diagonal matrix D and the multipliers used
c;*          to obtain the factor U or L (see below for further details).
c;*
c;*  LDA     (input) INTEGER
c;*          The leading dimension of the array A.  LDA >= max(1,N).
c;*
c;*  IPIV    (output) INTEGER array, dimension (N)
c;*          Details of the interchanges and the block structure of D.
c;*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
c;*          interchanged and D(k,k) is a 1-by-1 diagonal block.
c;*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
c;*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
c;*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
c;*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
c;*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
c;*
c;*  INFO    (output) INTEGER
c;*          = 0: successful exit
c;*          < 0: if INFO = -k, the k-th argument had an illegal value
c;*          > 0: if INFO = k, D(k,k) is exactly zero.  The factorization
c;*               has been completed, but the block diagonal matrix D is
c;*               exactly singular, and division by zero will occur if it
c;*               is used to solve a system of equations.
c;*
c;*  Further Details
c;*  ===============
c;*
c;*  If UPLO = 'U', then A = U*D*U', where
c;*     U = P(n)*U(n)* ... *P(k)U(k)* ...,
c;*  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
c;*  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
c;*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
c;*  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
c;*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
c;*
c;*             (   I    v    0   )   k-s
c;*     U(k) =  (   0    I    0   )   s
c;*             (   0    0    I   )   n-k
c;*                k-s   s   n-k
c;*
c;*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
c;*  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
c;*  and A(k,k), and v overwrites A(1:k-2,k-1:k).
c;*
c;*  If UPLO = 'L', then A = L*D*L', where
c;*     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
c;*  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
c;*  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
c;*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
c;*  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
c;*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
c;*
c;*             (   I    0     0   )  k-1
c;*     L(k) =  (   0    I     0   )  s
c;*             (   0    v     I   )  n-k-s+1
c;*                k-1   s  n-k-s+1
c;*
c;*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
c;*  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
c;*  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
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
c;      logical            upper
c;      integer            imax, jmax, k, kk, kp, kstep
c;      double precision   absakk, alpha, c, colmax, r1, r2, rowmax, s, t
c;*     ..
c;*     .. External Functions ..
c;      logical            lsame
c;      integer            idamax
c;      external           lsame, idamax
c;*     ..
c;*     .. External Subroutines ..
c;      external           dlaev2, drot, dscal, dswap, dsyr, xerbla
c;*     ..
c;*     .. Intrinsic Functions ..
c;      intrinsic          abs, max, sqrt
c;*     ..
c;*     .. Executable Statements ..
c;*
c;*     Test the input parameters.
c;*
c;      info = 0
c;      upper = lsame( uplo, 'U' )
c;      if( .not.upper .and. .not.lsame( uplo, 'L' ) ) then
c;         info = -1
c;      else if( n.lt.0 ) then
c;         info = -2
c;      else if( lda.lt.max( 1, n ) ) then
c;         info = -4
c;      end if
c;      if( info.ne.0 ) then
c;         call xerbla( 'DSYTF2', -info )
c;         return
c;      end if
c;*
c;*     Initialize ALPHA for use in choosing pivot block size.
c;*
c;      alpha = ( one+sqrt( sevten ) ) / eight
c;*
c;      if( upper ) then
c;*
c;*        Factorize A as U*D*U' using the upper triangle of A
c;*
c;*        K is the main loop index, decreasing from N to 1 in steps of
c;*        1 or 2
c;*
c;         k = n
c;   10    continue
c;*
c;*        If K < 1, exit from loop
c;*
c;         if( k.lt.1 )
c;     $      go to 30
c;         kstep = 1
c;*
c;*        Determine rows and columns to be interchanged and whether
c;*        a 1-by-1 or 2-by-2 pivot block will be used
c;*
c;         absakk = abs( a( k, k ) )
c;*
c;*        IMAX is the row-index of the largest off-diagonal element in
c;*        column K, and COLMAX is its absolute value
c;*
c;         if( k.gt.1 ) then
c;            imax = idamax( k-1, a( 1, k ), 1 )
c;            colmax = abs( a( imax, k ) )
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
c;*              JMAX is the column-index of the largest off-diagonal
c;*              element in row IMAX, and ROWMAX is its absolute value
c;*
c;               jmax = imax + idamax( k-imax, a( imax, imax+1 ), lda )
c;               rowmax = abs( a( imax, jmax ) )
c;               if( imax.gt.1 ) then
c;                  jmax = idamax( imax-1, a( 1, imax ), 1 )
c;                  rowmax = max( rowmax, abs( a( jmax, imax ) ) )
c;               end if
c;*
c;               if( absakk.ge.alpha*colmax*( colmax / rowmax ) ) then
c;*
c;*                 no interchange, use 1-by-1 pivot block
c;*
c;                  kp = k
c;               else if( abs( a( imax, imax ) ).ge.alpha*rowmax ) then
c;*
c;*                 interchange rows and columns K and IMAX, use 1-by-1
c;*                 pivot block
c;*
c;                  kp = imax
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
c;            if( kp.ne.kk ) then
c;*
c;*              Interchange rows and columns KK and KP in the leading
c;*              submatrix A(1:k,1:k)
c;*
c;               call dswap( kp-1, a( 1, kk ), 1, a( 1, kp ), 1 )
c;               call dswap( kk-kp-1, a( kp+1, kk ), 1, a( kp, kp+1 ),
c;     $                     lda )
c;               t = a( kk, kk )
c;               a( kk, kk ) = a( kp, kp )
c;               a( kp, kp ) = t
c;               if( kstep.eq.2 ) then
c;                  t = a( k-1, k )
c;                  a( k-1, k ) = a( kp, k )
c;                  a( kp, k ) = t
c;               end if
c;            end if
c;*
c;*           Update the leading submatrix
c;*
c;            if( kstep.eq.1 ) then
c;*
c;*              1-by-1 pivot block D(k): column k now holds
c;*
c;*              W(k) = U(k)*D(k)
c;*
c;*              where U(k) is the k-th column of U
c;*
c;*              Perform a rank-1 update of A(1:k-1,1:k-1) as
c;*
c;*              A := A - U(k)*D(k)*U(k)' = A - W(k)*1/D(k)*W(k)'
c;*
c;               r1 = one / a( k, k )
c;               call dsyr( uplo, k-1, -r1, a( 1, k ), 1, a, lda )
c;*
c;*              Store U(k) in column k
c;*
c;               call dscal( k-1, r1, a( 1, k ), 1 )
c;            else
c;*
c;*              2-by-2 pivot block D(k): columns k and k-1 now hold
c;*
c;*              ( W(k-1) W(k) ) = ( U(k-1) U(k) )*D(k)
c;*
c;*              where U(k) and U(k-1) are the k-th and (k-1)-th columns
c;*              of U
c;*
c;*              Perform a rank-2 update of A(1:k-2,1:k-2) as
c;*
c;*              A := A - ( U(k-1) U(k) )*D(k)*( U(k-1) U(k) )'
c;*                 = A - ( W(k-1) W(k) )*inv(D(k))*( W(k-1) W(k) )'
c;*
c;*              Convert this to two rank-1 updates by using the eigen-
c;*              decomposition of D(k)
c;*
c;               call dlaev2( a( k-1, k-1 ), a( k-1, k ), a( k, k ), r1,
c;     $                      r2, c, s )
c;               r1 = one / r1
c;               r2 = one / r2
c;               call drot( k-2, a( 1, k-1 ), 1, a( 1, k ), 1, c, s )
c;               call dsyr( uplo, k-2, -r1, a( 1, k-1 ), 1, a, lda )
c;               call dsyr( uplo, k-2, -r2, a( 1, k ), 1, a, lda )
c;*
c;*              Store U(k) and U(k-1) in columns k and k-1
c;*
c;               call dscal( k-2, r1, a( 1, k-1 ), 1 )
c;               call dscal( k-2, r2, a( 1, k ), 1 )
c;               call drot( k-2, a( 1, k-1 ), 1, a( 1, k ), 1, c, -s )
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
c;      else
c;*
c;*        Factorize A as L*D*L' using the lower triangle of A
c;*
c;*        K is the main loop index, increasing from 1 to N in steps of
c;*        1 or 2
c;*
c;         k = 1
c;   20    continue
c;*
c;*        If K > N, exit from loop
c;*
c;         if( k.gt.n )
c;     $      go to 30
c;         kstep = 1
c;*
c;*        Determine rows and columns to be interchanged and whether
c;*        a 1-by-1 or 2-by-2 pivot block will be used
c;*
c;         absakk = abs( a( k, k ) )
c;*
c;*        IMAX is the row-index of the largest off-diagonal element in
c;*        column K, and COLMAX is its absolute value
c;*
c;         if( k.lt.n ) then
c;            imax = k + idamax( n-k, a( k+1, k ), 1 )
c;            colmax = abs( a( imax, k ) )
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
c;*              JMAX is the column-index of the largest off-diagonal
c;*              element in row IMAX, and ROWMAX is its absolute value
c;*
c;               jmax = k - 1 + idamax( imax-k, a( imax, k ), lda )
c;               rowmax = abs( a( imax, jmax ) )
c;               if( imax.lt.n ) then
c;                  jmax = imax + idamax( n-imax, a( imax+1, imax ), 1 )
c;                  rowmax = max( rowmax, abs( a( jmax, imax ) ) )
c;               end if
c;*
c;               if( absakk.ge.alpha*colmax*( colmax / rowmax ) ) then
c;*
c;*                 no interchange, use 1-by-1 pivot block
c;*
c;                  kp = k
c;               else if( abs( a( imax, imax ) ).ge.alpha*rowmax ) then
c;*
c;*                 interchange rows and columns K and IMAX, use 1-by-1
c;*                 pivot block
c;*
c;                  kp = imax
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
c;            if( kp.ne.kk ) then
c;*
c;*              Interchange rows and columns KK and KP in the trailing
c;*              submatrix A(k:n,k:n)
c;*
c;               if( kp.lt.n )
c;     $            call dswap( n-kp, a( kp+1, kk ), 1, a( kp+1, kp ), 1 )
c;               call dswap( kp-kk-1, a( kk+1, kk ), 1, a( kp, kk+1 ),
c;     $                     lda )
c;               t = a( kk, kk )
c;               a( kk, kk ) = a( kp, kp )
c;               a( kp, kp ) = t
c;               if( kstep.eq.2 ) then
c;                  t = a( k+1, k )
c;                  a( k+1, k ) = a( kp, k )
c;                  a( kp, k ) = t
c;               end if
c;            end if
c;*
c;*           Update the trailing submatrix
c;*
c;            if( kstep.eq.1 ) then
c;*
c;*              1-by-1 pivot block D(k): column k now holds
c;*
c;*              W(k) = L(k)*D(k)
c;*
c;*              where L(k) is the k-th column of L
c;*
c;               if( k.lt.n ) then
c;*
c;*                 Perform a rank-1 update of A(k+1:n,k+1:n) as
c;*
c;*                 A := A - L(k)*D(k)*L(k)' = A - W(k)*(1/D(k))*W(k)'
c;*
c;                  r1 = one / a( k, k )
c;                  call dsyr( uplo, n-k, -r1, a( k+1, k ), 1,
c;     $                       a( k+1, k+1 ), lda )
c;*
c;*                 Store L(k) in column K
c;*
c;                  call dscal( n-k, r1, a( k+1, k ), 1 )
c;               end if
c;            else
c;*
c;*              2-by-2 pivot block D(k): columns K and K+1 now hold
c;*
c;*              ( W(k) W(k+1) ) = ( L(k) L(k+1) )*D(k)
c;*
c;*              where L(k) and L(k+1) are the k-th and (k+1)-th columns
c;*              of L
c;*
c;               if( k.lt.n-1 ) then
c;*
c;*                 Perform a rank-2 update of A(k+2:n,k+2:n) as
c;*
c;*                 A := A - ( L(k) L(k+1) )*D(k)*( L(k) L(k+1) )'
c;*                    = A - ( W(k) W(k+1) )*inv(D(k))*( W(k) W(k+1) )'
c;*
c;*                 Convert this to two rank-1 updates by using the eigen-
c;*                 decomposition of D(k)
c;*
c;                  call dlaev2( a( k, k ), a( k+1, k ), a( k+1, k+1 ),
c;     $                         r1, r2, c, s )
c;                  r1 = one / r1
c;                  r2 = one / r2
c;                  call drot( n-k-1, a( k+2, k ), 1, a( k+2, k+1 ), 1, c,
c;     $                       s )
c;                  call dsyr( uplo, n-k-1, -r1, a( k+2, k ), 1,
c;     $                       a( k+2, k+2 ), lda )
c;                  call dsyr( uplo, n-k-1, -r2, a( k+2, k+1 ), 1,
c;     $                       a( k+2, k+2 ), lda )
c;*
c;*                 Store L(k) and L(k+1) in columns k and k+1
c;*
c;                  call dscal( n-k-1, r1, a( k+2, k ), 1 )
c;                  call dscal( n-k-1, r2, a( k+2, k+1 ), 1 )
c;                  call drot( n-k-1, a( k+2, k ), 1, a( k+2, k+1 ), 1, c,
c;     $                       -s )
c;               end if
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
c;         go to 20
c;*
c;      end if
c;*
c;   30 continue
c;      return
c;*
c;*     End of DSYTF2
c;*
c;      end
#endif
#if !defined(HIB_CRAY) && !defined(HIB_UNIX_CONVEX)

c--------------------------------------------------------------------
      subroutine mxma(a,mcola,mrowa,b,mcolb,mrowb,
     1                r,mcolr,mrowr,ncol,nlink,nrow)
c--------------------------------------------------------------------
c.....unrolled version for ibm6000 and other risc machines
      implicit double precision (a-h,o-z)
      common /comxm/ ncache, mxmblk
      include "common/vax1"
      include "common/vax2"
      dimension r(*),a(*),b(*)
      if(ncol.eq.0.or.nrow.eq.0) return
      if(nlink.eq.0) then
        ijj=1
        do 6200 j=1,nrow
        ij=ijj
        do 6100 i=1,ncol
        r(ij)=0
6100    ij=ij+mcolr
6200    ijj=ijj+mrowr
        return
      else if(nlink.eq.1) then
        jj=1
        ijj=1
        do 6210 j=1,nrow
        ii=1
        ij=ijj
        do 6110 i=1,ncol
        r(ij)=a(ii)*b(jj)
        ii=ii+mcola
6110    ij=ij+mcolr
        ijj=ijj+mrowr
6210    jj=jj+mrowb
        return
      end if
      if(ncol.gt.5.or.nrow.gt.5) goto 6000
      goto (1000,2000,3000,4000,5000),nrow
c
c....nrow=1
1000  goto (1001,1002,1003,1004,1005),ncol
1001    ia1=1
        ib1=1
        s1=0
        do 1010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
1010    ib1=ib1+mcolb
        r(1)=s1
        return
c
1002    ia1=1
        ib1=1
        ia2=1+mcola
        s1=0
        s2=0
        do 1020 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
1020    ib1=ib1+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        return
c
1003    ia1=1
        ib1=1
        ia2=1+mcola
        ia3=ia2+mcola
        s1=0
        s2=0
        s3=0
        do 1030 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
1030    ib1=ib1+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        r(1+2*mcolr)=s3
        return
c
1004    ia1=1
        ib1=1
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        do 1040 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
        s4=a(ia4)*b(ib1)+s4
        ia4=ia4+mrowa
1040    ib1=ib1+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        r(1+2*mcolr)=s3
        r(1+3*mcolr)=s4
        return
c
1005    ia1=1
        ib1=1
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        do 1050 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
        s4=a(ia4)*b(ib1)+s4
        ia4=ia4+mrowa
        s5=a(ia5)*b(ib1)+s5
        ia5=ia5+mrowa
1050    ib1=ib1+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        r(1+2*mcolr)=s3
        r(1+3*mcolr)=s4
        r(1+4*mcolr)=s5
        return
c
c....nrow=2
2000  goto (2001,2002,2003,2004,2005),ncol
2001    ia1=1
        ib1=1
        ib2=1+mrowb
        s1=0
        s2=0
        do 2010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
2010    ia1=ia1+mrowa
        r(1)=s1
        r(1+mrowr)=s2
        return
c
2002    ia1=1
        ib1=1
        ia2=1+mcola
        ib2=1+mrowb
        s1=0
        s2=0
        s3=0
        s4=0
        do 2020 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia2)*b(ib1)+s2
        ib1=ib1+mcolb
        s3=a(ia1)*b(ib2)+s3
        ia1=ia1+mrowa
        s4=a(ia2)*b(ib2)+s4
        ib2=ib2+mcolb
2020    ia2=ia2+mrowa
        r(1)=s1
        r(1+mcolr)=s2
        r(1+mrowr)=s3
        r(1+mcolr+mrowr)=s4
        return
c
2003    ia1=1
        ib1=1
        ib2=1+mrowb
        ia2=1+mcola
        ia3=ia2+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        s6=0
        do 2030 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s4=a(ia1)*b(ib2)+s4
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        s5=a(ia2)*b(ib2)+s5
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ib1=ib1+mcolb
        s6=a(ia3)*b(ib2)+s6
        ia3=ia3+mrowa
2030    ib2=ib2+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        r(1+2*mcolr)=s3
        r(1+mrowr)=s4
        r(1+mrowr+mcolr)=s5
        r(1+mrowr+2*mcolr)=s6
        return
c
2004    ia1=1
        ib1=1
        ib2=1+mrowb
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        do 2040 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s12=a(ia1)*b(ib2)+s12
        ia1=ia1+mrowa
        s21=a(ia2)*b(ib1)+s21
        s22=a(ia2)*b(ib2)+s22
        ia2=ia2+mrowa
        s31=a(ia3)*b(ib1)+s31
        s32=a(ia3)*b(ib2)+s32
        ia3=ia3+mrowa
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s42=a(ia4)*b(ib2)+s42
        ia4=ia4+mrowa
2040    ib2=ib2+mcolb
        r(1)=s11
        r(1+mcolr)=s21
        r(1+2*mcolr)=s31
        r(1+3*mcolr)=s41
        r(1+mrowr)=s12
        r(1+mrowr+mcolr)=s22
        r(1+mrowr+2*mcolr)=s32
        r(1+mrowr+3*mcolr)=s42
        return
c
2005    ia1=1
        ib1=1
        ib2=1+mrowb
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s51=0
        s12=0
        s22=0
        s32=0
        s42=0
        s52=0
        do 2050 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s12=a(ia1)*b(ib2)+s12
        ia1=ia1+mrowa
        s21=a(ia2)*b(ib1)+s21
        s22=a(ia2)*b(ib2)+s22
        ia2=ia2+mrowa
        s31=a(ia3)*b(ib1)+s31
        s32=a(ia3)*b(ib2)+s32
        ia3=ia3+mrowa
        s41=a(ia4)*b(ib1)+s41
        s42=a(ia4)*b(ib2)+s42
        ia4=ia4+mrowa
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s52=a(ia5)*b(ib2)+s52
        ia5=ia5+mrowa
2050    ib2=ib2+mcolb
        r(1)=s11
        r(1+mcolr)=s21
        r(1+2*mcolr)=s31
        r(1+3*mcolr)=s41
        r(1+4*mcolr)=s51
        r(1+mrowr)=s12
        r(1+mrowr+mcolr)=s22
        r(1+mrowr+2*mcolr)=s32
        r(1+mrowr+3*mcolr)=s42
        r(1+mrowr+4*mcolr)=s52
        return
c
c....nrow=3
3000  goto (3001,3002,3003,3004,3005),ncol
3001    ia1=1
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        s1=0
        s2=0
        s3=0
        do 3010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
3010    ib3=ib3+mcolb
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        return
c
3002    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        s6=0
        do 3020 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s4=a(ia2)*b(ib1)+s4
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        s5=a(ia2)*b(ib2)+s5
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
        s6=a(ia2)*b(ib3)+s6
        ib3=ib3+mcolb
3020    ia2=ia2+mrowa
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+mcolr)=s4
        r(1+mcolr+mrowr)=s5
        r(1+mcolr+2*mrowr)=s6
        return
c
3003    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        s6=0
        s7=0
        s8=0
        s9=0
        do 3030 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
        s4=a(ia2)*b(ib1)+s4
        s5=a(ia2)*b(ib2)+s5
        s6=a(ia2)*b(ib3)+s6
        ia2=ia2+mrowa
        s7=a(ia3)*b(ib1)+s7
        ib1=ib1+mcolb
        s8=a(ia3)*b(ib2)+s8
        ib2=ib2+mcolb
        s9=a(ia3)*b(ib3)+s9
        ia3=ia3+mrowa
3030    ib3=ib3+mcolb
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+mcolr)=s4
        r(1+mcolr+mrowr)=s5
        r(1+mcolr+2*mrowr)=s6
        r(1+2*mcolr)=s7
        r(1+2*mcolr+mrowr)=s8
        r(1+2*mcolr+2*mrowr)=s9
        return
c
3004    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        s6=0
        s7=0
        s8=0
        s9=0
        s10=0
        s11=0
        s12=0
        do 3040 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
        s4=a(ia2)*b(ib1)+s4
        s5=a(ia2)*b(ib2)+s5
        s6=a(ia2)*b(ib3)+s6
        ia2=ia2+mrowa
        s7=a(ia3)*b(ib1)+s7
        s8=a(ia3)*b(ib2)+s8
        s9=a(ia3)*b(ib3)+s9
        ia3=ia3+mrowa
        s10=a(ia4)*b(ib1)+s10
        ib1=ib1+mcolb
        s11=a(ia4)*b(ib2)+s11
        ib2=ib2+mcolb
        s12=a(ia4)*b(ib3)+s12
        ia4=ia4+mrowa
3040    ib3=ib3+mcolb
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+mcolr)=s4
        r(1+mcolr+mrowr)=s5
        r(1+mcolr+2*mrowr)=s6
        r(1+2*mcolr)=s7
        r(1+2*mcolr+mrowr)=s8
        r(1+2*mcolr+2*mrowr)=s9
        r(1+3*mcolr)=s10
        r(1+3*mcolr+mrowr)=s11
        r(1+3*mcolr+2*mrowr)=s12
        return
c
3005    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s51=0
        s12=0
        s22=0
        s32=0
        s42=0
        s52=0
        s13=0
        s23=0
        s33=0
        s43=0
        s53=0
        do 3050 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s12=a(ia1)*b(ib2)+s12
        s13=a(ia1)*b(ib3)+s13
        ia1=ia1+mrowa
        s21=a(ia2)*b(ib1)+s21
        s22=a(ia2)*b(ib2)+s22
        s23=a(ia2)*b(ib3)+s23
        ia2=ia2+mrowa
        s31=a(ia3)*b(ib1)+s31
        s32=a(ia3)*b(ib2)+s32
        s33=a(ia3)*b(ib3)+s33
        ia3=ia3+mrowa
        s41=a(ia4)*b(ib1)+s41
        s42=a(ia4)*b(ib2)+s42
        s43=a(ia4)*b(ib3)+s43
        ia4=ia4+mrowa
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s52=a(ia5)*b(ib2)+s52
        ib2=ib2+mcolb
        s53=a(ia5)*b(ib3)+s53
        ia5=ia5+mrowa
3050    ib3=ib3+mcolb
        r(1)=s11
        r(1+mcolr)=s21
        r(1+2*mcolr)=s31
        r(1+3*mcolr)=s41
        r(1+4*mcolr)=s51
        r(1+mrowr)=s12
        r(1+mrowr+mcolr)=s22
        r(1+mrowr+2*mcolr)=s32
        r(1+mrowr+3*mcolr)=s42
        r(1+mrowr+4*mcolr)=s52
        r(1+2*mrowr)=s13
        r(1+2*mrowr+mcolr)=s23
        r(1+2*mrowr+2*mcolr)=s33
        r(1+2*mrowr+3*mcolr)=s43
        r(1+2*mrowr+4*mcolr)=s53
        return
c
c.....nrow=4
4000  goto(4001,4002,4003,4004,4005),ncol
4001    ia1=1
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        s1=0
        s2=0
        s3=0
        s4=0
        do 4010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ib3=ib3+mcolb
        s4=a(ia1)*b(ib4)+s4
        ib4=ib4+mcolb
4010    ia1=ia1+mrowa
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+3*mrowr)=s4
        return
c
4002    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ia1=1
        ia2=ia1+mcola
        s11=0
        s12=0
        s13=0
        s14=0
        s21=0
        s22=0
        s23=0
        s24=0
        do 4020 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        ia1=ia1+mrowa
        s24=a(ia2)*b(ib4)+s24
        ib4=ib4+mcolb
4020    ia2=ia2+mrowa
        r(1)=s11
        r(1+mrowr)=s12
        r(1+2*mrowr)=s13
        r(1+3*mrowr)=s14
        r(1+mcolr)=s21
        r(1+mcolr+mrowr)=s22
        r(1+mcolr+2*mrowr)=s23
        r(1+mcolr+3*mrowr)=s24
        return
c
4003    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0
        s12=0
        s13=0
        s14=0
        s21=0
        s22=0
        s23=0
        s24=0
        s31=0
        s32=0
        s33=0
        s34=0
        do 4030 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        ia1=ia1+mrowa
        s24=a(ia2)*b(ib4)+s24
        ia2=ia2+mrowa
        s34=a(ia3)*b(ib4)+s34
        ib4=ib4+mcolb
4030    ia3=ia3+mrowa
        r(1)=s11
        r(1+mrowr)=s12
        r(1+2*mrowr)=s13
        r(1+3*mrowr)=s14
        r(1+mcolr)=s21
        r(1+mcolr+mrowr)=s22
        r(1+mcolr+2*mrowr)=s23
        r(1+mcolr+3*mrowr)=s24
        r(1+2*mcolr)=s31
        r(1+2*mcolr+mrowr)=s32
        r(1+2*mcolr+2*mrowr)=s33
        r(1+2*mcolr+3*mrowr)=s34
        return
c
4004    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s12=0
        s13=0
        s14=0
        s21=0
        s22=0
        s23=0
        s24=0
        s31=0
        s32=0
        s33=0
        s34=0
        s41=0
        s42=0
        s43=0
        s44=0
        do 4040 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        ia1=ia1+mrowa
        s24=a(ia2)*b(ib4)+s24
        ia2=ia2+mrowa
        s34=a(ia3)*b(ib4)+s34
        ia3=ia3+mrowa
        s44=a(ia4)*b(ib4)+s44
        ib4=ib4+mcolb
4040    ia4=ia4+mrowa
        r(1)=s11
        r(1+mrowr)=s12
        r(1+2*mrowr)=s13
        r(1+3*mrowr)=s14
        r(1+mcolr)=s21
        r(1+mcolr+mrowr)=s22
        r(1+mcolr+2*mrowr)=s23
        r(1+mcolr+3*mrowr)=s24
        r(1+2*mcolr)=s31
        r(1+2*mcolr+mrowr)=s32
        r(1+2*mcolr+2*mrowr)=s33
        r(1+2*mcolr+3*mrowr)=s34
        r(1+3*mcolr)=s41
        r(1+3*mcolr+mrowr)=s42
        r(1+3*mcolr+2*mrowr)=s43
        r(1+3*mcolr+3*mrowr)=s44
        return
c
4005    ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        s1=0
        s2=0
        s3=0
        s4=0
        s6=0
        s7=0
        s8=0
        s9=0
        s11=0
        s12=0
        s13=0
        s14=0
        s16=0
        s17=0
        s18=0
        s19=0
        s21=0
        s22=0
        s23=0
        s24=0
        do 4050 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        s4=a(ia1)*b(ib4)+s4
        ia1=ia1+mrowa
        s6=a(ia2)*b(ib1)+s6
        s7=a(ia2)*b(ib2)+s7
        s8=a(ia2)*b(ib3)+s8
        s9=a(ia2)*b(ib4)+s9
        ia2=ia2+mrowa
        s11=a(ia3)*b(ib1)+s11
        s12=a(ia3)*b(ib2)+s12
        s13=a(ia3)*b(ib3)+s13
        s14=a(ia3)*b(ib4)+s14
        ia3=ia3+mrowa
        s16=a(ia4)*b(ib1)+s16
        s17=a(ia4)*b(ib2)+s17
        s18=a(ia4)*b(ib3)+s18
        s19=a(ia4)*b(ib4)+s19
        ia4=ia4+mrowa
        s21=a(ia5)*b(ib1)+s21
        s22=a(ia5)*b(ib2)+s22
        s23=a(ia5)*b(ib3)+s23
        s24=a(ia5)*b(ib4)+s24
        ia5=ia5+mrowa
        ib1=ib1+mcolb
        ib2=ib2+mcolb
        ib3=ib3+mcolb
        ib4=ib4+mcolb
4050    continue
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+3*mrowr)=s4
        r(1+mcolr)=s6
        r(1+mcolr+mrowr)=s7
        r(1+mcolr+2*mrowr)=s8
        r(1+mcolr+3*mrowr)=s9
        r(1+2*mcolr)=s11
        r(1+2*mcolr+mrowr)=s12
        r(1+2*mcolr+2*mrowr)=s13
        r(1+2*mcolr+3*mrowr)=s14
        r(1+3*mcolr)=s16
        r(1+3*mcolr+mrowr)=s17
        r(1+3*mcolr+2*mrowr)=s18
        r(1+3*mcolr+3*mrowr)=s19
        r(1+4*mcolr)=s21
        r(1+4*mcolr+mrowr)=s22
        r(1+4*mcolr+2*mrowr)=s23
        r(1+4*mcolr+3*mrowr)=s24
        return
c
c.....nrow=5
5000  goto(5001,5002,5003,5004,5004),ncol
5001    ia1=1
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        do 5010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ib3=ib3+mcolb
        s4=a(ia1)*b(ib4)+s4
        ib4=ib4+mcolb
        s5=a(ia1)*b(ib5)+s5
        ib5=ib5+mcolb
5010    ia1=ia1+mrowa
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+3*mrowr)=s4
        r(1+4*mrowr)=s5
        return
c
5002    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=1
        ia2=ia1+mcola
        s11=0
        s12=0
        s13=0
        s14=0
        s15=0
        s21=0
        s22=0
        s23=0
        s24=0
        s25=0
        do 5020 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        ia1=ia1+mrowa
        s25=a(ia2)*b(ib5)+s25
        ib5=ib5+mcolb
5020    ia2=ia2+mrowa
        r(1)=s11
        r(1+mrowr)=s12
        r(1+2*mrowr)=s13
        r(1+3*mrowr)=s14
        r(1+4*mrowr)=s15
        r(1+mcolr)=s21
        r(1+mcolr+mrowr)=s22
        r(1+mcolr+2*mrowr)=s23
        r(1+mcolr+3*mrowr)=s24
        r(1+mcolr+4*mrowr)=s25
        return
c
5003    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0
        s12=0
        s13=0
        s14=0
        s15=0
        s21=0
        s22=0
        s23=0
        s24=0
        s25=0
        s31=0
        s32=0
        s33=0
        s34=0
        s35=0
        do 5030 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        s34=a(ia3)*b(ib4)+s34
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        ia1=ia1+mrowa
        s25=a(ia2)*b(ib5)+s25
        ia2=ia2+mrowa
        s35=a(ia3)*b(ib5)+s35
        ib5=ib5+mcolb
5030    ia3=ia3+mrowa
        r(1)=s11
        r(1+mrowr)=s12
        r(1+2*mrowr)=s13
        r(1+3*mrowr)=s14
        r(1+4*mrowr)=s15
        r(1+mcolr)=s21
        r(1+mcolr+mrowr)=s22
        r(1+mcolr+2*mrowr)=s23
        r(1+mcolr+3*mrowr)=s24
        r(1+mcolr+4*mrowr)=s25
        r(1+2*mcolr)=s31
        r(1+2*mcolr+mrowr)=s32
        r(1+2*mcolr+2*mrowr)=s33
        r(1+2*mcolr+3*mrowr)=s34
        r(1+2*mcolr+4*mrowr)=s35
        return
c
5004    ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        s6=0
        s7=0
        s8=0
        s9=0
        s10=0
        s11=0
        s12=0
        s13=0
        s14=0
        s15=0
        s16=0
        s17=0
        s18=0
        s19=0
        s20=0
        s21=0
        s22=0
        s23=0
        s24=0
        s25=0
        do 5050 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        s4=a(ia1)*b(ib4)+s4
        s5=a(ia1)*b(ib5)+s5
        ia1=ia1+mrowa
        s6=a(ia2)*b(ib1)+s6
        s7=a(ia2)*b(ib2)+s7
        s8=a(ia2)*b(ib3)+s8
        s9=a(ia2)*b(ib4)+s9
        s10=a(ia2)*b(ib5)+s10
        ia2=ia2+mrowa
        s11=a(ia3)*b(ib1)+s11
        s12=a(ia3)*b(ib2)+s12
        s13=a(ia3)*b(ib3)+s13
        s14=a(ia3)*b(ib4)+s14
        s15=a(ia3)*b(ib5)+s15
        ia3=ia3+mrowa
        s16=a(ia4)*b(ib1)+s16
        s17=a(ia4)*b(ib2)+s17
        s18=a(ia4)*b(ib3)+s18
        s19=a(ia4)*b(ib4)+s19
        s20=a(ia4)*b(ib5)+s20
        ia4=ia4+mrowa
        if(ncol.eq.4) goto 5040
        s21=a(ia5)*b(ib1)+s21
        s22=a(ia5)*b(ib2)+s22
        s23=a(ia5)*b(ib3)+s23
        s24=a(ia5)*b(ib4)+s24
        s25=a(ia5)*b(ib5)+s25
        ia5=ia5+mrowa
5040    ib1=ib1+mcolb
        ib2=ib2+mcolb
        ib3=ib3+mcolb
        ib4=ib4+mcolb
5050    ib5=ib5+mcolb
        r(1)=s1
        r(1+mrowr)=s2
        r(1+2*mrowr)=s3
        r(1+3*mrowr)=s4
        r(1+4*mrowr)=s5
        r(1+mcolr)=s6
        r(1+mcolr+mrowr)=s7
        r(1+mcolr+2*mrowr)=s8
        r(1+mcolr+3*mrowr)=s9
        r(1+mcolr+4*mrowr)=s10
        r(1+2*mcolr)=s11
        r(1+2*mcolr+mrowr)=s12
        r(1+2*mcolr+2*mrowr)=s13
        r(1+2*mcolr+3*mrowr)=s14
        r(1+2*mcolr+4*mrowr)=s15
        r(1+3*mcolr)=s16
        r(1+3*mcolr+mrowr)=s17
        r(1+3*mcolr+2*mrowr)=s18
        r(1+3*mcolr+3*mrowr)=s19
        r(1+3*mcolr+4*mrowr)=s20
        if(ncol.eq.4) return
        r(1+4*mcolr)=s21
        r(1+4*mcolr+mrowr)=s22
        r(1+4*mcolr+2*mrowr)=s23
        r(1+4*mcolr+3*mrowr)=s24
        r(1+4*mcolr+4*mrowr)=s25
        return
c
6000  continue
      if(nlink*(ncol+nrow)+ncol*nrow.le.ncache) then
        call mxma4(a,mcola,mrowa,b,mcolb,mrowb,
     1             r,mcolr,mrowr,ncol,nlink,nrow)
        return
      end if
#endif
#if defined(HIB_UNIX_BLAS3)
c;      if(mcolr.eq.1) then
c;        if(mcola.eq.1.and.mcolb.eq.1) then
c;          call dgemm('N','N',ncol,nrow,nlink,1.0d0,a,mrowa,
c;     1                    b,max(nlink,mrowb),0.0d0,r,mrowr)
c;          return
c;        else if(mrowa.eq.1.and.mcolb.eq.1) then
c;          call dgemm('T','N',ncol,nrow,nlink,1.0d0,a,mcola,
c;     1                    b,max(nlink,mrowb),0.0d0,r,mrowr)
c;          return
c;        else if(mcola.eq.1.and.mrowb.eq.1) then
c;          call dgemm('N','T',ncol,nrow,nlink,1.0d0,a,mrowa,
c;     1                    b,mcolb,0.0d0,r,mrowr)
c;          return
c;        else if(mrowa.eq.1.and.mrowb.eq.1) then
c;          call dgemm('T','T',ncol,nrow,nlink,1.0d0,a,mcola,
c;     1                    b,mcolb,0.0d0,r,mrowr)
c;          return
c;        end if
c;        goto 6010
c;      else if(mrowr.eq.1) then
c;        if(mcola.eq.1.and.mcolb.eq.1) then
c;          call dgemm('T','T',nrow,ncol,nlink,1.0d0,b,max(nlink,mrowb),
c;     1                     a,mrowa,0.0d0,r,mcolr)
c;          return
c;        else if(mrowa.eq.1.and.mcolb.eq.1) then
c;          call dgemm('T','N',nrow,ncol,nlink,1.0d0,b,max(nlink,mrowb),
c;     1                     a,mcola,0.0d0,r,mcolr)
c;          return
c;        else if(mcola.eq.1.and.mrowb.eq.1) then
c;          call dgemm('N','T',nrow,ncol,nlink,1.0d0,b,mcolb,
c;     1                     a,mrowa,0.0d0,r,mcolr)
c;          return
c;        else if(mrowa.eq.1.and.mrowb.eq.1) then
c;          call dgemm('N','N',nrow,ncol,nlink,1.0d0,b,mcolb,
c;     1                     a,mcola,0.0d0,r,mcolr)
c;          return
c;        end if
c;      end if
#endif
#if !defined(HIB_CRAY) && !defined(HIB_UNIX_CONVEX)
6010  mxb=mxmblk
*      nkb=mxmbln
      nkb=108
      if(mrowa.ne.1.and.mcolb.ne.1) nkb=nkb/2
      ke=0
      if(ncol.gt.nrow) goto 6001
      do 60 k=1,nlink,nkb
      nr=nlink-ke
      nk=min(nkb,nr)
      if(nr.gt.nkb.and.nr.lt.nkb*2) nk=nr/2+1
      je=0
      do 70 j=1,nrow,mxb
      nr=nrow-je
      if(nr.eq.0) goto 70
      nj=min(mxb,nr)
      if(nr-nj.lt.4) nj=nr
      ie=0
      do 80 i=1,ncol,mxb
      nr=ncol-ie
      if(nr.eq.0) goto 80
      ni=min(mxb,nr)
      if(nr-ni.lt.4) ni=nr
      ia1=ke*mrowa+ie*mcola+1
      ib1=je*mrowb+ke*mcolb+1
      ir1=je*mrowr+ie*mcolr+1
      if(k.eq.1) then
        call mxma4(a(ia1),mcola,mrowa,b(ib1),mcolb,mrowb,r(ir1),
     1             mcolr, mrowr,ni,nk,nj)
      else
        call mxmb4(a(ia1),mcola,mrowa,b(ib1),mcolb,mrowb,r(ir1),
     1             mcolr,mrowr, ni,nk,nj)
      end if
80    ie=ie+ni
70    je=je+nj
60    ke=ke+nk
      return
6001  ke=0
      do 61 k=1,nlink,nkb
      nr=nlink-ke
      nk=min(nkb,nr)
      if(nr.gt.nkb.and.nr.lt.nkb*2) nk=nr/2+1
      ie=0
      do 81 i=1,ncol,mxb
      nr=ncol-ie
      if(nr.eq.0) goto 81
      ni=min(mxb,nr)
      if(nr-ni.lt.4) ni=nr
      je=0
      do 71 j=1,nrow,mxb
      nr=nrow-je
      if(nr.eq.0) goto 71
      nj=min(mxb,nr)
      if(nr-nj.lt.4) nj=nr
      ia1=ke*mrowa+ie*mcola+1
      ib1=je*mrowb+ke*mcolb+1
      ir1=je*mrowr+ie*mcolr+1
      if(k.eq.1) then
        call mxma4(a(ia1),mcola,mrowa,b(ib1),mcolb,mrowb,r(ir1),
     1             mcolr, mrowr,ni,nk,nj)
      else
        call mxmb4(a(ia1),mcola,mrowa,b(ib1),mcolb,mrowb,r(ir1),
     1             mcolr,mrowr, ni,nk,nj)
      end if
71    je=je+nj
81    ie=ie+ni
61    ke=ke+nk
      return
      end
c--------------------------------------------------------------------
      subroutine mxma4(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension r(1),a(1),b(1)
c... r(ncol,nrow)=a(ncol,nlink)*b(nlink,nrow) matrix mult
      iaa=1
      irr=1
      nrest=mod(nrow,4)
      if(nrest.eq.1.and.nrow.gt.1) nrest=5
      nre=nrow-nrest
      ncest=mod(ncol,4)
      if(ncest.eq.1.and.ncol.gt.1) ncest=5
      nce=ncol-ncest
      do 300 i=1,nce,4
      ibb=1
      ir1=irr
      ir2=ir1+mcolr
      ir3=ir2+mcolr
      ir4=ir3+mcolr
      irr=ir4+mcolr
      do 200 j=1,nre,4
      ib1=ibb
      ib2=ib1+mrowb
      ib3=ib2+mrowb
      ib4=ib3+mrowb
      ibb=ib4+mrowb
      ia1=iaa
      ia2=ia1+mcola
      ia3=ia2+mcola
      ia4=ia3+mcola
      s11=0
      s21=0
      s31=0
      s41=0
      s12=0
      s22=0
      s32=0
      s42=0
      s13=0
      s23=0
      s33=0
      s43=0
      s14=0
      s24=0
      s34=0
      s44=0
      do 100 k=1,nlink
      s11=a(ia1)*b(ib1)+s11
      s21=a(ia2)*b(ib1)+s21
      s31=a(ia3)*b(ib1)+s31
      s41=a(ia4)*b(ib1)+s41
      ib1=ib1+mcolb
      s12=a(ia1)*b(ib2)+s12
      s22=a(ia2)*b(ib2)+s22
      s32=a(ia3)*b(ib2)+s32
      s42=a(ia4)*b(ib2)+s42
      ib2=ib2+mcolb
      s13=a(ia1)*b(ib3)+s13
      s23=a(ia2)*b(ib3)+s23
      s33=a(ia3)*b(ib3)+s33
      s43=a(ia4)*b(ib3)+s43
      ib3=ib3+mcolb
      s14=a(ia1)*b(ib4)+s14
      s24=a(ia2)*b(ib4)+s24
      s34=a(ia3)*b(ib4)+s34
      s44=a(ia4)*b(ib4)+s44
      ib4=ib4+mcolb
      ia1=ia1+mrowa
      ia2=ia2+mrowa
      ia3=ia3+mrowa
100   ia4=ia4+mrowa
      r(ir1)=s11
      r(ir2)=s21
      r(ir3)=s31
      r(ir4)=s41
      r(ir1+mrowr)=s12
      r(ir2+mrowr)=s22
      r(ir3+mrowr)=s32
      r(ir4+mrowr)=s42
      r(ir1+2*mrowr)=s13
      r(ir2+2*mrowr)=s23
      r(ir3+2*mrowr)=s33
      r(ir4+2*mrowr)=s43
      r(ir1+3*mrowr)=s14
      r(ir2+3*mrowr)=s24
      r(ir3+3*mrowr)=s34
      r(ir4+3*mrowr)=s44
      ir1=ir1+4*mrowr
      ir2=ir2+4*mrowr
      ir3=ir3+4*mrowr
      ir4=ir4+4*mrowr
200   continue
      if(nrest.eq.0) goto 300
      goto (201,202,203,300,205),nrest
201     ib1=ibb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        do 101 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
101     ia4=ia4+mrowa
        r(ir1)=s11
        r(ir2)=s21
        r(ir3)=s31
        r(ir4)=s41
        goto 300
c
202     ib1=ibb
        ib2=ib1+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        do 102 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
102     ia4=ia4+mrowa
        r(ir1)=s11
        r(ir2)=s21
        r(ir3)=s31
        r(ir4)=s41
        r(ir1+mrowr)=s12
        r(ir2+mrowr)=s22
        r(ir3+mrowr)=s32
        r(ir4+mrowr)=s42
        goto 300
c
203     ib1=ibb
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        s13=0
        s23=0
        s33=0
        s43=0
        do 103 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
103     ia4=ia4+mrowa
        r(ir1)=s11
        r(ir2)=s21
        r(ir3)=s31
        r(ir4)=s41
        r(ir1+mrowr)=s12
        r(ir2+mrowr)=s22
        r(ir3+mrowr)=s32
        r(ir4+mrowr)=s42
        r(ir1+2*mrowr)=s13
        r(ir2+2*mrowr)=s23
        r(ir3+2*mrowr)=s33
        r(ir4+2*mrowr)=s43
        goto 300
c
205     ib1=ibb
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        s13=0
        s23=0
        s33=0
        s43=0
        s14=0
        s24=0
        s34=0
        s44=0
        s15=0
        s25=0
        s35=0
        s45=0
        do 105 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        s34=a(ia3)*b(ib4)+s34
        s44=a(ia4)*b(ib4)+s44
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        s25=a(ia2)*b(ib5)+s25
        s35=a(ia3)*b(ib5)+s35
        s45=a(ia4)*b(ib5)+s45
        ib5=ib5+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
105     ia4=ia4+mrowa
        r(ir1)=s11
        r(ir2)=s21
        r(ir3)=s31
        r(ir4)=s41
        r(ir1+mrowr)=s12
        r(ir2+mrowr)=s22
        r(ir3+mrowr)=s32
        r(ir4+mrowr)=s42
        r(ir1+2*mrowr)=s13
        r(ir2+2*mrowr)=s23
        r(ir3+2*mrowr)=s33
        r(ir4+2*mrowr)=s43
        r(ir1+3*mrowr)=s14
        r(ir2+3*mrowr)=s24
        r(ir3+3*mrowr)=s34
        r(ir4+3*mrowr)=s44
        r(ir1+4*mrowr)=s15
        r(ir2+4*mrowr)=s25
        r(ir3+4*mrowr)=s35
        r(ir4+4*mrowr)=s45
300   iaa=iaa+4*mcola
      if(ncest.eq.0) return
      goto (301,302,303,304,305),ncest
301     ibb=1
        ir1=irr
        do 2001 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        s11=0
        s12=0
        do 1001 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        ib2=ib2+mcolb
1001    ia1=ia1+mrowa
        r(ir1)=s11
        r(ir1+mrowr)=s12
        ir1=ir1+2*mrowr
2001    continue
        if(mod(nrow,2).eq.0) return
        s11=0
        do 1011 k=1,nlink
        s11=a(iaa)*b(ibb)+s11
        ibb=ibb+mcolb
1011    iaa=iaa+mrowa
        r(ir1)=s11
        return
c
302     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        do 2002 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        s11=0
        s21=0
        s12=0
        s22=0
        do 1002 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        ia1=ia1+mrowa
1002    ia2=ia2+mrowa
        r(ir1)=s11
        r(ir2)=s21
        r(ir1+mrowr)=s12
        r(ir2+mrowr)=s22
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
2002    continue
        if(mod(nrow,2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        s11=0
        s21=0
        do 1012 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        ia1=ia1+mrowa
        s21=a(ia2)*b(ibb)+s21
        ia2=ia2+mrowa
1012    ibb=ibb+mcolb
        r(ir1)=s11
        r(ir2)=s21
        return
c
303     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        ir3=ir2+mcolr
        do 2003 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0
        s21=0
        s31=0
        s12=0
        s22=0
        s32=0
        do 1003 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
1003    ia3=ia3+mrowa
        r(ir1)=s11
        r(ir2)=s21
        r(ir3)=s31
        r(ir1+mrowr)=s12
        r(ir2+mrowr)=s22
        r(ir3+mrowr)=s32
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
        ir3=ir3+2*mrowr
2003    continue
        if(mod(nrow,2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0
        s21=0
        s31=0
        do 1013 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        s21=a(ia2)*b(ibb)+s21
        s31=a(ia3)*b(ibb)+s31
        ibb=ibb+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
1013    ia3=ia3+mrowa
        r(ir1)=s11
        r(ir2)=s21
        r(ir3)=s31
304     return
305     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        ir3=ir2+mcolr
        ir4=ir3+mcolr
        ir5=ir4+mcolr
        do 2005 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s51=0
        s12=0
        s22=0
        s32=0
        s42=0
        s52=0
        do 1005 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        s52=a(ia5)*b(ib2)+s52
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
        ia4=ia4+mrowa
1005    ia5=ia5+mrowa
        r(ir1)=s11
        r(ir2)=s21
        r(ir3)=s31
        r(ir4)=s41
        r(ir5)=s51
        r(ir1+mrowr)=s12
        r(ir2+mrowr)=s22
        r(ir3+mrowr)=s32
        r(ir4+mrowr)=s42
        r(ir5+mrowr)=s52
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
        ir3=ir3+2*mrowr
        ir4=ir4+2*mrowr
        ir5=ir5+2*mrowr
2005    continue
        if(mod(nrow,2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s51=0
        do 1015 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        s21=a(ia2)*b(ibb)+s21
        s31=a(ia3)*b(ibb)+s31
        s41=a(ia4)*b(ibb)+s41
        s51=a(ia5)*b(ibb)+s51
        ibb=ibb+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
        ia4=ia4+mrowa
        ia5=ia5+mrowa
1015    continue
        r(ir1)=s11
        r(ir2)=s21
        r(ir3)=s31
        r(ir4)=s41
        r(ir5)=s51
        return
      end
c--------------------------------------------------------------------
      subroutine mxmb(a,mcola,mrowa,b,mcolb,mrowb,
     1                r,mcolr,mrowr,ncol,nlink,nrow)
c--------------------------------------------------------------------
c.....unrolled version for ibm6000 and other risc machines
      implicit double precision (a-h,o-z)
      common /comxm/ ncache, mxmblk
      include "common/vax1"
      dimension r(*),a(*),b(*)
      include "common/vax2"
      if(ncol.eq.0.or.nrow.eq.0.or.nlink.eq.0) return
      if(nlink.eq.1) then
        jj=1
        ijj=1
        do 6210 j=1,nrow
        ii=1
        ij=ijj
        do 6110 i=1,ncol
        r(ij)=r(ij)+a(ii)*b(jj)
        ii=ii+mcola
6110    ij=ij+mcolr
        ijj=ijj+mrowr
6210    jj=jj+mrowb
        return
      end if
      if(ncol.gt.5.or.nrow.gt.5) goto 6000
      goto (1000,2000,3000,4000,5000),nrow
c
c....nrow=1
1000  goto (1001,1002,1003,1004,1005),ncol
1001    ia1=1
        ib1=1
        s1=0
        do 1010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
1010    ib1=ib1+mcolb
        r(1)=r(1)+s1
        return
c
1002    ia1=1
        ib1=1
        ia2=1+mcola
        s1=0
        s2=0
        do 1020 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
1020    ib1=ib1+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        return
c
1003    ia1=1
        ib1=1
        ia2=1+mcola
        ia3=ia2+mcola
        s1=0
        s2=0
        s3=0
        do 1030 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
1030    ib1=ib1+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        r(1+2*mcolr)=r(1+2*mcolr)+s3
        return
c
1004    ia1=1
        ib1=1
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        do 1040 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
        s4=a(ia4)*b(ib1)+s4
        ia4=ia4+mrowa
1040    ib1=ib1+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        r(1+2*mcolr)=r(1+2*mcolr)+s3
        r(1+3*mcolr)=r(1+3*mcolr)+s4
        return
c
1005    ia1=1
        ib1=1
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        do 1050 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
        s4=a(ia4)*b(ib1)+s4
        ia4=ia4+mrowa
        s5=a(ia5)*b(ib1)+s5
        ia5=ia5+mrowa
1050    ib1=ib1+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        r(1+2*mcolr)=r(1+2*mcolr)+s3
        r(1+3*mcolr)=r(1+3*mcolr)+s4
        r(1+4*mcolr)=r(1+4*mcolr)+s5
        return
c
c....nrow=2
2000  goto (2001,2002,2003,2004,2005),ncol
2001    ia1=1
        ib1=1
        ib2=1+mrowb
        s1=0
        s2=0
        do 2010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
2010    ia1=ia1+mrowa
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        return
c
2002    ia1=1
        ib1=1
        ia2=1+mcola
        ib2=1+mrowb
        s1=0
        s2=0
        s3=0
        s4=0
        do 2020 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia2)*b(ib1)+s2
        ib1=ib1+mcolb
        s3=a(ia1)*b(ib2)+s3
        ia1=ia1+mrowa
        s4=a(ia2)*b(ib2)+s4
        ib2=ib2+mcolb
2020    ia2=ia2+mrowa
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        r(1+mrowr)=r(1+mrowr)+s3
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s4
        return
c
2003    ia1=1
        ib1=1
        ib2=1+mrowb
        ia2=1+mcola
        ia3=ia2+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        s6=0
        do 2030 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s4=a(ia1)*b(ib2)+s4
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        s5=a(ia2)*b(ib2)+s5
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ib1=ib1+mcolb
        s6=a(ia3)*b(ib2)+s6
        ia3=ia3+mrowa
2030    ib2=ib2+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        r(1+2*mcolr)=r(1+2*mcolr)+s3
        r(1+mrowr)=r(1+mrowr)+s4
        r(1+mrowr+mcolr)=r(1+mrowr+mcolr)+s5
        r(1+mrowr+2*mcolr)=r(1+mrowr+2*mcolr)+s6
        return
c
2004    ia1=1
        ib1=1
        ib2=1+mrowb
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        do 2040 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s12=a(ia1)*b(ib2)+s12
        ia1=ia1+mrowa
        s21=a(ia2)*b(ib1)+s21
        s22=a(ia2)*b(ib2)+s22
        ia2=ia2+mrowa
        s31=a(ia3)*b(ib1)+s31
        s32=a(ia3)*b(ib2)+s32
        ia3=ia3+mrowa
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s42=a(ia4)*b(ib2)+s42
        ia4=ia4+mrowa
2040    ib2=ib2+mcolb
        r(1)=r(1)+s11
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+2*mcolr)=r(1+2*mcolr)+s31
        r(1+3*mcolr)=r(1+3*mcolr)+s41
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+mrowr+mcolr)=r(1+mrowr+mcolr)+s22
        r(1+mrowr+2*mcolr)=r(1+mrowr+2*mcolr)+s32
        r(1+mrowr+3*mcolr)=r(1+mrowr+3*mcolr)+s42
        return
c
2005    ia1=1
        ib1=1
        ib2=1+mrowb
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s51=0
        s12=0
        s22=0
        s32=0
        s42=0
        s52=0
        do 2050 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s12=a(ia1)*b(ib2)+s12
        ia1=ia1+mrowa
        s21=a(ia2)*b(ib1)+s21
        s22=a(ia2)*b(ib2)+s22
        ia2=ia2+mrowa
        s31=a(ia3)*b(ib1)+s31
        s32=a(ia3)*b(ib2)+s32
        ia3=ia3+mrowa
        s41=a(ia4)*b(ib1)+s41
        s42=a(ia4)*b(ib2)+s42
        ia4=ia4+mrowa
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s52=a(ia5)*b(ib2)+s52
        ia5=ia5+mrowa
2050    ib2=ib2+mcolb
        r(1)=r(1)+s11
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+2*mcolr)=r(1+2*mcolr)+s31
        r(1+3*mcolr)=r(1+3*mcolr)+s41
        r(1+4*mcolr)=r(1+4*mcolr)+s51
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+mrowr+mcolr)=r(1+mrowr+mcolr)+s22
        r(1+mrowr+2*mcolr)=r(1+mrowr+2*mcolr)+s32
        r(1+mrowr+3*mcolr)=r(1+mrowr+3*mcolr)+s42
        r(1+mrowr+4*mcolr)=r(1+mrowr+4*mcolr)+s52
        return
c
c....nrow=3
3000  goto (3001,3002,3003,3004,3005),ncol
3001    ia1=1
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        s1=0
        s2=0
        s3=0
        do 3010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ib3=ib3+mcolb
3010    ia1=ia1+mrowa
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        return
c
3002    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        s6=0
        do 3020 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s4=a(ia2)*b(ib1)+s4
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        s5=a(ia2)*b(ib2)+s5
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
        s6=a(ia2)*b(ib3)+s6
        ib3=ib3+mcolb
3020    ia2=ia2+mrowa
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+mcolr)=r(1+mcolr)+s4
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s5
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s6
        return
c
3003    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        s6=0
        s7=0
        s8=0
        s9=0
        do 3030 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
        s4=a(ia2)*b(ib1)+s4
        s5=a(ia2)*b(ib2)+s5
        s6=a(ia2)*b(ib3)+s6
        ia2=ia2+mrowa
        s7=a(ia3)*b(ib1)+s7
        ib1=ib1+mcolb
        s8=a(ia3)*b(ib2)+s8
        ib2=ib2+mcolb
        s9=a(ia3)*b(ib3)+s9
        ia3=ia3+mrowa
3030    ib3=ib3+mcolb
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+mcolr)=r(1+mcolr)+s4
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s5
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s6
        r(1+2*mcolr)=r(1+2*mcolr)+s7
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s8
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s9
        return
c
3004    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        s6=0
        s7=0
        s8=0
        s9=0
        s10=0
        s11=0
        s12=0
        do 3040 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        ia1=ia1+mrowa
        s4=a(ia2)*b(ib1)+s4
        s5=a(ia2)*b(ib2)+s5
        s6=a(ia2)*b(ib3)+s6
        ia2=ia2+mrowa
        s7=a(ia3)*b(ib1)+s7
        s8=a(ia3)*b(ib2)+s8
        s9=a(ia3)*b(ib3)+s9
        ia3=ia3+mrowa
        s10=a(ia4)*b(ib1)+s10
        ib1=ib1+mcolb
        s11=a(ia4)*b(ib2)+s11
        ib2=ib2+mcolb
        s12=a(ia4)*b(ib3)+s12
        ia4=ia4+mrowa
3040    ib3=ib3+mcolb
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+mcolr)=r(1+mcolr)+s4
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s5
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s6
        r(1+2*mcolr)=r(1+2*mcolr)+s7
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s8
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s9
        r(1+3*mcolr)=r(1+3*mcolr)+s10
        r(1+3*mcolr+mrowr)=r(1+3*mcolr+mrowr)+s11
        r(1+3*mcolr+2*mrowr)=r(1+3*mcolr+2*mrowr)+s12
        return
c
3005    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s51=0
        s12=0
        s22=0
        s32=0
        s42=0
        s52=0
        s13=0
        s23=0
        s33=0
        s43=0
        s53=0
        do 3050 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s12=a(ia1)*b(ib2)+s12
        s13=a(ia1)*b(ib3)+s13
        ia1=ia1+mrowa
        s21=a(ia2)*b(ib1)+s21
        s22=a(ia2)*b(ib2)+s22
        s23=a(ia2)*b(ib3)+s23
        ia2=ia2+mrowa
        s31=a(ia3)*b(ib1)+s31
        s32=a(ia3)*b(ib2)+s32
        s33=a(ia3)*b(ib3)+s33
        ia3=ia3+mrowa
        s41=a(ia4)*b(ib1)+s41
        s42=a(ia4)*b(ib2)+s42
        s43=a(ia4)*b(ib3)+s43
        ia4=ia4+mrowa
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s52=a(ia5)*b(ib2)+s52
        ib2=ib2+mcolb
        s53=a(ia5)*b(ib3)+s53
        ia5=ia5+mrowa
3050    ib3=ib3+mcolb
        r(1)=r(1)+s11
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+2*mcolr)=r(1+2*mcolr)+s31
        r(1+3*mcolr)=r(1+3*mcolr)+s41
        r(1+4*mcolr)=r(1+4*mcolr)+s51
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+mrowr+mcolr)=r(1+mrowr+mcolr)+s22
        r(1+mrowr+2*mcolr)=r(1+mrowr+2*mcolr)+s32
        r(1+mrowr+3*mcolr)=r(1+mrowr+3*mcolr)+s42
        r(1+mrowr+4*mcolr)=r(1+mrowr+4*mcolr)+s52
        r(1+2*mrowr)=r(1+2*mrowr)+s13
        r(1+2*mrowr+mcolr)=r(1+2*mrowr+mcolr)+s23
        r(1+2*mrowr+2*mcolr)=r(1+2*mrowr+2*mcolr)+s33
        r(1+2*mrowr+3*mcolr)=r(1+2*mrowr+3*mcolr)+s43
        r(1+2*mrowr+4*mcolr)=r(1+2*mrowr+4*mcolr)+s53
        return
c
c.....nrow=4
4000  goto(4001,4002,4003,4004,4005),ncol
4001    ia1=1
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        s1=0
        s2=0
        s3=0
        s4=0
        do 4010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ib3=ib3+mcolb
        s4=a(ia1)*b(ib4)+s4
        ib4=ib4+mcolb
4010    ia1=ia1+mrowa
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+3*mrowr)=r(1+3*mrowr)+s4
        return
c
4002    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ia1=1
        ia2=ia1+mcola
        s11=0
        s12=0
        s13=0
        s14=0
        s21=0
        s22=0
        s23=0
        s24=0
        do 4020 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        ia1=ia1+mrowa
        s24=a(ia2)*b(ib4)+s24
        ib4=ib4+mcolb
4020    ia2=ia2+mrowa
        r(1)=r(1)+s11
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+2*mrowr)=r(1+2*mrowr)+s13
        r(1+3*mrowr)=r(1+3*mrowr)+s14
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s22
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s23
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s24
        return
c
4003    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0
        s12=0
        s13=0
        s14=0
        s21=0
        s22=0
        s23=0
        s24=0
        s31=0
        s32=0
        s33=0
        s34=0
        do 4030 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        ia1=ia1+mrowa
        s24=a(ia2)*b(ib4)+s24
        ia2=ia2+mrowa
        s34=a(ia3)*b(ib4)+s34
        ib4=ib4+mcolb
4030    ia3=ia3+mrowa
        r(1)=r(1)+s11
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+2*mrowr)=r(1+2*mrowr)+s13
        r(1+3*mrowr)=r(1+3*mrowr)+s14
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s22
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s23
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s24
        r(1+2*mcolr)=r(1+2*mcolr)+s31
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s32
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s33
        r(1+2*mcolr+3*mrowr)=r(1+2*mcolr+3*mrowr)+s34
        return
c
4004    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s12=0
        s13=0
        s14=0
        s21=0
        s22=0
        s23=0
        s24=0
        s31=0
        s32=0
        s33=0
        s34=0
        s41=0
        s42=0
        s43=0
        s44=0
        do 4040 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        ia1=ia1+mrowa
        s24=a(ia2)*b(ib4)+s24
        ia2=ia2+mrowa
        s34=a(ia3)*b(ib4)+s34
        ia3=ia3+mrowa
        s44=a(ia4)*b(ib4)+s44
        ib4=ib4+mcolb
4040    ia4=ia4+mrowa
        r(1)=r(1)+s11
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+2*mrowr)=r(1+2*mrowr)+s13
        r(1+3*mrowr)=r(1+3*mrowr)+s14
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s22
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s23
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s24
        r(1+2*mcolr)=r(1+2*mcolr)+s31
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s32
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s33
        r(1+2*mcolr+3*mrowr)=r(1+2*mcolr+3*mrowr)+s34
        r(1+3*mcolr)=r(1+3*mcolr)+s41
        r(1+3*mcolr+mrowr)=r(1+3*mcolr+mrowr)+s42
        r(1+3*mcolr+2*mrowr)=r(1+3*mcolr+2*mrowr)+s43
        r(1+3*mcolr+3*mrowr)=r(1+3*mcolr+3*mrowr)+s44
        return
c
4005    ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        s1=0
        s2=0
        s3=0
        s4=0
        s6=0
        s7=0
        s8=0
        s9=0
        s11=0
        s12=0
        s13=0
        s14=0
        s16=0
        s17=0
        s18=0
        s19=0
        s21=0
        s22=0
        s23=0
        s24=0
        do 4050 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        s4=a(ia1)*b(ib4)+s4
        ia1=ia1+mrowa
        s6=a(ia2)*b(ib1)+s6
        s7=a(ia2)*b(ib2)+s7
        s8=a(ia2)*b(ib3)+s8
        s9=a(ia2)*b(ib4)+s9
        ia2=ia2+mrowa
        s11=a(ia3)*b(ib1)+s11
        s12=a(ia3)*b(ib2)+s12
        s13=a(ia3)*b(ib3)+s13
        s14=a(ia3)*b(ib4)+s14
        ia3=ia3+mrowa
        s16=a(ia4)*b(ib1)+s16
        s17=a(ia4)*b(ib2)+s17
        s18=a(ia4)*b(ib3)+s18
        s19=a(ia4)*b(ib4)+s19
        ia4=ia4+mrowa
        s21=a(ia5)*b(ib1)+s21
        s22=a(ia5)*b(ib2)+s22
        s23=a(ia5)*b(ib3)+s23
        s24=a(ia5)*b(ib4)+s24
        ia5=ia5+mrowa
        ib1=ib1+mcolb
        ib2=ib2+mcolb
        ib3=ib3+mcolb
        ib4=ib4+mcolb
4050    continue
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+3*mrowr)=r(1+3*mrowr)+s4
        r(1+mcolr)=r(1+mcolr)+s6
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s7
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s8
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s9
        r(1+2*mcolr)=r(1+2*mcolr)+s11
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s12
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s13
        r(1+2*mcolr+3*mrowr)=r(1+2*mcolr+3*mrowr)+s14
        r(1+3*mcolr)=r(1+3*mcolr)+s16
        r(1+3*mcolr+mrowr)=r(1+3*mcolr+mrowr)+s17
        r(1+3*mcolr+2*mrowr)=r(1+3*mcolr+2*mrowr)+s18
        r(1+3*mcolr+3*mrowr)=r(1+3*mcolr+3*mrowr)+s19
        r(1+4*mcolr)=r(1+4*mcolr)+s21
        r(1+4*mcolr+mrowr)=r(1+4*mcolr+mrowr)+s22
        r(1+4*mcolr+2*mrowr)=r(1+4*mcolr+2*mrowr)+s23
        r(1+4*mcolr+3*mrowr)=r(1+4*mcolr+3*mrowr)+s24
        return
c
c.....nrow=5
5000  goto(5001,5002,5003,5004,5004),ncol
5001    ia1=1
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        do 5010 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ib1=ib1+mcolb
        s2=a(ia1)*b(ib2)+s2
        ib2=ib2+mcolb
        s3=a(ia1)*b(ib3)+s3
        ib3=ib3+mcolb
        s4=a(ia1)*b(ib4)+s4
        ib4=ib4+mcolb
        s5=a(ia1)*b(ib5)+s5
        ib5=ib5+mcolb
5010    ia1=ia1+mrowa
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+3*mrowr)=r(1+3*mrowr)+s4
        r(1+4*mrowr)=r(1+4*mrowr)+s5
        return
c
5002    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=1
        ia2=ia1+mcola
        s11=0
        s12=0
        s13=0
        s14=0
        s15=0
        s21=0
        s22=0
        s23=0
        s24=0
        s25=0
        do 5020 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        ia1=ia1+mrowa
        s25=a(ia2)*b(ib5)+s25
        ib5=ib5+mcolb
5020    ia2=ia2+mrowa
        r(1)=r(1)+s11
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+2*mrowr)=r(1+2*mrowr)+s13
        r(1+3*mrowr)=r(1+3*mrowr)+s14
        r(1+4*mrowr)=r(1+4*mrowr)+s15
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s22
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s23
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s24
        r(1+mcolr+4*mrowr)=r(1+mcolr+4*mrowr)+s25
        return
c
5003    ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0
        s12=0
        s13=0
        s14=0
        s15=0
        s21=0
        s22=0
        s23=0
        s24=0
        s25=0
        s31=0
        s32=0
        s33=0
        s34=0
        s35=0
        do 5030 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        s34=a(ia3)*b(ib4)+s34
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        ia1=ia1+mrowa
        s25=a(ia2)*b(ib5)+s25
        ia2=ia2+mrowa
        s35=a(ia3)*b(ib5)+s35
        ib5=ib5+mcolb
5030    ia3=ia3+mrowa
        r(1)=r(1)+s11
        r(1+mrowr)=r(1+mrowr)+s12
        r(1+2*mrowr)=r(1+2*mrowr)+s13
        r(1+3*mrowr)=r(1+3*mrowr)+s14
        r(1+4*mrowr)=r(1+4*mrowr)+s15
        r(1+mcolr)=r(1+mcolr)+s21
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s22
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s23
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s24
        r(1+mcolr+4*mrowr)=r(1+mcolr+4*mrowr)+s25
        r(1+2*mcolr)=r(1+2*mcolr)+s31
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s32
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s33
        r(1+2*mcolr+3*mrowr)=r(1+2*mcolr+3*mrowr)+s34
        r(1+2*mcolr+4*mrowr)=r(1+2*mcolr+4*mrowr)+s35
        return
c
5004    ia1=1
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        ib1=1
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        s1=0
        s2=0
        s3=0
        s4=0
        s5=0
        s6=0
        s7=0
        s8=0
        s9=0
        s10=0
        s11=0
        s12=0
        s13=0
        s14=0
        s15=0
        s16=0
        s17=0
        s18=0
        s19=0
        s20=0
        s21=0
        s22=0
        s23=0
        s24=0
        s25=0
        do 5050 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        s2=a(ia1)*b(ib2)+s2
        s3=a(ia1)*b(ib3)+s3
        s4=a(ia1)*b(ib4)+s4
        s5=a(ia1)*b(ib5)+s5
        ia1=ia1+mrowa
        s6=a(ia2)*b(ib1)+s6
        s7=a(ia2)*b(ib2)+s7
        s8=a(ia2)*b(ib3)+s8
        s9=a(ia2)*b(ib4)+s9
        s10=a(ia2)*b(ib5)+s10
        ia2=ia2+mrowa
        s11=a(ia3)*b(ib1)+s11
        s12=a(ia3)*b(ib2)+s12
        s13=a(ia3)*b(ib3)+s13
        s14=a(ia3)*b(ib4)+s14
        s15=a(ia3)*b(ib5)+s15
        ia3=ia3+mrowa
        s16=a(ia4)*b(ib1)+s16
        s17=a(ia4)*b(ib2)+s17
        s18=a(ia4)*b(ib3)+s18
        s19=a(ia4)*b(ib4)+s19
        s20=a(ia4)*b(ib5)+s20
        ia4=ia4+mrowa
        if(ncol.eq.4) goto 5040
        s21=a(ia5)*b(ib1)+s21
        s22=a(ia5)*b(ib2)+s22
        s23=a(ia5)*b(ib3)+s23
        s24=a(ia5)*b(ib4)+s24
        s25=a(ia5)*b(ib5)+s25
        ia5=ia5+mrowa
5040    ib1=ib1+mcolb
        ib2=ib2+mcolb
        ib3=ib3+mcolb
        ib4=ib4+mcolb
5050    ib5=ib5+mcolb
        r(1)=r(1)+s1
        r(1+mrowr)=r(1+mrowr)+s2
        r(1+2*mrowr)=r(1+2*mrowr)+s3
        r(1+3*mrowr)=r(1+3*mrowr)+s4
        r(1+4*mrowr)=r(1+4*mrowr)+s5
        r(1+mcolr)=r(1+mcolr)+s6
        r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s7
        r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s8
        r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s9
        r(1+mcolr+4*mrowr)=r(1+mcolr+4*mrowr)+s10
        r(1+2*mcolr)=r(1+2*mcolr)+s11
        r(1+2*mcolr+mrowr)=r(1+2*mcolr+mrowr)+s12
        r(1+2*mcolr+2*mrowr)=r(1+2*mcolr+2*mrowr)+s13
        r(1+2*mcolr+3*mrowr)=r(1+2*mcolr+3*mrowr)+s14
        r(1+2*mcolr+4*mrowr)=r(1+2*mcolr+4*mrowr)+s15
        r(1+3*mcolr)=r(1+3*mcolr)+s16
        r(1+3*mcolr+mrowr)=r(1+3*mcolr+mrowr)+s17
        r(1+3*mcolr+2*mrowr)=r(1+3*mcolr+2*mrowr)+s18
        r(1+3*mcolr+3*mrowr)=r(1+3*mcolr+3*mrowr)+s19
        r(1+3*mcolr+4*mrowr)=r(1+3*mcolr+4*mrowr)+s20
        if(ncol.eq.4) return
        r(1+4*mcolr)=r(1+4*mcolr)+s21
        r(1+4*mcolr+mrowr)=r(1+4*mcolr+mrowr)+s22
        r(1+4*mcolr+2*mrowr)=r(1+4*mcolr+2*mrowr)+s23
        r(1+4*mcolr+3*mrowr)=r(1+4*mcolr+3*mrowr)+s24
        r(1+4*mcolr+4*mrowr)=r(1+4*mcolr+4*mrowr)+s25
        return
c
6000  continue
      if(nlink*(ncol+nrow)+nrow*ncol.le.ncache) then
        call mxmb4(a,mcola,mrowa,b,mcolb,mrowb,
     1             r,mcolr,mrowr,ncol,nlink,nrow)
        return
      end if
#endif
#if defined(HIB_UNIX_BLAS3)
c;      if(mcolr.eq.1) then
c;        if(mcola.eq.1.and.mcolb.eq.1) then
c;          call dgemm('N','N',ncol,nrow,nlink,1.0d0,a,mrowa,
c;     1                     b,max(nlink,mrowb),1.0d0,r,mrowr)
c;          return
c;        else if(mrowa.eq.1.and.mcolb.eq.1) then
c;          call dgemm('T','N',ncol,nrow,nlink,1.0d0,a,mcola,
c;     1                     b,max(nlink,mrowb),1.0d0,r,mrowr)
c;          return
c;        else if(mcola.eq.1.and.mrowb.eq.1) then
c;          call dgemm('N','T',ncol,nrow,nlink,1.0d0,a,mrowa,
c;     1                     b,mcolb,1.0d0,r,mrowr)
c;          return
c;        else if(mrowa.eq.1.and.mrowb.eq.1) then
c;          call dgemm('T','T',ncol,nrow,nlink,1.0d0,a,mcola,
c;     1                     b,mcolb,1.0d0,r,mrowr)
c;          return
c;        end if
c;        goto 6010
c;      else if(mrowr.eq.1) then
c;        if(mcola.eq.1.and.mcolb.eq.1) then
c;          call dgemm('T','T',nrow,ncol,nlink,1.0d0,b,max(nlink,mrowb),
c;     1                     a,mrowa,1.0d0,r,mcolr)
c;          return
c;        else if(mrowa.eq.1.and.mcolb.eq.1) then
c;          call dgemm('T','N',nrow,ncol,nlink,1.0d0,b,max(nlink,mrowb),
c;     1                     a,mcola,1.0d0,r,mcolr)
c;          return
c;        else if(mcola.eq.1.and.mrowb.eq.1) then
c;          call dgemm('N','T',nrow,ncol,nlink,1.0d0,b,mcolb,
c;     1                     a,mrowa,1.0d0,r,mcolr)
c;          return
c;        else if(mrowa.eq.1.and.mrowb.eq.1) then
c;          call dgemm('N','N',nrow,ncol,nlink,1.0d0,b,mcolb,
c;     1                     a,mcola,1.0d0,r,mcolr)
c;          return
c;        end if
c;      end if
#endif
#if !defined(HIB_CRAY) && !defined(HIB_UNIX_CONVEX)

6010  mxb=mxmblk
*      nkb=mxmbln
      nkb=mxb
      if(mrowa.ne.1.and.mcolb.ne.1) nkb=nkb/2
      ke=0
      if(ncol.gt.nrow) goto 6001
      do 60 k=1,nlink,nkb
      nr=nlink-ke
      nk=min(nkb,nr)
      if(nr.gt.nkb.and.nr.lt.nkb*2) nk=nr/2+1
      je=0
      do 70 j=1,nrow,mxb
      nr=nrow-je
      if(nr.eq.0) goto 70
      nj=min(mxb,nr)
      if(nr-nj.lt.4) nj=nr
      ie=0
      do 80 i=1,ncol,mxb
      nr=ncol-ie
      if(nr.eq.0) goto 80
      ni=min(mxb,nr)
      if(nr-ni.lt.4) ni=nr
      ia1=ke*mrowa+ie*mcola+1
      ib1=je*mrowb+ke*mcolb+1
      ir1=je*mrowr+ie*mcolr+1
      call mxmb4(a(ia1),mcola,mrowa,b(ib1),mcolb,mrowb,r(ir1),
     1             mcolr,mrowr, ni,nk,nj)
80    ie=ie+ni
70    je=je+nj
60    ke=ke+nk
      return
6001  ke=0
      do 61 k=1,nlink,nkb
      nr=nlink-ke
      nk=min(nkb,nr)
      if(nr.gt.nkb.and.nr.lt.nkb*2) nk=nr/2+1
      ie=0
      do 81 i=1,ncol,mxb
      nr=ncol-ie
      if(nr.eq.0) goto 81
      ni=min(mxb,nr)
      if(nr-ni.lt.4) ni=nr
      je=0
      do 71 j=1,nrow,mxb
      nr=nrow-je
      if(nr.eq.0) goto 71
      nj=min(mxb,nr)
      if(nr-nj.lt.4) nj=nr
      ia1=ke*mrowa+ie*mcola+1
      ib1=je*mrowb+ke*mcolb+1
      ir1=je*mrowr+ie*mcolr+1
      call mxmb4(a(ia1),mcola,mrowa,b(ib1),mcolb,mrowb,r(ir1),
     1             mcolr,mrowr, ni,nk,nj)
71    je=je+nj
81    ie=ie+ni
61    ke=ke+nk
      return
      end
c--------------------------------------------------------------------
      subroutine mxmb4(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension r(1),a(1),b(1)
c... r(ncol,nrow)=a(ncol,nlink)*b(nlink,nrow) matrix mult
      iaa=1
      irr=1
      nrest=mod(nrow,4)
      if(nrest.eq.1.and.nrow.gt.1) nrest=5
      nre=nrow-nrest
      ncest=mod(ncol,4)
      if(ncest.eq.1.and.ncol.gt.1) ncest=5
      nce=ncol-ncest
      do 300 i=1,nce,4
      ibb=1
      ir1=irr
      ir2=ir1+mcolr
      ir3=ir2+mcolr
      ir4=ir3+mcolr
      irr=ir4+mcolr
      do 200 j=1,nre,4
      ib1=ibb
      ib2=ib1+mrowb
      ib3=ib2+mrowb
      ib4=ib3+mrowb
      ibb=ib4+mrowb
      ia1=iaa
      ia2=ia1+mcola
      ia3=ia2+mcola
      ia4=ia3+mcola
      s11=0
      s21=0
      s31=0
      s41=0
      s12=0
      s22=0
      s32=0
      s42=0
      s13=0
      s23=0
      s33=0
      s43=0
      s14=0
      s24=0
      s34=0
      s44=0
      do 100 k=1,nlink
      s11=a(ia1)*b(ib1)+s11
      s21=a(ia2)*b(ib1)+s21
      s31=a(ia3)*b(ib1)+s31
      s41=a(ia4)*b(ib1)+s41
      ib1=ib1+mcolb
      s12=a(ia1)*b(ib2)+s12
      s22=a(ia2)*b(ib2)+s22
      s32=a(ia3)*b(ib2)+s32
      s42=a(ia4)*b(ib2)+s42
      ib2=ib2+mcolb
      s13=a(ia1)*b(ib3)+s13
      s23=a(ia2)*b(ib3)+s23
      s33=a(ia3)*b(ib3)+s33
      s43=a(ia4)*b(ib3)+s43
      ib3=ib3+mcolb
      s14=a(ia1)*b(ib4)+s14
      s24=a(ia2)*b(ib4)+s24
      s34=a(ia3)*b(ib4)+s34
      s44=a(ia4)*b(ib4)+s44
      ib4=ib4+mcolb
      ia1=ia1+mrowa
      ia2=ia2+mrowa
      ia3=ia3+mrowa
100   ia4=ia4+mrowa
      r(ir1)=r(ir1)+s11
      r(ir2)=r(ir2)+s21
      r(ir3)=r(ir3)+s31
      r(ir4)=r(ir4)+s41
      r(ir1+mrowr)=r(ir1+mrowr)+s12
      r(ir2+mrowr)=r(ir2+mrowr)+s22
      r(ir3+mrowr)=r(ir3+mrowr)+s32
      r(ir4+mrowr)=r(ir4+mrowr)+s42
      r(ir1+2*mrowr)=r(ir1+2*mrowr)+s13
      r(ir2+2*mrowr)=r(ir2+2*mrowr)+s23
      r(ir3+2*mrowr)=r(ir3+2*mrowr)+s33
      r(ir4+2*mrowr)=r(ir4+2*mrowr)+s43
      r(ir1+3*mrowr)=r(ir1+3*mrowr)+s14
      r(ir2+3*mrowr)=r(ir2+3*mrowr)+s24
      r(ir3+3*mrowr)=r(ir3+3*mrowr)+s34
      r(ir4+3*mrowr)=r(ir4+3*mrowr)+s44
      ir1=ir1+4*mrowr
      ir2=ir2+4*mrowr
      ir3=ir3+4*mrowr
      ir4=ir4+4*mrowr
200   continue
      if(nrest.eq.0) goto 300
      goto (201,202,203,300,205),nrest
201     ib1=ibb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        do 101 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
101     ia4=ia4+mrowa
        r(ir1)=r(ir1)+s11
        r(ir2)=r(ir2)+s21
        r(ir3)=r(ir3)+s31
        r(ir4)=r(ir4)+s41
        goto 300
c
202     ib1=ibb
        ib2=ib1+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        do 102 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
102     ia4=ia4+mrowa
        r(ir1)=r(ir1)+s11
        r(ir2)=r(ir2)+s21
        r(ir3)=r(ir3)+s31
        r(ir4)=r(ir4)+s41
        r(ir1+mrowr)=r(ir1+mrowr)+s12
        r(ir2+mrowr)=r(ir2+mrowr)+s22
        r(ir3+mrowr)=r(ir3+mrowr)+s32
        r(ir4+mrowr)=r(ir4+mrowr)+s42
        goto 300
c
203     ib1=ibb
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        s13=0
        s23=0
        s33=0
        s43=0
        do 103 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
103     ia4=ia4+mrowa
        r(ir1)=r(ir1)+s11
        r(ir2)=r(ir2)+s21
        r(ir3)=r(ir3)+s31
        r(ir4)=r(ir4)+s41
        r(ir1+mrowr)=r(ir1+mrowr)+s12
        r(ir2+mrowr)=r(ir2+mrowr)+s22
        r(ir3+mrowr)=r(ir3+mrowr)+s32
        r(ir4+mrowr)=r(ir4+mrowr)+s42
        r(ir1+2*mrowr)=r(ir1+2*mrowr)+s13
        r(ir2+2*mrowr)=r(ir2+2*mrowr)+s23
        r(ir3+2*mrowr)=r(ir3+2*mrowr)+s33
        r(ir4+2*mrowr)=r(ir4+2*mrowr)+s43
        goto 300
c
205     ib1=ibb
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        s13=0
        s23=0
        s33=0
        s43=0
        s14=0
        s24=0
        s34=0
        s44=0
        s15=0
        s25=0
        s35=0
        s45=0
        do 105 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        s34=a(ia3)*b(ib4)+s34
        s44=a(ia4)*b(ib4)+s44
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        s25=a(ia2)*b(ib5)+s25
        s35=a(ia3)*b(ib5)+s35
        s45=a(ia4)*b(ib5)+s45
        ib5=ib5+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
105     ia4=ia4+mrowa
        r(ir1)=r(ir1)+s11
        r(ir2)=r(ir2)+s21
        r(ir3)=r(ir3)+s31
        r(ir4)=r(ir4)+s41
        r(ir1+mrowr)=r(ir1+mrowr)+s12
        r(ir2+mrowr)=r(ir2+mrowr)+s22
        r(ir3+mrowr)=r(ir3+mrowr)+s32
        r(ir4+mrowr)=r(ir4+mrowr)+s42
        r(ir1+2*mrowr)=r(ir1+2*mrowr)+s13
        r(ir2+2*mrowr)=r(ir2+2*mrowr)+s23
        r(ir3+2*mrowr)=r(ir3+2*mrowr)+s33
        r(ir4+2*mrowr)=r(ir4+2*mrowr)+s43
        r(ir1+3*mrowr)=r(ir1+3*mrowr)+s14
        r(ir2+3*mrowr)=r(ir2+3*mrowr)+s24
        r(ir3+3*mrowr)=r(ir3+3*mrowr)+s34
        r(ir4+3*mrowr)=r(ir4+3*mrowr)+s44
        r(ir1+4*mrowr)=r(ir1+4*mrowr)+s15
        r(ir2+4*mrowr)=r(ir2+4*mrowr)+s25
        r(ir3+4*mrowr)=r(ir3+4*mrowr)+s35
        r(ir4+4*mrowr)=r(ir4+4*mrowr)+s45
300   iaa=iaa+4*mcola
      if(ncest.eq.0) return
      goto (301,302,303,304,305),ncest
301     ibb=1
        ir1=irr
        do 2001 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        s11=0
        s12=0
        do 1001 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        ib2=ib2+mcolb
1001    ia1=ia1+mrowa
        r(ir1)=r(ir1)+s11
        r(ir1+mrowr)=r(ir1+mrowr)+s12
        ir1=ir1+2*mrowr
2001    continue
        if(mod(nrow,2).eq.0) return
        s11=0
        do 1011 k=1,nlink
        s11=a(iaa)*b(ibb)+s11
        ibb=ibb+mcolb
1011    iaa=iaa+mrowa
        r(ir1)=r(ir1)+s11
        return
c
302     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        do 2002 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        s11=0
        s21=0
        s12=0
        s22=0
        do 1002 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        ia1=ia1+mrowa
1002    ia2=ia2+mrowa
        r(ir1)=r(ir1)+s11
        r(ir2)=r(ir2)+s21
        r(ir1+mrowr)=r(ir1+mrowr)+s12
        r(ir2+mrowr)=r(ir2+mrowr)+s22
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
2002    continue
        if(mod(nrow,2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        s11=0
        s21=0
        do 1012 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        ia1=ia1+mrowa
        s21=a(ia2)*b(ibb)+s21
        ia2=ia2+mrowa
1012    ibb=ibb+mcolb
        r(ir1)=r(ir1)+s11
        r(ir2)=r(ir2)+s21
        return
c
303     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        ir3=ir2+mcolr
        do 2003 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0
        s21=0
        s31=0
        s12=0
        s22=0
        s32=0
        do 1003 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
1003    ia3=ia3+mrowa
        r(ir1)=r(ir1)+s11
        r(ir2)=r(ir2)+s21
        r(ir3)=r(ir3)+s31
        r(ir1+mrowr)=r(ir1+mrowr)+s12
        r(ir2+mrowr)=r(ir2+mrowr)+s22
        r(ir3+mrowr)=r(ir3+mrowr)+s32
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
        ir3=ir3+2*mrowr
2003    continue
        if(mod(nrow,2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0
        s21=0
        s31=0
        do 1013 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        s21=a(ia2)*b(ibb)+s21
        s31=a(ia3)*b(ibb)+s31
        ibb=ibb+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
1013    ia3=ia3+mrowa
        r(ir1)=r(ir1)+s11
        r(ir2)=r(ir2)+s21
        r(ir3)=r(ir3)+s31
304     return
305     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        ir3=ir2+mcolr
        ir4=ir3+mcolr
        ir5=ir4+mcolr
        do 2005 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s51=0
        s12=0
        s22=0
        s32=0
        s42=0
        s52=0
        do 1005 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        s52=a(ia5)*b(ib2)+s52
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
        ia4=ia4+mrowa
1005    ia5=ia5+mrowa
        r(ir1)=r(ir1)+s11
        r(ir2)=r(ir2)+s21
        r(ir3)=r(ir3)+s31
        r(ir4)=r(ir4)+s41
        r(ir5)=r(ir5)+s51
        r(ir1+mrowr)=r(ir1+mrowr)+s12
        r(ir2+mrowr)=r(ir2+mrowr)+s22
        r(ir3+mrowr)=r(ir3+mrowr)+s32
        r(ir4+mrowr)=r(ir4+mrowr)+s42
        r(ir5+mrowr)=r(ir5+mrowr)+s52
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
        ir3=ir3+2*mrowr
        ir4=ir4+2*mrowr
        ir5=ir5+2*mrowr
2005    continue
        if(mod(nrow,2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s51=0
        do 1015 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        s21=a(ia2)*b(ibb)+s21
        s31=a(ia3)*b(ibb)+s31
        s41=a(ia4)*b(ibb)+s41
        s51=a(ia5)*b(ibb)+s51
        ibb=ibb+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
        ia4=ia4+mrowa
        ia5=ia5+mrowa
1015    continue
        r(ir1)=r(ir1)+s11
        r(ir2)=r(ir2)+s21
        r(ir3)=r(ir3)+s31
        r(ir4)=r(ir4)+s41
        r(ir5)=r(ir5)+s51
        return
      end
c--------------------------------------------------------------------
      subroutine mxman(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension r(1),a(1),b(1)
c... r(ncol,nrow)=-a(ncol,nlink)*b(nlink,nrow) matrix mult
      iaa=1
      irr=1
      nrest=mod(nrow,4)
      if(nrest.eq.1.and.nrow.gt.1) nrest=5
      nre=nrow-nrest
      ncest=mod(ncol,4)
      if(ncest.eq.1.and.ncol.gt.1) ncest=5
      nce=ncol-ncest
      do 300 i=1,nce,4
      ibb=1
      ir1=irr
      ir2=ir1+mcolr
      ir3=ir2+mcolr
      ir4=ir3+mcolr
      irr=ir4+mcolr
      do 200 j=1,nre,4
      ib1=ibb
      ib2=ib1+mrowb
      ib3=ib2+mrowb
      ib4=ib3+mrowb
      ibb=ib4+mrowb
      ia1=iaa
      ia2=ia1+mcola
      ia3=ia2+mcola
      ia4=ia3+mcola
      s11=0
      s21=0
      s31=0
      s41=0
      s12=0
      s22=0
      s32=0
      s42=0
      s13=0
      s23=0
      s33=0
      s43=0
      s14=0
      s24=0
      s34=0
      s44=0
      do 100 k=1,nlink
      s11=a(ia1)*b(ib1)+s11
      s21=a(ia2)*b(ib1)+s21
      s31=a(ia3)*b(ib1)+s31
      s41=a(ia4)*b(ib1)+s41
      ib1=ib1+mcolb
      s12=a(ia1)*b(ib2)+s12
      s22=a(ia2)*b(ib2)+s22
      s32=a(ia3)*b(ib2)+s32
      s42=a(ia4)*b(ib2)+s42
      ib2=ib2+mcolb
      s13=a(ia1)*b(ib3)+s13
      s23=a(ia2)*b(ib3)+s23
      s33=a(ia3)*b(ib3)+s33
      s43=a(ia4)*b(ib3)+s43
      ib3=ib3+mcolb
      s14=a(ia1)*b(ib4)+s14
      s24=a(ia2)*b(ib4)+s24
      s34=a(ia3)*b(ib4)+s34
      s44=a(ia4)*b(ib4)+s44
      ib4=ib4+mcolb
      ia1=ia1+mrowa
      ia2=ia2+mrowa
      ia3=ia3+mrowa
100   ia4=ia4+mrowa
      r(ir1)=-s11
      r(ir2)=-s21
      r(ir3)=-s31
      r(ir4)=-s41
      r(ir1+mrowr)=-s12
      r(ir2+mrowr)=-s22
      r(ir3+mrowr)=-s32
      r(ir4+mrowr)=-s42
      r(ir1+2*mrowr)=-s13
      r(ir2+2*mrowr)=-s23
      r(ir3+2*mrowr)=-s33
      r(ir4+2*mrowr)=-s43
      r(ir1+3*mrowr)=-s14
      r(ir2+3*mrowr)=-s24
      r(ir3+3*mrowr)=-s34
      r(ir4+3*mrowr)=-s44
      ir1=ir1+4*mrowr
      ir2=ir2+4*mrowr
      ir3=ir3+4*mrowr
      ir4=ir4+4*mrowr
200   continue
      if(nrest.eq.0) goto 300
      goto (201,202,203,300,205),nrest
201     ib1=ibb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        do 101 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
101     ia4=ia4+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir4)=-s41
        goto 300
c
202     ib1=ibb
        ib2=ib1+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        do 102 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
102     ia4=ia4+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir4)=-s41
        r(ir1+mrowr)=-s12
        r(ir2+mrowr)=-s22
        r(ir3+mrowr)=-s32
        r(ir4+mrowr)=-s42
        goto 300
c
203     ib1=ibb
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        s13=0
        s23=0
        s33=0
        s43=0
        do 103 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
103     ia4=ia4+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir4)=-s41
        r(ir1+mrowr)=-s12
        r(ir2+mrowr)=-s22
        r(ir3+mrowr)=-s32
        r(ir4+mrowr)=-s42
        r(ir1+2*mrowr)=-s13
        r(ir2+2*mrowr)=-s23
        r(ir3+2*mrowr)=-s33
        r(ir4+2*mrowr)=-s43
        goto 300
c
205     ib1=ibb
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        s13=0
        s23=0
        s33=0
        s43=0
        s14=0
        s24=0
        s34=0
        s44=0
        s15=0
        s25=0
        s35=0
        s45=0
        do 105 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        s34=a(ia3)*b(ib4)+s34
        s44=a(ia4)*b(ib4)+s44
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        s25=a(ia2)*b(ib5)+s25
        s35=a(ia3)*b(ib5)+s35
        s45=a(ia4)*b(ib5)+s45
        ib5=ib5+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
105     ia4=ia4+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir4)=-s41
        r(ir1+mrowr)=-s12
        r(ir2+mrowr)=-s22
        r(ir3+mrowr)=-s32
        r(ir4+mrowr)=-s42
        r(ir1+2*mrowr)=-s13
        r(ir2+2*mrowr)=-s23
        r(ir3+2*mrowr)=-s33
        r(ir4+2*mrowr)=-s43
        r(ir1+3*mrowr)=-s14
        r(ir2+3*mrowr)=-s24
        r(ir3+3*mrowr)=-s34
        r(ir4+3*mrowr)=-s44
        r(ir1+4*mrowr)=-s15
        r(ir2+4*mrowr)=-s25
        r(ir3+4*mrowr)=-s35
        r(ir4+4*mrowr)=-s45
300   iaa=iaa+4*mcola
      if(ncest.eq.0) return
      goto (301,302,303,304,305),ncest
301     ibb=1
        ir1=irr
        do 2001 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        s11=0
        s12=0
        do 1001 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        ib2=ib2+mcolb
1001    ia1=ia1+mrowa
        r(ir1)=-s11
        r(ir1+mrowr)=-s12
        ir1=ir1+2*mrowr
2001    continue
        if(mod(nrow,2).eq.0) return
        s11=0
        do 1011 k=1,nlink
        s11=a(iaa)*b(ibb)+s11
        ibb=ibb+mcolb
1011    iaa=iaa+mrowa
        r(ir1)=-s11
        return
c
302     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        do 2002 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        s11=0
        s21=0
        s12=0
        s22=0
        do 1002 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        ia1=ia1+mrowa
1002    ia2=ia2+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir1+mrowr)=-s12
        r(ir2+mrowr)=-s22
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
2002    continue
        if(mod(nrow,2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        s11=0
        s21=0
        do 1012 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        ia1=ia1+mrowa
        s21=a(ia2)*b(ibb)+s21
        ia2=ia2+mrowa
1012    ibb=ibb+mcolb
        r(ir1)=-s11
        r(ir2)=-s21
        return
c
303     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        ir3=ir2+mcolr
        do 2003 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0
        s21=0
        s31=0
        s12=0
        s22=0
        s32=0
        do 1003 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
1003    ia3=ia3+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir1+mrowr)=-s12
        r(ir2+mrowr)=-s22
        r(ir3+mrowr)=-s32
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
        ir3=ir3+2*mrowr
2003    continue
        if(mod(nrow,2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0
        s21=0
        s31=0
        do 1013 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        s21=a(ia2)*b(ibb)+s21
        s31=a(ia3)*b(ibb)+s31
        ibb=ibb+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
1013    ia3=ia3+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
304     return
305     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        ir3=ir2+mcolr
        ir4=ir3+mcolr
        ir5=ir4+mcolr
        do 2005 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s51=0
        s12=0
        s22=0
        s32=0
        s42=0
        s52=0
        do 1005 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        s52=a(ia5)*b(ib2)+s52
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
        ia4=ia4+mrowa
1005    ia5=ia5+mrowa
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir4)=-s41
        r(ir5)=-s51
        r(ir1+mrowr)=-s12
        r(ir2+mrowr)=-s22
        r(ir3+mrowr)=-s32
        r(ir4+mrowr)=-s42
        r(ir5+mrowr)=-s52
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
        ir3=ir3+2*mrowr
        ir4=ir4+2*mrowr
        ir5=ir5+2*mrowr
2005    continue
        if(mod(nrow,2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s51=0
        do 1015 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        s21=a(ia2)*b(ibb)+s21
        s31=a(ia3)*b(ibb)+s31
        s41=a(ia4)*b(ibb)+s41
        s51=a(ia5)*b(ibb)+s51
        ibb=ibb+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
        ia4=ia4+mrowa
        ia5=ia5+mrowa
1015    continue
        r(ir1)=-s11
        r(ir2)=-s21
        r(ir3)=-s31
        r(ir4)=-s41
        r(ir5)=-s51
        return
      end
c--------------------------------------------------------------------
      subroutine mxmbn(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension r(1),a(1),b(1)
c... r(ncol,nrow)=r(ncol,nrow)-a(ncol,nlink)*b(nlink,nrow) matrix mult
      iaa=1
      irr=1
      nrest=mod(nrow,4)
      if(nrest.eq.1.and.nrow.gt.1) nrest=5
      nre=nrow-nrest
      ncest=mod(ncol,4)
      if(ncest.eq.1.and.ncol.gt.1) ncest=5
      nce=ncol-ncest
      do 300 i=1,nce,4
      ibb=1
      ir1=irr
      ir2=ir1+mcolr
      ir3=ir2+mcolr
      ir4=ir3+mcolr
      irr=ir4+mcolr
      do 200 j=1,nre,4
      ib1=ibb
      ib2=ib1+mrowb
      ib3=ib2+mrowb
      ib4=ib3+mrowb
      ibb=ib4+mrowb
      ia1=iaa
      ia2=ia1+mcola
      ia3=ia2+mcola
      ia4=ia3+mcola
      s11=0
      s21=0
      s31=0
      s41=0
      s12=0
      s22=0
      s32=0
      s42=0
      s13=0
      s23=0
      s33=0
      s43=0
      s14=0
      s24=0
      s34=0
      s44=0
      do 100 k=1,nlink
      s11=a(ia1)*b(ib1)+s11
      s21=a(ia2)*b(ib1)+s21
      s31=a(ia3)*b(ib1)+s31
      s41=a(ia4)*b(ib1)+s41
      ib1=ib1+mcolb
      s12=a(ia1)*b(ib2)+s12
      s22=a(ia2)*b(ib2)+s22
      s32=a(ia3)*b(ib2)+s32
      s42=a(ia4)*b(ib2)+s42
      ib2=ib2+mcolb
      s13=a(ia1)*b(ib3)+s13
      s23=a(ia2)*b(ib3)+s23
      s33=a(ia3)*b(ib3)+s33
      s43=a(ia4)*b(ib3)+s43
      ib3=ib3+mcolb
      s14=a(ia1)*b(ib4)+s14
      s24=a(ia2)*b(ib4)+s24
      s34=a(ia3)*b(ib4)+s34
      s44=a(ia4)*b(ib4)+s44
      ib4=ib4+mcolb
      ia1=ia1+mrowa
      ia2=ia2+mrowa
      ia3=ia3+mrowa
100   ia4=ia4+mrowa
      r(ir1)=r(ir1)-s11
      r(ir2)=r(ir2)-s21
      r(ir3)=r(ir3)-s31
      r(ir4)=r(ir4)-s41
      r(ir1+mrowr)=r(ir1+mrowr)-s12
      r(ir2+mrowr)=r(ir2+mrowr)-s22
      r(ir3+mrowr)=r(ir3+mrowr)-s32
      r(ir4+mrowr)=r(ir4+mrowr)-s42
      r(ir1+2*mrowr)=r(ir1+2*mrowr)-s13
      r(ir2+2*mrowr)=r(ir2+2*mrowr)-s23
      r(ir3+2*mrowr)=r(ir3+2*mrowr)-s33
      r(ir4+2*mrowr)=r(ir4+2*mrowr)-s43
      r(ir1+3*mrowr)=r(ir1+3*mrowr)-s14
      r(ir2+3*mrowr)=r(ir2+3*mrowr)-s24
      r(ir3+3*mrowr)=r(ir3+3*mrowr)-s34
      r(ir4+3*mrowr)=r(ir4+3*mrowr)-s44
      ir1=ir1+4*mrowr
      ir2=ir2+4*mrowr
      ir3=ir3+4*mrowr
      ir4=ir4+4*mrowr
200   continue
      if(nrest.eq.0) goto 300
      goto (201,202,203,300,205),nrest
201     ib1=ibb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        do 101 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
101     ia4=ia4+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir4)=r(ir4)-s41
        goto 300
c
202     ib1=ibb
        ib2=ib1+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        do 102 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
102     ia4=ia4+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir4)=r(ir4)-s41
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        r(ir2+mrowr)=r(ir2+mrowr)-s22
        r(ir3+mrowr)=r(ir3+mrowr)-s32
        r(ir4+mrowr)=r(ir4+mrowr)-s42
        goto 300
c
203     ib1=ibb
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        s13=0
        s23=0
        s33=0
        s43=0
        do 103 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
103     ia4=ia4+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir4)=r(ir4)-s41
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        r(ir2+mrowr)=r(ir2+mrowr)-s22
        r(ir3+mrowr)=r(ir3+mrowr)-s32
        r(ir4+mrowr)=r(ir4+mrowr)-s42
        r(ir1+2*mrowr)=r(ir1+2*mrowr)-s13
        r(ir2+2*mrowr)=r(ir2+2*mrowr)-s23
        r(ir3+2*mrowr)=r(ir3+2*mrowr)-s33
        r(ir4+2*mrowr)=r(ir4+2*mrowr)-s43
        goto 300
c
205     ib1=ibb
        ib2=ib1+mrowb
        ib3=ib2+mrowb
        ib4=ib3+mrowb
        ib5=ib4+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s12=0
        s22=0
        s32=0
        s42=0
        s13=0
        s23=0
        s33=0
        s43=0
        s14=0
        s24=0
        s34=0
        s44=0
        s15=0
        s25=0
        s35=0
        s45=0
        do 105 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        ib2=ib2+mcolb
        s13=a(ia1)*b(ib3)+s13
        s23=a(ia2)*b(ib3)+s23
        s33=a(ia3)*b(ib3)+s33
        s43=a(ia4)*b(ib3)+s43
        ib3=ib3+mcolb
        s14=a(ia1)*b(ib4)+s14
        s24=a(ia2)*b(ib4)+s24
        s34=a(ia3)*b(ib4)+s34
        s44=a(ia4)*b(ib4)+s44
        ib4=ib4+mcolb
        s15=a(ia1)*b(ib5)+s15
        s25=a(ia2)*b(ib5)+s25
        s35=a(ia3)*b(ib5)+s35
        s45=a(ia4)*b(ib5)+s45
        ib5=ib5+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
105     ia4=ia4+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir4)=r(ir4)-s41
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        r(ir2+mrowr)=r(ir2+mrowr)-s22
        r(ir3+mrowr)=r(ir3+mrowr)-s32
        r(ir4+mrowr)=r(ir4+mrowr)-s42
        r(ir1+2*mrowr)=r(ir1+2*mrowr)-s13
        r(ir2+2*mrowr)=r(ir2+2*mrowr)-s23
        r(ir3+2*mrowr)=r(ir3+2*mrowr)-s33
        r(ir4+2*mrowr)=r(ir4+2*mrowr)-s43
        r(ir1+3*mrowr)=r(ir1+3*mrowr)-s14
        r(ir2+3*mrowr)=r(ir2+3*mrowr)-s24
        r(ir3+3*mrowr)=r(ir3+3*mrowr)-s34
        r(ir4+3*mrowr)=r(ir4+3*mrowr)-s44
        r(ir1+4*mrowr)=r(ir1+4*mrowr)-s15
        r(ir2+4*mrowr)=r(ir2+4*mrowr)-s25
        r(ir3+4*mrowr)=r(ir3+4*mrowr)-s35
        r(ir4+4*mrowr)=r(ir4+4*mrowr)-s45
300   iaa=iaa+4*mcola
      if(ncest.eq.0) return
      goto (301,302,303,304,305),ncest
301     ibb=1
        ir1=irr
        do 2001 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        s11=0
        s12=0
        do 1001 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        ib2=ib2+mcolb
1001    ia1=ia1+mrowa
        r(ir1)=r(ir1)-s11
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        ir1=ir1+2*mrowr
2001    continue
        if(mod(nrow,2).eq.0) return
        s11=0
        do 1011 k=1,nlink
        s11=a(iaa)*b(ibb)+s11
        ibb=ibb+mcolb
1011    iaa=iaa+mrowa
        r(ir1)=r(ir1)-s11
        return
c
302     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        do 2002 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        s11=0
        s21=0
        s12=0
        s22=0
        do 1002 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        ib2=ib2+mcolb
        ia1=ia1+mrowa
1002    ia2=ia2+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        r(ir2+mrowr)=r(ir2+mrowr)-s22
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
2002    continue
        if(mod(nrow,2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        s11=0
        s21=0
        do 1012 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        ia1=ia1+mrowa
        s21=a(ia2)*b(ibb)+s21
        ia2=ia2+mrowa
1012    ibb=ibb+mcolb
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        return
c
303     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        ir3=ir2+mcolr
        do 2003 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0
        s21=0
        s31=0
        s12=0
        s22=0
        s32=0
        do 1003 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
1003    ia3=ia3+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        r(ir2+mrowr)=r(ir2+mrowr)-s22
        r(ir3+mrowr)=r(ir3+mrowr)-s32
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
        ir3=ir3+2*mrowr
2003    continue
        if(mod(nrow,2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        s11=0
        s21=0
        s31=0
        do 1013 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        s21=a(ia2)*b(ibb)+s21
        s31=a(ia3)*b(ibb)+s31
        ibb=ibb+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
1013    ia3=ia3+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
304     return
305     ibb=1
        ir1=irr
        ir2=ir1+mcolr
        ir3=ir2+mcolr
        ir4=ir3+mcolr
        ir5=ir4+mcolr
        do 2005 j=1,nrow-1,2
        ib1=ibb
        ib2=ib1+mrowb
        ibb=ib2+mrowb
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s51=0
        s12=0
        s22=0
        s32=0
        s42=0
        s52=0
        do 1005 k=1,nlink
        s11=a(ia1)*b(ib1)+s11
        s21=a(ia2)*b(ib1)+s21
        s31=a(ia3)*b(ib1)+s31
        s41=a(ia4)*b(ib1)+s41
        s51=a(ia5)*b(ib1)+s51
        ib1=ib1+mcolb
        s12=a(ia1)*b(ib2)+s12
        s22=a(ia2)*b(ib2)+s22
        s32=a(ia3)*b(ib2)+s32
        s42=a(ia4)*b(ib2)+s42
        s52=a(ia5)*b(ib2)+s52
        ib2=ib2+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
        ia4=ia4+mrowa
1005    ia5=ia5+mrowa
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir4)=r(ir4)-s41
        r(ir5)=r(ir5)-s51
        r(ir1+mrowr)=r(ir1+mrowr)-s12
        r(ir2+mrowr)=r(ir2+mrowr)-s22
        r(ir3+mrowr)=r(ir3+mrowr)-s32
        r(ir4+mrowr)=r(ir4+mrowr)-s42
        r(ir5+mrowr)=r(ir5+mrowr)-s52
        ir1=ir1+2*mrowr
        ir2=ir2+2*mrowr
        ir3=ir3+2*mrowr
        ir4=ir4+2*mrowr
        ir5=ir5+2*mrowr
2005    continue
        if(mod(nrow,2).eq.0) return
        ia1=iaa
        ia2=ia1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        ia5=ia4+mcola
        s11=0
        s21=0
        s31=0
        s41=0
        s51=0
        do 1015 k=1,nlink
        s11=a(ia1)*b(ibb)+s11
        s21=a(ia2)*b(ibb)+s21
        s31=a(ia3)*b(ibb)+s31
        s41=a(ia4)*b(ibb)+s41
        s51=a(ia5)*b(ibb)+s51
        ibb=ibb+mcolb
        ia1=ia1+mrowa
        ia2=ia2+mrowa
        ia3=ia3+mrowa
        ia4=ia4+mrowa
        ia5=ia5+mrowa
1015    continue
        r(ir1)=r(ir1)-s11
        r(ir2)=r(ir2)-s21
        r(ir3)=r(ir3)-s31
        r(ir4)=r(ir4)-s41
        r(ir5)=r(ir5)-s51
        return
      end
c--------------------------------------------------------------------
      subroutine mxva(a,mcola,mrowa,b,mcolb,
     1                r,mcolr,ncol,nlink)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension r(1),a(1),b(1)
c...  r(ncol)=a(ncol,nlink)*b(nlink) matrix mult
* mcola is spacing between adjacent column elements of a (column stride)
* mrowa is spacing between adjacent row elements of a (row stride)
* mcolb is spacing between adjacent elements of b
* mcolr is spacing between elements of product vector
* ncol is number of rows in a and number of elements in r
* nlink is number of columns in a and number of elements in b
      if(ncol.gt.4) goto 50
        ia1=1
        ib1=1
        if(ncol.eq.1) then
        s1=0
        do 10 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
10      ib1=ib1+mcolb
        r(1)=s1
        return
      else if(ncol.eq.2) then
        ia2=1+mcola
        s1=0
        s2=0
        do 20 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
20      ib1=ib1+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        return
      else if(ncol.eq.3) then
        ia2=1+mcola
        ia3=ia2+mcola
        s1=0
        s2=0
        s3=0
        do 30 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
30      ib1=ib1+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        r(1+2*mcolr)=s3
        return
      else
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        do 40 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
        s4=a(ia4)*b(ib1)+s4
        ia4=ia4+mrowa
40      ib1=ib1+mcolb
        r(1)=s1
        r(1+mcolr)=s2
        r(1+2*mcolr)=s3
        r(1+3*mcolr)=s4
        return
      end if
c
50    mcola4=4*mcola
      mcolr4=4*mcolr
      iaa1=1
      iaa2=1+mcola
      iaa3=iaa2+mcola
      iaa4=iaa3+mcola
      ir1=1
      ir2=1+mcolr
      ir3=ir2+mcolr
      ir4=ir3+mcolr
      ncoll=(ncol/4)*4
      do 200 i=1,ncoll,4
      ia1=iaa1
      ia2=iaa2
      ia3=iaa3
      ia4=iaa4
      ib1=1
      s1=0
      s2=0
      s3=0
      s4=0
      do 100 k=1,nlink
      s1=a(ia1)*b(ib1)+s1
      ia1=ia1+mrowa
      s2=a(ia2)*b(ib1)+s2
      ia2=ia2+mrowa
      s3=a(ia3)*b(ib1)+s3
      ia3=ia3+mrowa
      s4=a(ia4)*b(ib1)+s4
      ia4=ia4+mrowa
100   ib1=ib1+mcolb
110   r(ir1)=s1
      ir1=ir1+mcolr4
      r(ir2)=s2
      ir2=ir2+mcolr4
      r(ir3)=s3
      ir3=ir3+mcolr4
      r(ir4)=s4
      ir4=ir4+mcolr4
      iaa1=iaa1+mcola4
      iaa2=iaa2+mcola4
      iaa3=iaa3+mcola4
200   iaa4=iaa4+mcola4
      ncoll=ncol-ncoll
      if(ncoll.eq.0) then
        return
      else if(ncoll.eq.1) then
        ia1=iaa1
        ib1=1
        s1=0
        do 101 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
101     ib1=ib1+mcolb
        r(ir1)=s1
      else if(ncoll.eq.2) then
        ia1=iaa1
        ia2=iaa2
        ib1=1
        s1=0
        s2=0
        do 102 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
102     ib1=ib1+mcolb
        r(ir1)=s1
        r(ir2)=s2
      else if(ncoll.eq.3) then
        ia1=iaa1
        ia2=iaa2
        ia3=iaa3
        ib1=1
        s1=0
        s2=0
        s3=0
        do 103 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
103     ib1=ib1+mcolb
        r(ir1)=s1
        r(ir2)=s2
        r(ir3)=s3
      end if
      return
      end
c--------------------------------------------------------------------
      subroutine mxvb(a,mcola,mrowa,b,mcolb,
     1                r,mcolr,ncol,nlink)
c--------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension r(1),a(1),b(1)
c...  r(ncol)=a(ncol,nlink)*b(nlink) matrix mult
      if(ncol.gt.4) goto 50
        ia1=1
        ib1=1
        if(ncol.eq.1) then
        s1=0
        do 10 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
10      ib1=ib1+mcolb
        r(1)=r(1)+s1
        return
      else if(ncol.eq.2) then
        ia2=1+mcola
        s1=0
        s2=0
        do 20 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
20      ib1=ib1+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        return
      else if(ncol.eq.3) then
        ia2=1+mcola
        ia3=ia2+mcola
        s1=0
        s2=0
        s3=0
        do 30 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
30      ib1=ib1+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        r(1+2*mcolr)=r(1+2*mcolr)+s3
        return
      else
        ia2=1+mcola
        ia3=ia2+mcola
        ia4=ia3+mcola
        s1=0
        s2=0
        s3=0
        s4=0
        do 40 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
        s4=a(ia4)*b(ib1)+s4
        ia4=ia4+mrowa
40      ib1=ib1+mcolb
        r(1)=r(1)+s1
        r(1+mcolr)=r(1+mcolr)+s2
        r(1+2*mcolr)=r(1+2*mcolr)+s3
        r(1+3*mcolr)=r(1+3*mcolr)+s4
        return
      end if
c
50    mcola4=4*mcola
      mcolr4=4*mcolr
      iaa1=1
      iaa2=1+mcola
      iaa3=iaa2+mcola
      iaa4=iaa3+mcola
      ir1=1
      ir2=1+mcolr
      ir3=ir2+mcolr
      ir4=ir3+mcolr
      ncoll=(ncol/4)*4
      do 200 i=1,ncoll,4
      ia1=iaa1
      ia2=iaa2
      ia3=iaa3
      ia4=iaa4
      ib1=1
      s1=0
      s2=0
      s3=0
      s4=0
      do 100 k=1,nlink
      s1=a(ia1)*b(ib1)+s1
      ia1=ia1+mrowa
      s2=a(ia2)*b(ib1)+s2
      ia2=ia2+mrowa
      s3=a(ia3)*b(ib1)+s3
      ia3=ia3+mrowa
      s4=a(ia4)*b(ib1)+s4
      ia4=ia4+mrowa
100   ib1=ib1+mcolb
110   r(ir1)=r(ir1)+s1
      ir1=ir1+mcolr4
      r(ir2)=r(ir2)+s2
      ir2=ir2+mcolr4
      r(ir3)=r(ir3)+s3
      ir3=ir3+mcolr4
      r(ir4)=r(ir4)+s4
      ir4=ir4+mcolr4
      iaa1=iaa1+mcola4
      iaa2=iaa2+mcola4
      iaa3=iaa3+mcola4
200   iaa4=iaa4+mcola4
      ncoll=ncol-ncoll
      if(ncoll.eq.0) then
        return
      else if(ncoll.eq.1) then
        ia1=iaa1
        ib1=1
        s1=0
        do 101 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
101     ib1=ib1+mcolb
        r(ir1)=r(ir1)+s1
      else if(ncoll.eq.2) then
        ia1=iaa1
        ia2=iaa2
        ib1=1
        s1=0
        s2=0
        do 102 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
102     ib1=ib1+mcolb
        r(ir1)=r(ir1)+s1
        r(ir2)=r(ir2)+s2
      else if(ncoll.eq.3) then
        ia1=iaa1
        ia2=iaa2
        ia3=iaa3
        ib1=1
        s1=0
        s2=0
        s3=0
        do 103 k=1,nlink
        s1=a(ia1)*b(ib1)+s1
        ia1=ia1+mrowa
        s2=a(ia2)*b(ib1)+s2
        ia2=ia2+mrowa
        s3=a(ia3)*b(ib1)+s3
        ia3=ia3+mrowa
103     ib1=ib1+mcolb
        r(ir1)=r(ir1)+s1
        r(ir2)=r(ir2)+s2
        r(ir3)=r(ir3)+s3
      end if
      return
      end
#endif
c -------------------------
#if (defined(HIB_UNIX_AIX) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)) && !defined(HIB_UNIX_HP800)
c;      subroutine dsymv ( uplo, n, alpha, a, lda, x, incx,
c;     $                   beta, y, incy )
c;*     .. Scalar Arguments ..
c;      double precision   alpha, beta
c;      integer            incx, incy, lda, n
c;      character*1        uplo
c;*     .. Array Arguments ..
c;      double precision   a( lda, * ), x( * ), y( * )
c;*     ..
c;*
c;*  Purpose
c;*  =======
c;*
c;*  DSYMV  performs the matrix-vector  operation
c;*
c;*     y := alpha*A*x + beta*y,
c;*
c;*  where alpha and beta are scalars, x and y are n element vectors and
c;*  A is an n by n symmetric matrix.
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
c;*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
c;*           Before entry with  UPLO = 'U' or 'u', the leading n by n
c;*           upper triangular part of the array A must contain the upper
c;*           triangular part of the symmetric matrix and the strictly
c;*           lower triangular part of A is not referenced.
c;*           Before entry with UPLO = 'L' or 'l', the leading n by n
c;*           lower triangular part of the array A must contain the lower
c;*           triangular part of the symmetric matrix and the strictly
c;*           upper triangular part of A is not referenced.
c;*           Unchanged on exit.
c;*
c;*  LDA    - INTEGER.
c;*           On entry, LDA specifies the first dimension of A as declared
c;*           in the calling (sub) program. LDA must be at least
c;*           max( 1, n ).
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
c;*  BETA   - DOUBLE PRECISION.
c;*           On entry, BETA specifies the scalar beta. When BETA is
c;*           supplied as zero then Y need not be set on input.
c;*           Unchanged on exit.
c;*
c;*  Y      - DOUBLE PRECISION array of dimension at least
c;*           ( 1 + ( n - 1 )*abs( INCY ) ).
c;*           Before entry, the incremented array Y must contain the n
c;*           element vector y. On exit, Y is overwritten by the updated
c;*           vector y.
c;*
c;*  INCY   - INTEGER.
c;*           On entry, INCY specifies the increment for the elements of
c;*           Y. INCY must not be zero.
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
c;      double precision   one         , zero
c;      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c;*     .. Local Scalars ..
c;      double precision   temp1, temp2
c;      integer            i, info, ix, iy, j, jx, jy, kx, ky
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
c;      else if( lda.lt.max( 1, n ) )then
c;         info = 5
c;      else if( incx.eq.0 )then
c;         info = 7
c;      else if( incy.eq.0 )then
c;         info = 10
c;      end if
c;      if( info.ne.0 )then
c;         call xerbla( 'DSYMV ', info )
c;         return
c;      end if
c;*
c;*     Quick return if possible.
c;*
c;      if( ( n.eq.0 ).or.( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
c;     $   return
c;*
c;*     Set up the start points in  X  and  Y.
c;*
c;      if( incx.gt.0 )then
c;         kx = 1
c;      else
c;         kx = 1 - ( n - 1 )*incx
c;      end if
c;      if( incy.gt.0 )then
c;         ky = 1
c;      else
c;         ky = 1 - ( n - 1 )*incy
c;      end if
c;*
c;*     Start the operations. In this version the elements of A are
c;*     accessed sequentially with one pass through the triangular part
c;*     of A.
c;*
c;*     First form  y := beta*y.
c;*
c;      if( beta.ne.one )then
c;         if( incy.eq.1 )then
c;            if( beta.eq.zero )then
c;               do 10, i = 1, n
c;                  y( i ) = zero
c;   10          continue
c;            else
c;               do 20, i = 1, n
c;                  y( i ) = beta*y( i )
c;   20          continue
c;            end if
c;         else
c;            iy = ky
c;            if( beta.eq.zero )then
c;               do 30, i = 1, n
c;                  y( iy ) = zero
c;                  iy      = iy   + incy
c;   30          continue
c;            else
c;               do 40, i = 1, n
c;                  y( iy ) = beta*y( iy )
c;                  iy      = iy           + incy
c;   40          continue
c;            end if
c;         end if
c;      end if
c;      if( alpha.eq.zero )
c;     $   return
c;      if( lsame( uplo, 'U' ) )then
c;*
c;*        Form  y  when A is stored in upper triangle.
c;*
c;         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
c;            do 60, j = 1, n
c;               temp1 = alpha*x( j )
c;               temp2 = zero
c;               do 50, i = 1, j - 1
c;                  y( i ) = y( i ) + temp1*a( i, j )
c;                  temp2  = temp2  + a( i, j )*x( i )
c;   50          continue
c;               y( j ) = y( j ) + temp1*a( j, j ) + alpha*temp2
c;   60       continue
c;         else
c;            jx = kx
c;            jy = ky
c;            do 80, j = 1, n
c;               temp1 = alpha*x( jx )
c;               temp2 = zero
c;               ix    = kx
c;               iy    = ky
c;               do 70, i = 1, j - 1
c;                  y( iy ) = y( iy ) + temp1*a( i, j )
c;                  temp2   = temp2   + a( i, j )*x( ix )
c;                  ix      = ix      + incx
c;                  iy      = iy      + incy
c;   70          continue
c;               y( jy ) = y( jy ) + temp1*a( j, j ) + alpha*temp2
c;               jx      = jx      + incx
c;               jy      = jy      + incy
c;   80       continue
c;         end if
c;      else
c;*
c;*        Form  y  when A is stored in lower triangle.
c;*
c;         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
c;            do 100, j = 1, n
c;               temp1  = alpha*x( j )
c;               temp2  = zero
c;               y( j ) = y( j )       + temp1*a( j, j )
c;               do 90, i = j + 1, n
c;                  y( i ) = y( i ) + temp1*a( i, j )
c;                  temp2  = temp2  + a( i, j )*x( i )
c;   90          continue
c;               y( j ) = y( j ) + alpha*temp2
c;  100       continue
c;         else
c;            jx = kx
c;            jy = ky
c;            do 120, j = 1, n
c;               temp1   = alpha*x( jx )
c;               temp2   = zero
c;               y( jy ) = y( jy )       + temp1*a( j, j )
c;               ix      = jx
c;               iy      = jy
c;               do 110, i = j + 1, n
c;                  ix      = ix      + incx
c;                  iy      = iy      + incy
c;                  y( iy ) = y( iy ) + temp1*a( i, j )
c;                  temp2   = temp2   + a( i, j )*x( ix )
c;  110          continue
c;               y( jy ) = y( jy ) + alpha*temp2
c;               jx      = jx      + incx
c;               jy      = jy      + incy
c;  120       continue
c;         end if
c;      end if
c;*
c;      return
c;*
c;*     End of DSYMV .
c;*
c;      end
c;      subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx,
c;     $                   beta, y, incy )
c;*     .. Scalar Arguments ..
c;      double precision   alpha, beta
c;      integer            incx, incy, lda, m, n
c;      character*1        trans
c;*     .. Array Arguments ..
c;      double precision   a( lda, * ), x( * ), y( * )
c;*     ..
c;*
c;*  Purpose
c;*  =======
c;*
c;*  DGEMV  performs one of the matrix-vector operations
c;*
c;*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
c;*
c;*  where alpha and beta are scalars, x and y are vectors and A is an
c;*  m by n matrix.
c;*
c;*  Parameters
c;*  ==========
c;*
c;*  TRANS  - CHARACTER*1.
c;*           On entry, TRANS specifies the operation to be performed as
c;*           follows:
c;*
c;*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
c;*
c;*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
c;*
c;*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
c;*
c;*           Unchanged on exit.
c;*
c;*  M      - INTEGER.
c;*           On entry, M specifies the number of rows of the matrix A.
c;*           M must be at least zero.
c;*           Unchanged on exit.
c;*
c;*  N      - INTEGER.
c;*           On entry, N specifies the number of columns of the matrix A.
c;*           N must be at least zero.
c;*           Unchanged on exit.
c;*
c;*  ALPHA  - DOUBLE PRECISION.
c;*           On entry, ALPHA specifies the scalar alpha.
c;*           Unchanged on exit.
c;*
c;*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
c;*           Before entry, the leading m by n part of the array A must
c;*           contain the matrix of coefficients.
c;*           Unchanged on exit.
c;*
c;*  LDA    - INTEGER.
c;*           On entry, LDA specifies the first dimension of A as declared
c;*           in the calling (sub) program. LDA must be at least
c;*           max( 1, m ).
c;*           Unchanged on exit.
c;*
c;*  X      - DOUBLE PRECISION array of DIMENSION at least
c;*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
c;*           and at least
c;*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
c;*           Before entry, the incremented array X must contain the
c;*           vector x.
c;*           Unchanged on exit.
c;*
c;*  INCX   - INTEGER.
c;*           On entry, INCX specifies the increment for the elements of
c;*           X. INCX must not be zero.
c;*           Unchanged on exit.
c;*
c;*  BETA   - DOUBLE PRECISION.
c;*           On entry, BETA specifies the scalar beta. When BETA is
c;*           supplied as zero then Y need not be set on input.
c;*           Unchanged on exit.
c;*
c;*  Y      - DOUBLE PRECISION array of DIMENSION at least
c;*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
c;*           and at least
c;*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
c;*           Before entry with BETA non-zero, the incremented array Y
c;*           must contain the vector y. On exit, Y is overwritten by the
c;*           updated vector y.
c;*
c;*  INCY   - INTEGER.
c;*           On entry, INCY specifies the increment for the elements of
c;*           Y. INCY must not be zero.
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
c;      double precision   one         , zero
c;      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c;*     .. Local Scalars ..
c;      double precision   temp
c;      integer            i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
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
c;      if     ( .not.lsame( trans, 'N' ).and.
c;     $         .not.lsame( trans, 'T' ).and.
c;     $         .not.lsame( trans, 'C' )      )then
c;         info = 1
c;      else if( m.lt.0 )then
c;         info = 2
c;      else if( n.lt.0 )then
c;         info = 3
c;      else if( lda.lt.max( 1, m ) )then
c;         info = 6
c;      else if( incx.eq.0 )then
c;         info = 8
c;      else if( incy.eq.0 )then
c;         info = 11
c;      end if
c;      if( info.ne.0 )then
c;         call xerbla( 'DGEMV ', info )
c;         return
c;      end if
c;*
c;*     Quick return if possible.
c;*
c;      if( ( m.eq.0 ).or.( n.eq.0 ).or.
c;     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
c;     $   return
c;*
c;*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
c;*     up the start points in  X  and  Y.
c;*
c;      if( lsame( trans, 'N' ) )then
c;         lenx = n
c;         leny = m
c;      else
c;         lenx = m
c;         leny = n
c;      end if
c;      if( incx.gt.0 )then
c;         kx = 1
c;      else
c;         kx = 1 - ( lenx - 1 )*incx
c;      end if
c;      if( incy.gt.0 )then
c;         ky = 1
c;      else
c;         ky = 1 - ( leny - 1 )*incy
c;      end if
c;*
c;*     Start the operations. In this version the elements of A are
c;*     accessed sequentially with one pass through A.
c;*
c;*     First form  y := beta*y.
c;*
c;      if( beta.ne.one )then
c;         if( incy.eq.1 )then
c;            if( beta.eq.zero )then
c;               do 10, i = 1, leny
c;                  y( i ) = zero
c;   10          continue
c;            else
c;               do 20, i = 1, leny
c;                  y( i ) = beta*y( i )
c;   20          continue
c;            end if
c;         else
c;            iy = ky
c;            if( beta.eq.zero )then
c;               do 30, i = 1, leny
c;                  y( iy ) = zero
c;                  iy      = iy   + incy
c;   30          continue
c;            else
c;               do 40, i = 1, leny
c;                  y( iy ) = beta*y( iy )
c;                  iy      = iy           + incy
c;   40          continue
c;            end if
c;         end if
c;      end if
c;      if( alpha.eq.zero )
c;     $   return
c;      if( lsame( trans, 'N' ) )then
c;*
c;*        Form  y := alpha*A*x + y.
c;*
c;         jx = kx
c;         if( incy.eq.1 )then
c;            do 60, j = 1, n
c;               if( x( jx ).ne.zero )then
c;                  temp = alpha*x( jx )
c;                  do 50, i = 1, m
c;                     y( i ) = y( i ) + temp*a( i, j )
c;   50             continue
c;               end if
c;               jx = jx + incx
c;   60       continue
c;         else
c;            do 80, j = 1, n
c;               if( x( jx ).ne.zero )then
c;                  temp = alpha*x( jx )
c;                  iy   = ky
c;                  do 70, i = 1, m
c;                     y( iy ) = y( iy ) + temp*a( i, j )
c;                     iy      = iy      + incy
c;   70             continue
c;               end if
c;               jx = jx + incx
c;   80       continue
c;         end if
c;      else
c;*
c;*        Form  y := alpha*A'*x + y.
c;*
c;         jy = ky
c;         if( incx.eq.1 )then
c;            do 100, j = 1, n
c;               temp = zero
c;               do 90, i = 1, m
c;                  temp = temp + a( i, j )*x( i )
c;   90          continue
c;               y( jy ) = y( jy ) + alpha*temp
c;               jy      = jy      + incy
c;  100       continue
c;         else
c;            do 120, j = 1, n
c;               temp = zero
c;               ix   = kx
c;               do 110, i = 1, m
c;                  temp = temp + a( i, j )*x( ix )
c;                  ix   = ix   + incx
c;  110          continue
c;               y( jy ) = y( jy ) + alpha*temp
c;               jy      = jy      + incy
c;  120       continue
c;         end if
c;      end if
c;*
c;      return
c;*
c;*     End of DGEMV .
c;*
c;      end
c;      subroutine dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb,
c;     $                   beta, c, ldc )
c;*     .. Scalar Arguments ..
c;      character*1        transa, transb
c;      integer            m, n, k, lda, ldb, ldc
c;      double precision   alpha, beta
c;*     .. Array Arguments ..
c;      double precision   a( lda, * ), b( ldb, * ), c( ldc, * )
c;*     ..
c;*
c;*  Purpose
c;*  =======
c;*
c;*  DGEMM  performs one of the matrix-matrix operations
c;*
c;*     C := alpha*op( A )*op( B ) + beta*C,
c;*
c;*  where  op( X ) is one of
c;*
c;*     op( X ) = X   or   op( X ) = X',
c;*
c;*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
c;*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
c;*
c;*  Parameters
c;*  ==========
c;*
c;*  TRANSA - CHARACTER*1.
c;*           On entry, TRANSA specifies the form of op( A ) to be used in
c;*           the matrix multiplication as follows:
c;*
c;*              TRANSA = 'N' or 'n',  op( A ) = A.
c;*
c;*              TRANSA = 'T' or 't',  op( A ) = A'.
c;*
c;*              TRANSA = 'C' or 'c',  op( A ) = A'.
c;*
c;*           Unchanged on exit.
c;*
c;*  TRANSB - CHARACTER*1.
c;*           On entry, TRANSB specifies the form of op( B ) to be used in
c;*           the matrix multiplication as follows:
c;*
c;*              TRANSB = 'N' or 'n',  op( B ) = B.
c;*
c;*              TRANSB = 'T' or 't',  op( B ) = B'.
c;*
c;*              TRANSB = 'C' or 'c',  op( B ) = B'.
c;*
c;*           Unchanged on exit.
c;*
c;*  M      - INTEGER.
c;*           On entry,  M  specifies  the number  of rows  of the  matrix
c;*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
c;*           Unchanged on exit.
c;*
c;*  N      - INTEGER.
c;*           On entry,  N  specifies the number  of columns of the matrix
c;*           op( B ) and the number of columns of the matrix C. N must be
c;*           at least zero.
c;*           Unchanged on exit.
c;*
c;*  K      - INTEGER.
c;*           On entry,  K  specifies  the number of columns of the matrix
c;*           op( A ) and the number of rows of the matrix op( B ). K must
c;*           be at least  zero.
c;*           Unchanged on exit.
c;*
c;*  ALPHA  - DOUBLE PRECISION.
c;*           On entry, ALPHA specifies the scalar alpha.
c;*           Unchanged on exit.
c;*
c;*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
c;*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
c;*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
c;*           part of the array  A  must contain the matrix  A,  otherwise
c;*           the leading  k by m  part of the array  A  must contain  the
c;*           matrix A.
c;*           Unchanged on exit.
c;*
c;*  LDA    - INTEGER.
c;*           On entry, LDA specifies the first dimension of A as declared
c;*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
c;*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c;*           least  max( 1, k ).
c;*           Unchanged on exit.
c;*
c;*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
c;*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
c;*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
c;*           part of the array  B  must contain the matrix  B,  otherwise
c;*           the leading  n by k  part of the array  B  must contain  the
c;*           matrix B.
c;*           Unchanged on exit.
c;*
c;*  LDB    - INTEGER.
c;*           On entry, LDB specifies the first dimension of B as declared
c;*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
c;*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
c;*           least  max( 1, n ).
c;*           Unchanged on exit.
c;*
c;*  BETA   - DOUBLE PRECISION.
c;*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c;*           supplied as zero then C need not be set on input.
c;*           Unchanged on exit.
c;*
c;*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
c;*           Before entry, the leading  m by n  part of the array  C must
c;*           contain the matrix  C,  except when  beta  is zero, in which
c;*           case C need not be set on entry.
c;*           On exit, the array  C  is overwritten by the  m by n  matrix
c;*           ( alpha*op( A )*op( B ) + beta*C ).
c;*
c;*  LDC    - INTEGER.
c;*           On entry, LDC specifies the first dimension of C as declared
c;*           in  the  calling  (sub)  program.   LDC  must  be  at  least
c;*           max( 1, m ).
c;*           Unchanged on exit.
c;*
c;*
c;*  Level 3 Blas routine.
c;*
c;*  -- Written on 8-February-1989.
c;*     Jack Dongarra, Argonne National Laboratory.
c;*     Iain Duff, AERE Harwell.
c;*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c;*     Sven Hammarling, Numerical Algorithms Group Ltd.
c;*
c;*
c;*     .. External Functions ..
c;      logical            lsame
c;      external           lsame
c;*     .. External Subroutines ..
c;      external           xerbla
c;*     .. Intrinsic Functions ..
c;      intrinsic          max
c;*     .. Local Scalars ..
c;      logical            nota, notb
c;      integer            i, info, j, l, ncola, nrowa, nrowb
c;      double precision   temp
c;*     .. Parameters ..
c;      double precision   one         , zero
c;      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c;*     ..
c;*     .. Executable Statements ..
c;*
c;*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
c;*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
c;*     and  columns of  A  and the  number of  rows  of  B  respectively.
c;*
c;      nota  = lsame( transa, 'N' )
c;      notb  = lsame( transb, 'N' )
c;      if( nota )then
c;         nrowa = m
c;         ncola = k
c;      else
c;         nrowa = k
c;         ncola = m
c;      end if
c;      if( notb )then
c;         nrowb = k
c;      else
c;         nrowb = n
c;      end if
c;*
c;*     Test the input parameters.
c;*
c;      info = 0
c;      if(      ( .not.nota                 ).and.
c;     $         ( .not.lsame( transa, 'C' ) ).and.
c;     $         ( .not.lsame( transa, 'T' ) )      )then
c;         info = 1
c;      else if( ( .not.notb                 ).and.
c;     $         ( .not.lsame( transb, 'C' ) ).and.
c;     $         ( .not.lsame( transb, 'T' ) )      )then
c;         info = 2
c;      else if( m  .lt.0               )then
c;         info = 3
c;      else if( n  .lt.0               )then
c;         info = 4
c;      else if( k  .lt.0               )then
c;         info = 5
c;      else if( lda.lt.max( 1, nrowa ) )then
c;         info = 8
c;      else if( ldb.lt.max( 1, nrowb ) )then
c;         info = 10
c;      else if( ldc.lt.max( 1, m     ) )then
c;         info = 13
c;      end if
c;      if( info.ne.0 )then
c;         call xerbla( 'DGEMM ', info )
c;         return
c;      end if
c;*
c;*     Quick return if possible.
c;*
c;      if( ( m.eq.0 ).or.( n.eq.0 ).or.
c;     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
c;     $   return
c;*
c;*     And if  alpha.eq.zero.
c;*
c;      if( alpha.eq.zero )then
c;         if( beta.eq.zero )then
c;            do 20, j = 1, n
c;               do 10, i = 1, m
c;                  c( i, j ) = zero
c;   10          continue
c;   20       continue
c;         else
c;            do 40, j = 1, n
c;               do 30, i = 1, m
c;                  c( i, j ) = beta*c( i, j )
c;   30          continue
c;   40       continue
c;         end if
c;         return
c;      end if
c;*
c;*     Start the operations.
c;*
c;      if( notb )then
c;         if( nota )then
c;*
c;*           Form  C := alpha*A*B + beta*C.
c;*
c;            do 90, j = 1, n
c;               if( beta.eq.zero )then
c;                  do 50, i = 1, m
c;                     c( i, j ) = zero
c;   50             continue
c;               else if( beta.ne.one )then
c;                  do 60, i = 1, m
c;                     c( i, j ) = beta*c( i, j )
c;   60             continue
c;               end if
c;               do 80, l = 1, k
c;                  if( b( l, j ).ne.zero )then
c;                     temp = alpha*b( l, j )
c;                     do 70, i = 1, m
c;                        c( i, j ) = c( i, j ) + temp*a( i, l )
c;   70                continue
c;                  end if
c;   80          continue
c;   90       continue
c;         else
c;*
c;*           Form  C := alpha*A'*B + beta*C
c;*
c;            do 120, j = 1, n
c;               do 110, i = 1, m
c;                  temp = zero
c;                  do 100, l = 1, k
c;                     temp = temp + a( l, i )*b( l, j )
c;  100             continue
c;                  if( beta.eq.zero )then
c;                     c( i, j ) = alpha*temp
c;                  else
c;                     c( i, j ) = alpha*temp + beta*c( i, j )
c;                  end if
c;  110          continue
c;  120       continue
c;         end if
c;      else
c;         if( nota )then
c;*
c;*           Form  C := alpha*A*B' + beta*C
c;*
c;            do 170, j = 1, n
c;               if( beta.eq.zero )then
c;                  do 130, i = 1, m
c;                     c( i, j ) = zero
c;  130             continue
c;               else if( beta.ne.one )then
c;                  do 140, i = 1, m
c;                     c( i, j ) = beta*c( i, j )
c;  140             continue
c;               end if
c;               do 160, l = 1, k
c;                  if( b( j, l ).ne.zero )then
c;                     temp = alpha*b( j, l )
c;                     do 150, i = 1, m
c;                        c( i, j ) = c( i, j ) + temp*a( i, l )
c;  150                continue
c;                  end if
c;  160          continue
c;  170       continue
c;         else
c;*
c;*           Form  C := alpha*A'*B' + beta*C
c;*
c;            do 200, j = 1, n
c;               do 190, i = 1, m
c;                  temp = zero
c;                  do 180, l = 1, k
c;                     temp = temp + a( l, i )*b( j, l )
c;  180             continue
c;                  if( beta.eq.zero )then
c;                     c( i, j ) = alpha*temp
c;                  else
c;                     c( i, j ) = alpha*temp + beta*c( i, j )
c;                  end if
c;  190          continue
c;  200       continue
c;         end if
c;      end if
c;*
c;      return
c;*
c;*     End of DGEMM .
c;*
c;      end
#endif

