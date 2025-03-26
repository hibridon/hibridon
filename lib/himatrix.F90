module mod_himatrix

!   5. transp         subroutine to carry out in place transposition    *
!                     of n x n matrix a                                 *

contains
! NB #if defined(HIB_UNIX_AIX) rather than unix-ibm for fortran not essl routines
#if defined(HIB_UNIX_IBM)
subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
!
integer n,nm,ierr,matz
double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
!
!     this subroutine calls the essl routine dspev to get the
!     eigenvalues and eigenvectors of a real symmetric matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.


!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1  and  fv2  are temporary storage arrays.
!
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
ione = 1
naux=2*nm
! compress matrix into lower triangle
call triang(a,nm,n)
! call essl routine
call dspev(ione,a,w,z,nm,n,fv1,naux)
return
end
subroutine triang(a,nrow,n)
! -------------------
! subroutine to pack lower triangle of symmetric matrix, stored as
! column form
! input
!   a:         on input: full matrix stored column by column
!              on return: packed lower triangle
!   nrow:      maximum row dimension of matrix
!   n:         order of matrix
! written by:  millard alexander
! current revision date;  1-aug-91
! -------------------
implicit double precision (a-h,o-z)
dimension a(25)
if (n .gt. 1) then
  indnew=n+1
  indold=nrow+2
  do  50 i = 2, n
    nelt=n-i+1
    call dcopy(nelt,a(indold),1,a(indnew),1)
    indold=indold+nrow+1
    indnew=indnew+nelt
50   continue
endif
return
end
#endif
!#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
#if defined(HIB_UNIX) && !defined(HIB_UNIX_IBM) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
subroutine rsg(nm,n,a,b,w,matz,z,fv1,fv2,ierr)
!
integer n,nm,ierr,matz
double precision a(nm,n),b(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the real symmetric generalized eigenproblem  ax = (lambda)bx.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrices  a  and  b.
!
!        a  contains a real symmetric matrix.
!
!        b  contains a positive definite real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1  and  fv2  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
if (n .le. nm) go to 10
ierr = 10 * n
go to 50
!
10 call  reduc(nm,n,a,b,fv2,ierr)
if (ierr .ne. 0) go to 50
if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
call  tred1(nm,n,a,w,fv1,fv2)
call  tqlrat(n,w,fv2,ierr)
go to 50
!     .......... find both eigenvalues and eigenvectors ..........
20 call  tred2(nm,n,a,w,fv1,z)
call  tql2(nm,n,w,fv1,z,ierr)
if (ierr .ne. 0) go to 50
call  rebak(nm,n,b,fv2,n,z)
50 return
end
subroutine rebak(nm,n,b,dl,m,z)
!
integer i,j,k,m,n,i1,ii,nm
double precision b(nm,n),dl(n),z(nm,m)
double precision x
!
!     this subroutine is a translation of the algol procedure rebaka,
!     num. math. 11, 99-110(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
!
!     this subroutine forms the eigenvectors of a generalized
!     symmetric eigensystem by back transforming those of the
!     derived symmetric matrix determined by  reduc.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix system.
!
!        b contains information about the similarity transformation
!          (cholesky decomposition) used in the reduction by  reduc
!          in its strict lower triangle.
!
!        dl contains further information about the transformation.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
if (m .eq. 0) go to 200
!
do 100 j = 1, m
!     .......... for i=n step -1 until 1 do -- ..........
   do 100 ii = 1, n
      i = n + 1 - ii
      i1 = i + 1
      x = z(i,j)
      if (i .eq. n) go to 80
!
      do 60 k = i1, n
60       x = x - b(k,i) * z(k,j)
!
80       z(i,j) = x / dl(i)
100 continue
!
200 return
end
subroutine reduc(nm,n,a,b,dl,ierr)
!
integer i,j,k,n,i1,j1,nm,nn,ierr
double precision a(nm,n),b(nm,n),dl(n)
double precision x,y
!
!     this subroutine is a translation of the algol procedure reduc1,
!     num. math. 11, 99-110(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
!
!     this subroutine reduces the generalized symmetric eigenproblem
!     ax=(lambda)bx, where b is positive definite, to the standard
!     symmetric eigenproblem using the cholesky factorization of b.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices a and b.  if the cholesky
!          factor l of b is already available, n should be prefixed
!          with a minus sign.
!
!        a and b contain the real symmetric input matrices.  only the
!          full upper triangles of the matrices need be supplied.  if
!          n is negative, the strict lower triangle of b contains,
!          instead, the strict lower triangle of its cholesky factor l.
!
!        dl contains, if n is negative, the diagonal elements of l.
!
!     on output
!
!        a contains in its full lower triangle the full lower triangle
!          of the symmetric matrix derived from the reduction to the
!          standard form.  the strict upper triangle of a is unaltered.
!
!        b contains in its strict lower triangle the strict lower
!          triangle of its cholesky factor l.  the full upper
!          triangle of b is unaltered.
!
!        dl contains the diagonal elements of l.
!
!        ierr is set to
!          zero       for normal return,
!          7*n+1      if b is not positive definite.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
ierr = 0
nn = iabs(n)
if (n .lt. 0) go to 100
!     .......... form l in the arrays b and dl ..........
do 80 i = 1, n
   i1 = i - 1
!
   do 80 j = i, n
      x = b(i,j)
      if (i .eq. 1) go to 40
!
      do 20 k = 1, i1
20       x = x - b(i,k) * b(j,k)
!
40       if (j .ne. i) go to 60
      if (x .le. 0.0d0) go to 1000
      y = dsqrt(x)
      dl(i) = y
      go to 80
60       b(j,i) = x / y
80 continue
!     .......... form the transpose of the upper triangle of inv(l)*a
!                in the lower triangle of the array a ..........
100 do 200 i = 1, nn
   i1 = i - 1
   y = dl(i)
!
   do 200 j = i, nn
      x = a(i,j)
      if (i .eq. 1) go to 180
!
      do 160 k = 1, i1
160       x = x - b(i,k) * a(j,k)
!
180       a(j,i) = x / y
200 continue
!     .......... pre-multiply by inv(l) and overwrite ..........
do 300 j = 1, nn
   j1 = j - 1
!
   do 300 i = j, nn
      x = a(i,j)
      if (i .eq. j) go to 240
      i1 = i - 1
!
      do 220 k = j, i1
220       x = x - a(k,j) * b(i,k)
!
240       if (j .eq. 1) go to 280
!
      do 260 k = 1, j1
260       x = x - a(j,k) * b(i,k)
!
280       a(i,j) = x / dl(i)
300 continue
!
go to 1001
!     .......... set error -- b is not positive definite ..........
1000 ierr = 7 * n + 1
1001 return
end
!** from netlib, tue may 24 12:02:59 cdt 1988 ***
!  eispack eigenvalue package
subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
!
integer n,nm,ierr,matz
double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a real symmetric matrix.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrix  a.
!
!        a  contains the real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1  and  fv2  are temporary storage arrays.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
if (n .le. nm) go to 10
ierr = 10 * n
go to 50
!
10 if (matz .ne. 0) go to 20
!     .......... find eigenvalues only ..........
call  tred1(nm,n,a,w,fv1,fv2)
!  tqlrat encounters catastrophic underflow on the vax
call  tqlrat(n,w,fv2,ierr)
!     call  tql1(n,w,fv1,ierr)
go to 50
!     .......... find both eigenvalues and eigenvectors ..........
20 call  tred2(nm,n,a,w,fv1,z)
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_IBM) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
call  tql2(nm,n,w,fv1,z,ierr)
50 return
end
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
subroutine tql2(nm,n,d,e,z,ierr)
!
integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
double precision d(n),e(n),z(nm,n)
double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag
!
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
ierr = 0
if (n .eq. 1) go to 1001
!
do 100 i = 2, n
100 e(i-1) = e(i)
!
f = 0.0d0
tst1 = 0.0d0
e(n) = 0.0d0
!
do 240 l = 1, n
   j = 0
   h = dabs(d(l)) + dabs(e(l))
   if (tst1 .lt. h) tst1 = h
!     .......... look for small sub-diagonal element ..........
   do 110 m = l, n
      tst2 = tst1 + dabs(e(m))
      if (tst2 .eq. tst1) go to 120
!     .......... e(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
110    continue
!
120    if (m .eq. l) go to 220
130    if (j .eq. 30) go to 1000
   j = j + 1
!     .......... form shift ..........
   l1 = l + 1
   l2 = l1 + 1
   g = d(l)
   p = (d(l1) - g) / (2.0d0 * e(l))
   r = pythag(p,1.0d0)
   d(l) = e(l) / (p + dsign(r,p))
   d(l1) = e(l) * (p + dsign(r,p))
   dl1 = d(l1)
   h = g - d(l)
   if (l2 .gt. n) go to 145
!
   do 140 i = l2, n
140    d(i) = d(i) - h
!
145    f = f + h
!     .......... ql transformation ..........
   p = d(m)
   c = 1.0d0
   c2 = c
   el1 = e(l1)
   s = 0.0d0
   mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
!        do 200 ii = 1, mml
   do 200 i=m-1,l,-1
      c3 = c2
      c2 = c
      s2 = s
!           i = m - ii
      g = c * e(i)
      h = c * p
      r = pythag(p,e(i))
      e(i+1) = s * r
      s = e(i) / r
      c = p / r
      p = c * d(i) - s * g
      d(i+1) = h + s * (c * g + s * d(i))
!  replace this inner loop with blas call (millard alexander 12/24/90)
!     .......... form vector ..........
!            do 180 k = 1, n
!               h = z(k,i+1)
!               hh = z(k, i)
!               z(k,i+1) = s * hh + c * h
!               z(k,i) = c * hh - s * h
!  180       continue
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
     call drot (n, z(1,i+1), 1, z(1,i), 1, c, s)
#endif
#if defined(HIB_UNIX_CONVEX)
     call srot (n, z(1,i+1), 1, z(1,i), 1, c, s)
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
200    continue
!
   p = -s * s2 * c3 * el1 * e(l) / dl1
   e(l) = s * p
   d(l) = c * p
   tst2 = tst1 + dabs(e(l))
   if (tst2 .gt. tst1) go to 130
220    d(l) = d(l) + f
240 continue
!     .......... order eigenvalues and eigenvectors ..........
do 300 ii = 2, n
   i = ii - 1
   k = i
   p = d(i)
!
   do 260 j = ii, n
      if (d(j) .ge. p) go to 260
      k = j
      p = d(j)
260    continue
!
   if (k .eq. i) go to 300
   d(k) = d(i)
   d(i) = p
!
   do 280 j = 1, n
      p = z(j,i)
      z(j,i) = z(j,k)
      z(j,k) = p
280    continue
!
300 continue
!
go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 return
end
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
subroutine tqlrat(n,d,e2,ierr)
!
integer i,j,l,m,n,l1,mml,ierr
double precision d(n),e2(n)
double precision b,c,f,g,h,p,r,s,t,epslon,pythag
!
!     this subroutine is a translation of the algol procedure tqlrat,
!     algorithm 464, comm. acm 16, 689(1973) by reinsch.
!
!     this subroutine finds the eigenvalues of a symmetric
!     tridiagonal matrix by the rational ql method.
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e2 contains the squares of the subdiagonal elements of the
!          input matrix in its last n-1 positions.  e2(1) is arbitrary.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        e2 has been destroyed.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1987.
!     modified by c. moler to fix underflow/overflow difficulties,
!     especially on the vax and other machines where epslon(1.0d0)**2
!     nearly underflows.  see the loop involving statement 102 and
!     the two statements just before statement 200.
!
!     ------------------------------------------------------------------
!
ierr = 0
if (n .eq. 1) go to 1001
!
do 100 i = 2, n
100 e2(i-1) = e2(i)
!
f = 0.0d0
t = 0.0d0
e2(n) = 0.0d0
!
do 290 l = 1, n
   j = 0
   h = dabs(d(l)) + dsqrt(e2(l))
   if (t .gt. h) go to 105
   t = h
   b = epslon(t)
   c = b * b
   if (c .ne. 0.0d0) go to 105
!        spliting tolerance underflowed.  look for larger value.
   do 102 i = l, n
      h = dabs(d(i)) + dsqrt(e2(i))
      if (h .gt. t) t = h
102    continue
   b = epslon(t)
   c = b * b
!     .......... look for small squared sub-diagonal element ..........
105    do 110 m = l, n
      if (e2(m) .le. c) go to 120
!     .......... e2(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
110    continue
!
120    if (m .eq. l) go to 210
130    if (j .eq. 30) go to 1000
   j = j + 1
!     .......... form shift ..........
   l1 = l + 1
   s = dsqrt(e2(l))
   g = d(l)
   p = (d(l1) - g) / (2.0d0 * s)
   r = pythag(p,1.0d0)
   d(l) = s / (p + dsign(r,p))
   h = g - d(l)
!
   do 140 i = l1, n
140    d(i) = d(i) - h
!
   f = f + h
!     .......... rational ql transformation ..........
   g = d(m)
   if (g .eq. 0.0d0) g = b
   h = g
   s = 0.0d0
   mml = m - l
!     .......... for i=m-1 step -1 until l do -- ..........
!           do 200 ii = 1, mml
!           i = m - ii
      do 200 i=m-1,l,-1
      p = g * h
      r = p + e2(i)
      e2(i+1) = s * r
      s = e2(i) / r
      d(i+1) = h + s * (h + d(i))
      g = d(i) - e2(i) / g
!           avoid division by zero on next pass
      if (g .eq. 0.0d0) g = epslon(d(i))
      h = g * (p / r)
200    continue
!
   e2(l) = s * g
   d(l) = h
!     .......... guard against underflow in convergence test ..........
   if (h .eq. 0.0d0) go to 210
   if (dabs(e2(l)) .le. dabs(c/h)) go to 210
   e2(l) = h * e2(l)
   if (e2(l) .ne. 0.0d0) go to 130
210    p = d(l) + f
!     .......... order eigenvalues ..........
   if (l .eq. 1) go to 250
!     .......... for i=l step -1 until 2 do -- ..........
!        do 230 ii = 2, l
!           i = l + 2 - ii
      do 230 i=l,2,-1
      if (p .ge. d(i-1)) go to 270
      d(i) = d(i-1)
230    continue
!
250    i = 1
270    d(i) = p
290 continue
!
go to 1001
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 return
end
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
subroutine tred1(nm,n,a,d,e,e2)
!
integer i,j,k,l,n,nm,jp1
double precision a(nm,n),d(n),e(n),e2(n)
double precision f,g,h,scale
!
!     this subroutine is a translation of the algol procedure tred1,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix
!     to a symmetric tridiagonal matrix using
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        a contains information about the orthogonal trans-
!          formations used in the reduction in its strict lower
!          triangle.  the full upper triangle of a is unaltered.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
do 100 i = 1, n
   d(i) = a(n,i)
   a(n,i) = a(i,i)
100 continue
!     .......... for i=n step -1 until 1 do -- ..........
!     do 300 ii = 1, n
!        i = n + 1 - ii
   do 300 i=n,1,-1
   l = i - 1
   h = 0.0d0
   scale = 0.0d0
   if (l .lt. 1) go to 130
!     .......... scale row (algol tol then not needed) ..........
   do 120 k = 1, l
120    scale = scale + dabs(d(k))
!
   if (scale .ne. 0.0d0) go to 140
!
   do 125 j = 1, l
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = 0.0d0
125    continue
!
130    e(i) = 0.0d0
   e2(i) = 0.0d0
   go to 300
!
140    do 150 k = 1, l
      d(k) = d(k) / scale
      h = h + d(k) * d(k)
150    continue
!
   e2(i) = scale * scale * h
   f = d(l)
   g = -dsign(dsqrt(h),f)
   e(i) = scale * g
   h = h - f * g
   d(l) = f - g
   if (l .eq. 1) go to 285
!     .......... form a*u ..........
   do 170 j = 1, l
170    e(j) = 0.0d0
!
   do 240 j = 1, l
      f = d(j)
      g = e(j) + a(j,j) * f
      jp1 = j + 1
      if (l .lt. jp1) go to 220
!
      do 200 k = jp1, l
         g = g + a(k,j) * d(k)
         e(k) = e(k) + a(k,j) * f
200       continue
!
220       e(j) = g
240    continue
!     .......... form p ..........
   f = 0.0d0
!
   do 245 j = 1, l
      e(j) = e(j) / h
      f = f + e(j) * d(j)
245    continue
!
   h = f / (h + h)
!     .......... form q ..........
   do 250 j = 1, l
250    e(j) = e(j) - h * d(j)
!     .......... form reduced a ..........
   do 280 j = 1, l
      f = d(j)
      g = e(j)
!
      do 260 k = j, l
260       a(k,j) = a(k,j) - f * e(k) - g * d(k)
!
280    continue
!
285    do 290 j = 1, l
      f = d(j)
      d(j) = a(l,j)
      a(l,j) = a(i,j)
      a(i,j) = f * scale
290    continue
!
300 continue
!
return
end
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
subroutine tred2(nm,n,a,d,e,z)
!
integer i,j,k,l,n,nm,jp1
double precision a(nm,n),d(n),e(n),z(nm,n)
double precision f,g,h,hh,scale
!
!     this subroutine is a translation of the algol procedure tred2,
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix to a
!     symmetric tridiagonal matrix using and accumulating
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e(1) is set to zero.
!
!        z contains the orthogonal transformation matrix
!          produced in the reduction.
!
!        a and z may coincide.  if distinct, a is unaltered.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!
!     ------------------------------------------------------------------
!
do 100 i = 1, n
!
   do 80 j = i, n
80    z(j,i) = a(j,i)
!
   d(i) = a(n,i)
100 continue
!
if (n .eq. 1) go to 510
!     .......... for i=n step -1 until 2 do -- ..........
!     do 300 ii = 2, n
!        i = n + 2 - ii
do 300 i=n,2,-1
   l = i - 1
   h = 0.0d0
   scale = 0.0d0
   if (l .lt. 2) go to 130
!     .......... scale row (algol tol then not needed) ..........
   do 120 k = 1, l
120    scale = scale + dabs(d(k))
!
   if (scale .ne. 0.0d0) go to 140
130    e(i) = d(l)
!
   do 135 j = 1, l
      d(j) = z(l,j)
      z(i,j) = 0.0d0
      z(j,i) = 0.0d0
135    continue
!
   go to 290
!
140    do 150 k = 1, l
      d(k) = d(k) / scale
      h = h + d(k) * d(k)
150    continue
!
   f = d(l)
   g = -dsign(dsqrt(h),f)
   e(i) = scale * g
   h = h - f * g
   d(l) = f - g
!     .......... form a*u ..........
   do 170 j = 1, l
170    e(j) = 0.0d0
!
   do 240 j = 1, l
      f = d(j)
      z(j,i) = f
      g = e(j) + z(j,j) * f
      jp1 = j + 1
      if (l .lt. jp1) go to 220
!
      do 200 k = jp1, l
         g = g + z(k,j) * d(k)
         e(k) = e(k) + z(k,j) * f
200       continue
!
220       e(j) = g
240    continue
!     .......... form p ..........
   f = 0.0d0
!
   do 245 j = 1, l
      e(j) = e(j) / h
      f = f + e(j) * d(j)
245    continue
!
   hh = f / (h + h)
!     .......... form q ..........
   do 250 j = 1, l
250    e(j) = e(j) - hh * d(j)
!     .......... form reduced a ..........
   do 280 j = 1, l
      f = d(j)
      g = e(j)
!
      do 260 k = j, l
260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
!
      d(j) = z(l,j)
      z(i,j) = 0.0d0
280    continue
!
290    d(i) = h
300 continue
!     .......... accumulation of transformation matrices ..........
do 500 i = 2, n
   l = i - 1
   z(n,l) = z(l,l)
   z(l,l) = 1.0d0
   h = d(i)
   if (h .eq. 0.0d0) go to 380
!
   do 330 k = 1, l
330    d(k) = z(k,i) / h
!
   do 360 j = 1, l
      g = 0.0d0
!
      do 340 k = 1, l
340       g = g + z(k,i) * z(k,j)
!
      do 360 k = 1, l
         z(k,j) = z(k,j) - g * d(k)
360    continue
!
380    do 400 k = 1, l
400    z(k,i) = 0.0d0
!
500 continue
!
510 do 520 i = 1, n
   d(i) = z(n,i)
   z(n,i) = 0.0d0
520 continue
!
z(n,n) = 1.0d0
e(1) = 0.0d0
return
end
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
double precision function epslon (x)
double precision x
!
!     estimate unit roundoff in quantities of size x.
!
double precision a,b,c,eps
!
!     this program should function properly on all systems
!     satisfying the following two assumptions,
!        1.  the base used in representing floating point
!            numbers is not a power of three.
!        2.  the quantity  a  in statement 10 is represented to
!            the accuracy used in floating point variables
!            that are stored in memory.
!     the statement number 10 and the go to 10 are intended to
!     force optimizing compilers to generate code satisfying
!     assumption 2.
!     under these assumptions, it should be true that,
!            a  is not exactly equal to four-thirds,
!            b  has a zero for its last bit or digit,
!            c  is not exactly equal to one,
!            eps  measures the separation of 1.0 from
!                 the next larger floating point number.
!     the developers of eispack would appreciate being informed
!     about any systems where these assumptions do not hold.
!
!     this version dated 4/6/83.
!
a = 4.0d0/3.0d0
10 b = a - 1.0d0
c = b + b + b
eps = dabs(c-1.0d0)
if (eps .eq. 0.0d0) go to 10
epslon = eps*dabs(x)
return
end
double precision function pythag(a,b)
double precision a,b
!
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!
double precision p,r,s,t,u
p = dmax1(dabs(a),dabs(b))
if (p .eq. 0.0d0) go to 20
r = (dmin1(dabs(a),dabs(b))/p)**2
10 continue
   t = 4.0d0 + r
   if (t .eq. 4.0d0) go to 20
   s = r/t
   u = 1.0d0 + 2.0d0*s
   p = u*p
   r = (s/u)**2 * r
go to 10
20 pythag = p
return
end
#endif
!**************************************************************************
!                                                                         *
!                          linpack library                                *
!                                                                         *
!**************************************************************************
!                         routines included:                              *
!                                                                         *
!   1.  ssidi     computes determinant, inertia and inverse of a real     *
!                 symmetric matrix using the factors from ssifa           *
!                                                                         *
!   2.  ssifa     factors a real symmetric matrix by elimination          *
!                 with symmetric pivoting.                                *
!                                                                         *
!   3.  sgefa     factors a real matrix by gaussian elimination.          *
!                                                                         *
!   4.  sgesl     solves the real system a * x = b  or  trans(a) * x = b  *
!                 using the factors computed by sgeco or sgefa.           *
!                                                                         *
!**************************************************************************
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
! ------------------------------------------------------------
subroutine sgefa(a,lda,n,ipvt,info)
use mod_hiutil, only: daxpy_wrapper
implicit double precision (a-h,o-z)
integer lda,n,ipvt(1),info
dimension a(lda,1)
!      real a(lda,1)
!
!     sgefa factors a real matrix by gaussian elimination.
!
!     sgefa is usually called by sgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
!
!     on entry
!
!        a       real(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that sgesl or sgedi will divide by zero
!                     if called.  use  rcond  in sgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas saxpy,sscal,isamax
!
!     internal variables
!
!      real t
integer idamax,j,k,kp1,l,nm1,idummy
idummy=1
!
!
!     gaussian elimination with partial pivoting
!
info = 0
nm1 = n - 1
if (nm1 .lt. 1) go to 70
do 60 k = 1, nm1
   kp1 = k + 1
!
!        find l = pivot index
!
   l = idamax(n-k+1,a(k,k),1) + k - 1
   ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
   if (a(l,k) .eq. 0.0) go to 40
!
!           interchange if necessary
!
      if (l .eq. k) go to 10
         t = a(l,k)
         a(l,k) = a(k,k)
         a(k,k) = t
10       continue
!
!           compute multipliers
!
      t = -1.0/a(k,k)
      call dscal(n-k,t,a(k+1,k),idummy)
!
!           row elimination with column indexing
!
      do 30 j = kp1, n
         t = a(l,j)
         if (l .eq. k) go to 20
            a(l,j) = a(k,j)
            a(k,j) = t
20          continue
       call daxpy_wrapper(n-k,t,a(k+1,k),idummy,a(k+1,j),idummy)
30       continue
   go to 50
40    continue
      info = k
50    continue
60 continue
70 continue
ipvt(n) = n
if (a(n,n) .eq. 0.0) info = n
return
end
#endif
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
! ------------------------------------------------
subroutine sgesl(a,lda,n,ipvt,b,job)
! current revision date: 29/09/87
use mod_hiutil, only: daxpy_wrapper
implicit double precision (a-h,o-z)
integer lda,n,ipvt(1),job
dimension a(lda,1),b(1)
!      real a(lda,1),b(1)
!
!     sgesl solves the real system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by sgeco or sgefa.
!
!     on entry
!
!        a       real(lda, n)
!                the output from sgeco or sgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from sgeco or sgefa.
!
!        b       real(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b  where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if sgeco has set rcond .gt. 0.0
!        or sgefa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call sgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call sgesl(a,lda,n,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas saxpy,sdot
!
!     internal variables
!
!      real sdot,t
integer k,kb,l,nm1,idummy
idummy=1
!
nm1 = n - 1
if (job .ne. 0) go to 50
!
!        job = 0 , solve  a * x = b
!        first solve  l*y = b
!
   if (nm1 .lt. 1) go to 30
   do 20 k = 1, nm1
      l = ipvt(k)
      t = b(l)
      if (l .eq. k) go to 10
         b(l) = b(k)
         b(k) = t
10       continue
      call daxpy_wrapper(n-k,t,a(k+1,k),idummy,b(k+1),idummy)
20    continue
30    continue
!
!        now solve  u*x = y
!
   do 40 kb = 1, n
      k = n + 1 - kb
      b(k) = b(k)/a(k,k)
      t = -b(k)
      call daxpy_wrapper(k-1,t,a(1,k),idummy,b(1),idummy)
40    continue
go to 100
50 continue
!
!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b
!
   do 60 k = 1, n
      t = ddot(k-1,a(1,k),idummy,b(1),idummy)
      b(k) = (b(k) - t)/a(k,k)
60    continue
!
!        now solve trans(l)*x = y
!
   if (nm1 .lt. 1) go to 90
   do 80 kb = 1, nm1
      k = n - kb
      b(k) = b(k) + ddot(n-k,a(k+1,k),idummy,b(k+1),idummy)
      l = ipvt(k)
      if (l .eq. k) go to 70
         t = b(l)
         b(l) = b(k)
         b(k) = t
70       continue
80    continue
90    continue
100 continue
return
end
#endif
! -----------------------------------------------------------------------
#if defined(HIB_UNIX_DARWIN) || defined(HIB_UNIX_X86)
subroutine syminv (a,lda,n,ierr)
use mod_hiutil, only: dgetri_wrapper
implicit double precision (a-h,o-z)
!
!     -----------------------------------------------------------------
!     This subroutine uses LAPACK DGETRF and DGETRI
!     to invert a real symmetric indefinite matrix.
!     latest revision:  18-feb-2008 by millard alexander (replace calls to "sy" with "ge"
!         lapack routines
!     -----------------------------------------------------------------
!
dimension a(lda,n)
dimension ipiv(4*n),work(256*n)

!
lwork = 256*n
!  form full matrix
do j = 2,n
   do i = 1,j-1
      a(i,j) = a(j,i)
   enddo
enddo


!      write(6,444) lda,n
!444   format ('calling dgetrf:  lda=',i8,'  n=',i8)


call dgetrf (n,n,a,lda,ipiv,ierr)
if (ierr .ne. 0) stop 'error in dgetrf'
call dgetri_wrapper(n,a,lda,ipiv,work,lwork,ierr)
if (ierr .ne. 0) stop 'error in dgetri'
return
end
#endif
#if (defined(HIB_UNIX) || defined(HIB_CRAY)) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
 subroutine smxinv (a, nmax, n, scr, kpvt, ierr)
!  subroutine to invert the real symmetric matrix a using
!  lapack routines dsytrf and dsytri
!  written by:  millard alexander
!  latest revision date:  13-nov-1995
! -----------------------------------------------------------------------
!  variables in call list:
!  nmax:   maximum row dimension of matrices
!  n:      actual order of matrices
!  a:      symmetric matrix of order n x n, stored in packed column form
!          on input:  lower triangle contains matrix to be inverted
!          on return: contains inverse(a), both lower and upper triangles
!  scr:    scratch vector of length n
!  kpvt:   scratch vector of length n
!  ierr:   on return:  set equal to 0 if normal return
!                      set equal to nn if singularity due to row nn of matri
!  subroutines used:
!  dsytrf, dsytri:   lapack routines to factor and to invert a symmetric, re
!                  matrix
! -----------------------------------------------------------------------
#if defined(HIB_UNIX) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
use mod_cosc11, only: sc11
#endif
#if defined(HIB_UNIX_IBM)
use mod_cokaux, only kaux => naux
#endif
implicit double  precision (a-h,o-z)
integer icol, icolpt, ierr, ione, irowpt, izero, n, ncol, &
        nmax, nmaxp1
integer kpvt
dimension a(1), scr(1), kpvt(1)
#endif
#if defined(HIB_UNIX_IBM)
dimension det(2)
#endif
#if (defined(HIB_UNIX) || defined(HIB_CRAY)) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
data ione, izero  /1, 0/
#endif
#if defined(HIB_UNIX_IBM)
naux=kaux
#endif
#if (defined(HIB_CRAY) || defined(HIB_UNIX)) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
ierr = izero
!  first fill in upper half of original matrix
nmaxp1 = nmax + 1
icolpt = 2
irowpt = nmaxp1
do  50  icol = 1, n - 1
!  icolpt points to first sub-diagonal element in column icol
!  irowpt points to first super-diagonal element in row icol
!  ncol is number of subdiagonal elements in column icol
  ncol = n - icol
  call dcopy (ncol, a(icolpt), 1, a(irowpt), nmax)
  icolpt = icolpt + nmaxp1
  irowpt = irowpt + nmaxp1
50 continue
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
 lwork=64*nmax
 lwork=-1
 call dsytrf('U',n,a,nmax,kpvt,sc11,lwork,info)
 print *, ' n, info, work(1):  ', n, info, sc11(1)
 stop 'after dsytrf'
 if (info .lt. 0) then
   write (9, 54) info
   write (6, 54) info
54    format (' *** ERROR IN',i4, &
          'TH VARIABLE IN DSYTRF ***')
   stop
 endif
 if (info .gt. 0) then
   write (9, 55) info
   write (6, 55) info
55    format (' *** POSSIBLE SINGULARITY IN',i4, &
          'TH ROW OF MATRIX TO BE INVERTED ***')
   stop
 end if
#endif
#if defined(HIB_UNIX_IBM)
call dgeicd(a, nmax, n,izero,rcont,det,sc11,naux)
#endif
! old linpack calls
#if defined(HIB_CRAY)
call ssifa (a, nmax, n, kpvt, ierr)
#endif
#if defined(HIB_UNIX_CONVEX)
call dsifa (a, nmax, n, kpvt, ierr)
#endif
#if defined(HIB_CRAY) || defined(HIB_UNIX_CONVEX)
if (ierr .ne. 0) then
  write (9, 60) ierr
  write (6, 60) ierr
60   format (' *** POSSIBLE SINGULARITY IN',i4, &
          'TH ROW OF MATRIX TO BE INVERTED ***')
  ierr = 0
end if
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)
 call dsytri('U',n,a,nmax,kpvt,sc11,info)
 if (info .lt. 0) then
   write (9, 64) info
   write (6, 64) info
64    format (' *** ERROR IN',i4, &
          'TH VARIABLE IN DSYTRI ***')
   stop
 endif
 if (info .gt. 0) then
   write (9, 65) info
   write (6, 65) info
65    format (' *** POSSIBLE SINGULARITY IN',i4, &
          'TH ROW OF MATRIX TO BE INVERTED ***')
   stop
 end if
#endif
#if defined(HIB_CRAY)
call ssidi (a, nmax, n, kpvt, det, inert, scr, ione)
#endif
#if defined(HIB_UNIX_CONVEX)
call dsidi (a, nmax, n, kpvt, det, inert, scr, ione)
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_CRAY) || defined(HIB_UNIX_CONVEX) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_IRIS)  || defined(HIB_UNIX_S)
!  inverse is now in upper triangle of matrix a
!  copy upper triangle of inverse back into lower triangle of a
icolpt = 2
irowpt = nmaxp1
do 100 icol = 1, n - 1
!  icolpt points to first sub-diagonal element in column icol
!  irowpt points to first super-diagonal element in row icol
!  ncol is number of subdiagonal elements in column icol
  ncol = n - icol
#endif
#if defined(HIB_UNIX_CONVEX)
call scopy (ncol, a(irowpt), nmax, a(icolpt), 1)
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_CRAY) || defined(HIB_UNIX_CONVEX) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_IRIS)  || defined(HIB_UNIX_S)
  call dcopy (ncol, a(irowpt), nmax, a(icolpt), 1)
#endif
#if defined(HIB_UNIX_AIX) || defined(HIB_CRAY) || defined(HIB_UNIX_CONVEX) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_IRIS)  || defined(HIB_UNIX_S)
  icolpt = icolpt + nmaxp1
  irowpt = irowpt + nmaxp1
100 continue
#endif
#if (defined(HIB_CRAY) || defined(HIB_UNIX)) && !defined(HIB_UNIX_DARWIN) && !defined(HIB_UNIX_X86)
  return
  end
#endif
!-----------------------------
#if defined(HIB_NONE)
subroutine dsytri( uplo, n, a, lda, ipiv, work, info )
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
external           dcopy, dswap, dsymv, xerbla
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
integer            idamax
external           lsame, idamax
!     ..
!     .. External Subroutines ..
external           dlaev2, drot, dscal, dswap, dsyr, xerbla
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
#if !defined(HIB_CRAY) && !defined(HIB_UNIX_CONVEX)

!--------------------------------------------------------------------
subroutine mxma(a,mcola,mrowa,b,mcolb,mrowb, &
                r,mcolr,mrowr,ncol,nlink,nrow)
!--------------------------------------------------------------------
!.....unrolled version for ibm6000 and other risc machines
! subroutine mxma(a,iac,iar, b,ibc,ibr, c,icc,icr, nar,nac,nbc)
! fortran version of cray mxma to calculate c=a*b.
! note: *** this code should not be modified ***
!
! arguments:
! a = first input matrix with effective dimensions a(nar,nac).
! iac = spacing between consecutive elements in a column of a.
! iar = spacing between consecutive elements in a row of a.
! b = second input matrix with effective dimensions b(nac,nbc).
! ibc = spacing between consecutive elements in a column of b.
! ibr = spacing between consecutive elements in a row of b.
! c = output matrix with effective dimensions c(nar,nbc).
! icc = spacing between consecutive elements in a column of c.
! icr = spacing between consecutive elements in a row of c.
! nar = number of rows in a and c.
! nac = number of columns in a and rows in b.
! nbc = number of columns in b and c.
!
! besides allowing arbitrary sub-blocking and element spacing of
! any of the matrices, this routine also allows the transposition
! of any of the matrix arguments.
!
! for example, suppose that z(*) is dimensioned in the calling
! program as:
!
! dimension z(n,*)
!
! then the following call will calculate a matrix product
! involving z(*):
!
! call mxma(...z,1,n, ...)
!
! while a matrix product involving z(transpose) may be calculated
! with:
! call mxma(...z,n,1, ...)
!
! written 04-jan-85 by ron shepard.
!
use mod_comxm, only: ncache, mxmblk
implicit double precision (a-h,o-z)
#include "common/vax1.F90"
#include "common/vax2.F90"
real(8), dimension(ncol*nlink*mcola*mrowa), intent(in) :: a
integer, intent(in) :: mcola
integer, intent(in) :: mrowa
real(8), dimension(nlink*nrow*mcolb*mrowb), intent(in) :: b
integer, intent(in) :: mcolb
integer, intent(in) :: mrowb
real(8), dimension(ncol*nrow*mcolr*mrowr), intent(out) :: r
integer, intent(in) :: mcolr
integer, intent(in) :: mrowr
integer, intent(in) :: ncol
integer, intent(in) :: nlink
integer, intent(in) :: nrow
if(ncol.eq.0.or.nrow.eq.0) return
if(nlink.eq.0) then
  ijj=1
  do 6200 j=1,nrow
  ij=ijj
  do 6100 i=1,ncol
  r(ij)=0
6100   ij=ij+mcolr
6200   ijj=ijj+mrowr
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
6110   ij=ij+mcolr
  ijj=ijj+mrowr
6210   jj=jj+mrowb
  return
end if
if(ncol.gt.5.or.nrow.gt.5) goto 6000
goto (1000,2000,3000,4000,5000),nrow
!
!....nrow=1
1000 goto (1001,1002,1003,1004,1005),ncol
1001   ia1=1
  ib1=1
  s1=0
  do 1010 k=1,nlink
  s1=a(ia1)*b(ib1)+s1
  ia1=ia1+mrowa
1010   ib1=ib1+mcolb
  r(1)=s1
  return
!
1002   ia1=1
  ib1=1
  ia2=1+mcola
  s1=0
  s2=0
  do 1020 k=1,nlink
  s1=a(ia1)*b(ib1)+s1
  ia1=ia1+mrowa
  s2=a(ia2)*b(ib1)+s2
  ia2=ia2+mrowa
1020   ib1=ib1+mcolb
  r(1)=s1
  r(1+mcolr)=s2
  return
!
1003   ia1=1
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
1030   ib1=ib1+mcolb
  r(1)=s1
  r(1+mcolr)=s2
  r(1+2*mcolr)=s3
  return
!
1004   ia1=1
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
1040   ib1=ib1+mcolb
  r(1)=s1
  r(1+mcolr)=s2
  r(1+2*mcolr)=s3
  r(1+3*mcolr)=s4
  return
!
1005   ia1=1
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
1050   ib1=ib1+mcolb
  r(1)=s1
  r(1+mcolr)=s2
  r(1+2*mcolr)=s3
  r(1+3*mcolr)=s4
  r(1+4*mcolr)=s5
  return
!
!....nrow=2
2000 goto (2001,2002,2003,2004,2005),ncol
2001   ia1=1
  ib1=1
  ib2=1+mrowb
  s1=0
  s2=0
  do 2010 k=1,nlink
  s1=a(ia1)*b(ib1)+s1
  ib1=ib1+mcolb
  s2=a(ia1)*b(ib2)+s2
  ib2=ib2+mcolb
2010   ia1=ia1+mrowa
  r(1)=s1
  r(1+mrowr)=s2
  return
!
2002   ia1=1
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
2020   ia2=ia2+mrowa
  r(1)=s1
  r(1+mcolr)=s2
  r(1+mrowr)=s3
  r(1+mcolr+mrowr)=s4
  return
!
2003   ia1=1
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
2030   ib2=ib2+mcolb
  r(1)=s1
  r(1+mcolr)=s2
  r(1+2*mcolr)=s3
  r(1+mrowr)=s4
  r(1+mrowr+mcolr)=s5
  r(1+mrowr+2*mcolr)=s6
  return
!
2004   ia1=1
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
2040   ib2=ib2+mcolb
  r(1)=s11
  r(1+mcolr)=s21
  r(1+2*mcolr)=s31
  r(1+3*mcolr)=s41
  r(1+mrowr)=s12
  r(1+mrowr+mcolr)=s22
  r(1+mrowr+2*mcolr)=s32
  r(1+mrowr+3*mcolr)=s42
  return
!
2005   ia1=1
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
2050   ib2=ib2+mcolb
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
!
!....nrow=3
3000 goto (3001,3002,3003,3004,3005),ncol
3001   ia1=1
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
3010   ib3=ib3+mcolb
  r(1)=s1
  r(1+mrowr)=s2
  r(1+2*mrowr)=s3
  return
!
3002   ib1=1
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
3020   ia2=ia2+mrowa
  r(1)=s1
  r(1+mrowr)=s2
  r(1+2*mrowr)=s3
  r(1+mcolr)=s4
  r(1+mcolr+mrowr)=s5
  r(1+mcolr+2*mrowr)=s6
  return
!
3003   ib1=1
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
3030   ib3=ib3+mcolb
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
!
3004   ib1=1
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
3040   ib3=ib3+mcolb
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
!
3005   ib1=1
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
3050   ib3=ib3+mcolb
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
!
!.....nrow=4
4000 goto(4001,4002,4003,4004,4005),ncol
4001   ia1=1
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
4010   ia1=ia1+mrowa
  r(1)=s1
  r(1+mrowr)=s2
  r(1+2*mrowr)=s3
  r(1+3*mrowr)=s4
  return
!
4002   ib1=1
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
4020   ia2=ia2+mrowa
  r(1)=s11
  r(1+mrowr)=s12
  r(1+2*mrowr)=s13
  r(1+3*mrowr)=s14
  r(1+mcolr)=s21
  r(1+mcolr+mrowr)=s22
  r(1+mcolr+2*mrowr)=s23
  r(1+mcolr+3*mrowr)=s24
  return
!
4003   ib1=1
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
4030   ia3=ia3+mrowa
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
!
4004   ib1=1
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
4040   ia4=ia4+mrowa
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
!
4005   ia1=1
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
4050   continue
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
!
!.....nrow=5
5000 goto(5001,5002,5003,5004,5004),ncol
5001   ia1=1
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
5010   ia1=ia1+mrowa
  r(1)=s1
  r(1+mrowr)=s2
  r(1+2*mrowr)=s3
  r(1+3*mrowr)=s4
  r(1+4*mrowr)=s5
  return
!
5002   ib1=1
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
5020   ia2=ia2+mrowa
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
!
5003   ib1=1
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
5030   ia3=ia3+mrowa
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
!
5004   ia1=1
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
5040   ib1=ib1+mcolb
  ib2=ib2+mcolb
  ib3=ib3+mcolb
  ib4=ib4+mcolb
5050   ib5=ib5+mcolb
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
!
6000 continue
if(nlink*(ncol+nrow)+ncol*nrow.le.ncache) then
  call mxma4(a,mcola,mrowa,b,mcolb,mrowb, &
             r,mcolr,mrowr,ncol,nlink,nrow)
  return
end if
#endif
#if defined(HIB_UNIX_BLAS3)
if(mcolr.eq.1) then
  if(mcola.eq.1.and.mcolb.eq.1) then
    call dgemm('N','N',ncol,nrow,nlink,1.0d0,a,mrowa, &
                    b,max(nlink,mrowb),0.0d0,r,mrowr)
    return
  else if(mrowa.eq.1.and.mcolb.eq.1) then
    call dgemm('T','N',ncol,nrow,nlink,1.0d0,a,mcola, &
                    b,max(nlink,mrowb),0.0d0,r,mrowr)
    return
  else if(mcola.eq.1.and.mrowb.eq.1) then
    call dgemm('N','T',ncol,nrow,nlink,1.0d0,a,mrowa, &
                    b,mcolb,0.0d0,r,mrowr)
    return
  else if(mrowa.eq.1.and.mrowb.eq.1) then
    call dgemm('T','T',ncol,nrow,nlink,1.0d0,a,mcola, &
                    b,mcolb,0.0d0,r,mrowr)
    return
  end if
else if(mrowr.eq.1) then
  if(mcola.eq.1.and.mcolb.eq.1) then
    call dgemm('T','T',nrow,ncol,nlink,1.0d0,b,max(nlink,mrowb), &
                     a,mrowa,0.0d0,r,mcolr)
    return
  else if(mrowa.eq.1.and.mcolb.eq.1) then
    call dgemm('T','N',nrow,ncol,nlink,1.0d0,b,max(nlink,mrowb), &
                     a,mcola,0.0d0,r,mcolr)
    return
  else if(mcola.eq.1.and.mrowb.eq.1) then
    call dgemm('N','T',nrow,ncol,nlink,1.0d0,b,mcolb, &
                     a,mrowa,0.0d0,r,mcolr)
    return
  else if(mrowa.eq.1.and.mrowb.eq.1) then
    call dgemm('N','N',nrow,ncol,nlink,1.0d0,b,mcolb, &
                     a,mcola,0.0d0,r,mcolr)
    return
  end if
end if
#endif
#if !defined(HIB_CRAY) && !defined(HIB_UNIX_CONVEX)
mxb=mxmblk
!      nkb=mxmbln
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
  call mxma4(a(ia1),mcola,mrowa,b(ib1),mcolb,mrowb,r(ir1), &
             mcolr, mrowr,ni,nk,nj)
else
  call mxmb4(a(ia1),mcola,mrowa,b(ib1),mcolb,mrowb,r(ir1), &
             mcolr,mrowr, ni,nk,nj)
end if
80 ie=ie+ni
70 je=je+nj
60 ke=ke+nk
return
6001 ke=0
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
  call mxma4(a(ia1),mcola,mrowa,b(ib1),mcolb,mrowb,r(ir1), &
             mcolr, mrowr,ni,nk,nj)
else
  call mxmb4(a(ia1),mcola,mrowa,b(ib1),mcolb,mrowb,r(ir1), &
             mcolr,mrowr, ni,nk,nj)
end if
71 je=je+nj
81 ie=ie+ni
61 ke=ke+nk
return
end
!--------------------------------------------------------------------
subroutine mxma4(a,mcola,mrowa,b,mcolb,mrowb, &
r,mcolr,mrowr,ncol,nlink,nrow)
!--------------------------------------------------------------------
implicit double precision (a-h,o-z)
dimension r(1),a(1),b(1)
!... r(ncol,nrow)=a(ncol,nlink)*b(nlink,nrow) matrix mult
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
100 ia4=ia4+mrowa
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
200 continue
if(nrest.eq.0) goto 300
goto (201,202,203,300,205),nrest
201   ib1=ibb
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
101   ia4=ia4+mrowa
  r(ir1)=s11
  r(ir2)=s21
  r(ir3)=s31
  r(ir4)=s41
  goto 300
!
202   ib1=ibb
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
102   ia4=ia4+mrowa
  r(ir1)=s11
  r(ir2)=s21
  r(ir3)=s31
  r(ir4)=s41
  r(ir1+mrowr)=s12
  r(ir2+mrowr)=s22
  r(ir3+mrowr)=s32
  r(ir4+mrowr)=s42
  goto 300
!
203   ib1=ibb
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
103   ia4=ia4+mrowa
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
!
205   ib1=ibb
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
105   ia4=ia4+mrowa
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
300 iaa=iaa+4*mcola
if(ncest.eq.0) return
goto (301,302,303,304,305),ncest
301   ibb=1
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
1001   ia1=ia1+mrowa
  r(ir1)=s11
  r(ir1+mrowr)=s12
  ir1=ir1+2*mrowr
2001   continue
  if(mod(nrow,2).eq.0) return
  s11=0
  do 1011 k=1,nlink
  s11=a(iaa)*b(ibb)+s11
  ibb=ibb+mcolb
1011   iaa=iaa+mrowa
  r(ir1)=s11
  return
!
302   ibb=1
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
1002   ia2=ia2+mrowa
  r(ir1)=s11
  r(ir2)=s21
  r(ir1+mrowr)=s12
  r(ir2+mrowr)=s22
  ir1=ir1+2*mrowr
  ir2=ir2+2*mrowr
2002   continue
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
1012   ibb=ibb+mcolb
  r(ir1)=s11
  r(ir2)=s21
  return
!
303   ibb=1
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
1003   ia3=ia3+mrowa
  r(ir1)=s11
  r(ir2)=s21
  r(ir3)=s31
  r(ir1+mrowr)=s12
  r(ir2+mrowr)=s22
  r(ir3+mrowr)=s32
  ir1=ir1+2*mrowr
  ir2=ir2+2*mrowr
  ir3=ir3+2*mrowr
2003   continue
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
1013   ia3=ia3+mrowa
  r(ir1)=s11
  r(ir2)=s21
  r(ir3)=s31
304   return
305   ibb=1
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
1005   ia5=ia5+mrowa
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
2005   continue
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
1015   continue
  r(ir1)=s11
  r(ir2)=s21
  r(ir3)=s31
  r(ir4)=s41
  r(ir5)=s51
  return
end
!--------------------------------------------------------------------
subroutine mxmb(a,mcola,mrowa,b,mcolb,mrowb, &
                r,mcolr,mrowr,ncol,nlink,nrow)
!--------------------------------------------------------------------
!.....unrolled version for ibm6000 and other risc machines
use mod_comxm, only: ncache, mxmblk
implicit double precision (a-h,o-z)
#include "common/vax1.F90"
dimension r(*),a(*),b(*)
#include "common/vax2.F90"
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
6110   ij=ij+mcolr
  ijj=ijj+mrowr
6210   jj=jj+mrowb
  return
end if
if(ncol.gt.5.or.nrow.gt.5) goto 6000
goto (1000,2000,3000,4000,5000),nrow
!
!....nrow=1
1000 goto (1001,1002,1003,1004,1005),ncol
1001   ia1=1
  ib1=1
  s1=0
  do 1010 k=1,nlink
  s1=a(ia1)*b(ib1)+s1
  ia1=ia1+mrowa
1010   ib1=ib1+mcolb
  r(1)=r(1)+s1
  return
!
1002   ia1=1
  ib1=1
  ia2=1+mcola
  s1=0
  s2=0
  do 1020 k=1,nlink
  s1=a(ia1)*b(ib1)+s1
  ia1=ia1+mrowa
  s2=a(ia2)*b(ib1)+s2
  ia2=ia2+mrowa
1020   ib1=ib1+mcolb
  r(1)=r(1)+s1
  r(1+mcolr)=r(1+mcolr)+s2
  return
!
1003   ia1=1
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
1030   ib1=ib1+mcolb
  r(1)=r(1)+s1
  r(1+mcolr)=r(1+mcolr)+s2
  r(1+2*mcolr)=r(1+2*mcolr)+s3
  return
!
1004   ia1=1
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
1040   ib1=ib1+mcolb
  r(1)=r(1)+s1
  r(1+mcolr)=r(1+mcolr)+s2
  r(1+2*mcolr)=r(1+2*mcolr)+s3
  r(1+3*mcolr)=r(1+3*mcolr)+s4
  return
!
1005   ia1=1
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
1050   ib1=ib1+mcolb
  r(1)=r(1)+s1
  r(1+mcolr)=r(1+mcolr)+s2
  r(1+2*mcolr)=r(1+2*mcolr)+s3
  r(1+3*mcolr)=r(1+3*mcolr)+s4
  r(1+4*mcolr)=r(1+4*mcolr)+s5
  return
!
!....nrow=2
2000 goto (2001,2002,2003,2004,2005),ncol
2001   ia1=1
  ib1=1
  ib2=1+mrowb
  s1=0
  s2=0
  do 2010 k=1,nlink
  s1=a(ia1)*b(ib1)+s1
  ib1=ib1+mcolb
  s2=a(ia1)*b(ib2)+s2
  ib2=ib2+mcolb
2010   ia1=ia1+mrowa
  r(1)=r(1)+s1
  r(1+mrowr)=r(1+mrowr)+s2
  return
!
2002   ia1=1
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
2020   ia2=ia2+mrowa
  r(1)=r(1)+s1
  r(1+mcolr)=r(1+mcolr)+s2
  r(1+mrowr)=r(1+mrowr)+s3
  r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s4
  return
!
2003   ia1=1
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
2030   ib2=ib2+mcolb
  r(1)=r(1)+s1
  r(1+mcolr)=r(1+mcolr)+s2
  r(1+2*mcolr)=r(1+2*mcolr)+s3
  r(1+mrowr)=r(1+mrowr)+s4
  r(1+mrowr+mcolr)=r(1+mrowr+mcolr)+s5
  r(1+mrowr+2*mcolr)=r(1+mrowr+2*mcolr)+s6
  return
!
2004   ia1=1
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
2040   ib2=ib2+mcolb
  r(1)=r(1)+s11
  r(1+mcolr)=r(1+mcolr)+s21
  r(1+2*mcolr)=r(1+2*mcolr)+s31
  r(1+3*mcolr)=r(1+3*mcolr)+s41
  r(1+mrowr)=r(1+mrowr)+s12
  r(1+mrowr+mcolr)=r(1+mrowr+mcolr)+s22
  r(1+mrowr+2*mcolr)=r(1+mrowr+2*mcolr)+s32
  r(1+mrowr+3*mcolr)=r(1+mrowr+3*mcolr)+s42
  return
!
2005   ia1=1
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
2050   ib2=ib2+mcolb
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
!
!....nrow=3
3000 goto (3001,3002,3003,3004,3005),ncol
3001   ia1=1
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
3010   ia1=ia1+mrowa
  r(1)=r(1)+s1
  r(1+mrowr)=r(1+mrowr)+s2
  r(1+2*mrowr)=r(1+2*mrowr)+s3
  return
!
3002   ib1=1
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
3020   ia2=ia2+mrowa
  r(1)=r(1)+s1
  r(1+mrowr)=r(1+mrowr)+s2
  r(1+2*mrowr)=r(1+2*mrowr)+s3
  r(1+mcolr)=r(1+mcolr)+s4
  r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s5
  r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s6
  return
!
3003   ib1=1
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
3030   ib3=ib3+mcolb
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
!
3004   ib1=1
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
3040   ib3=ib3+mcolb
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
!
3005   ib1=1
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
3050   ib3=ib3+mcolb
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
!
!.....nrow=4
4000 goto(4001,4002,4003,4004,4005),ncol
4001   ia1=1
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
4010   ia1=ia1+mrowa
  r(1)=r(1)+s1
  r(1+mrowr)=r(1+mrowr)+s2
  r(1+2*mrowr)=r(1+2*mrowr)+s3
  r(1+3*mrowr)=r(1+3*mrowr)+s4
  return
!
4002   ib1=1
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
4020   ia2=ia2+mrowa
  r(1)=r(1)+s11
  r(1+mrowr)=r(1+mrowr)+s12
  r(1+2*mrowr)=r(1+2*mrowr)+s13
  r(1+3*mrowr)=r(1+3*mrowr)+s14
  r(1+mcolr)=r(1+mcolr)+s21
  r(1+mcolr+mrowr)=r(1+mcolr+mrowr)+s22
  r(1+mcolr+2*mrowr)=r(1+mcolr+2*mrowr)+s23
  r(1+mcolr+3*mrowr)=r(1+mcolr+3*mrowr)+s24
  return
!
4003   ib1=1
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
4030   ia3=ia3+mrowa
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
!
4004   ib1=1
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
4040   ia4=ia4+mrowa
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
!
4005   ia1=1
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
4050   continue
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
!
!.....nrow=5
5000 goto(5001,5002,5003,5004,5004),ncol
5001   ia1=1
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
5010   ia1=ia1+mrowa
  r(1)=r(1)+s1
  r(1+mrowr)=r(1+mrowr)+s2
  r(1+2*mrowr)=r(1+2*mrowr)+s3
  r(1+3*mrowr)=r(1+3*mrowr)+s4
  r(1+4*mrowr)=r(1+4*mrowr)+s5
  return
!
5002   ib1=1
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
5020   ia2=ia2+mrowa
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
!
5003   ib1=1
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
5030   ia3=ia3+mrowa
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
!
5004   ia1=1
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
5040   ib1=ib1+mcolb
  ib2=ib2+mcolb
  ib3=ib3+mcolb
  ib4=ib4+mcolb
5050   ib5=ib5+mcolb
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
!
6000 continue
if(nlink*(ncol+nrow)+nrow*ncol.le.ncache) then
  call mxmb4(a,mcola,mrowa,b,mcolb,mrowb, &
             r,mcolr,mrowr,ncol,nlink,nrow)
  return
end if
#endif
#if defined(HIB_UNIX_BLAS3)
if(mcolr.eq.1) then
  if(mcola.eq.1.and.mcolb.eq.1) then
    call dgemm('N','N',ncol,nrow,nlink,1.0d0,a,mrowa, &
                     b,max(nlink,mrowb),1.0d0,r,mrowr)
    return
  else if(mrowa.eq.1.and.mcolb.eq.1) then
    call dgemm('T','N',ncol,nrow,nlink,1.0d0,a,mcola, &
                     b,max(nlink,mrowb),1.0d0,r,mrowr)
    return
  else if(mcola.eq.1.and.mrowb.eq.1) then
    call dgemm('N','T',ncol,nrow,nlink,1.0d0,a,mrowa, &
                     b,mcolb,1.0d0,r,mrowr)
    return
  else if(mrowa.eq.1.and.mrowb.eq.1) then
    call dgemm('T','T',ncol,nrow,nlink,1.0d0,a,mcola, &
                     b,mcolb,1.0d0,r,mrowr)
    return
  end if
else if(mrowr.eq.1) then
  if(mcola.eq.1.and.mcolb.eq.1) then
    call dgemm('T','T',nrow,ncol,nlink,1.0d0,b,max(nlink,mrowb), &
                     a,mrowa,1.0d0,r,mcolr)
    return
  else if(mrowa.eq.1.and.mcolb.eq.1) then
    call dgemm('T','N',nrow,ncol,nlink,1.0d0,b,max(nlink,mrowb), &
                     a,mcola,1.0d0,r,mcolr)
    return
  else if(mcola.eq.1.and.mrowb.eq.1) then
    call dgemm('N','T',nrow,ncol,nlink,1.0d0,b,mcolb, &
                     a,mrowa,1.0d0,r,mcolr)
    return
  else if(mrowa.eq.1.and.mrowb.eq.1) then
    call dgemm('N','N',nrow,ncol,nlink,1.0d0,b,mcolb, &
                     a,mcola,1.0d0,r,mcolr)
    return
  end if
end if
#endif
#if !defined(HIB_CRAY) && !defined(HIB_UNIX_CONVEX)

mxb=mxmblk
!      nkb=mxmbln
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
call mxmb4(a(ia1),mcola,mrowa,b(ib1),mcolb,mrowb,r(ir1), &
             mcolr,mrowr, ni,nk,nj)
80 ie=ie+ni
70 je=je+nj
60 ke=ke+nk
return
6001 ke=0
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
call mxmb4(a(ia1),mcola,mrowa,b(ib1),mcolb,mrowb,r(ir1), &
             mcolr,mrowr, ni,nk,nj)
71 je=je+nj
81 ie=ie+ni
61 ke=ke+nk
return
end
!--------------------------------------------------------------------
subroutine mxmb4(a,mcola,mrowa,b,mcolb,mrowb, &
r,mcolr,mrowr,ncol,nlink,nrow)
!--------------------------------------------------------------------
implicit double precision (a-h,o-z)
dimension r(1),a(1),b(1)
!... r(ncol,nrow)=a(ncol,nlink)*b(nlink,nrow) matrix mult
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
100 ia4=ia4+mrowa
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
200 continue
if(nrest.eq.0) goto 300
goto (201,202,203,300,205),nrest
201   ib1=ibb
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
101   ia4=ia4+mrowa
  r(ir1)=r(ir1)+s11
  r(ir2)=r(ir2)+s21
  r(ir3)=r(ir3)+s31
  r(ir4)=r(ir4)+s41
  goto 300
!
202   ib1=ibb
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
102   ia4=ia4+mrowa
  r(ir1)=r(ir1)+s11
  r(ir2)=r(ir2)+s21
  r(ir3)=r(ir3)+s31
  r(ir4)=r(ir4)+s41
  r(ir1+mrowr)=r(ir1+mrowr)+s12
  r(ir2+mrowr)=r(ir2+mrowr)+s22
  r(ir3+mrowr)=r(ir3+mrowr)+s32
  r(ir4+mrowr)=r(ir4+mrowr)+s42
  goto 300
!
203   ib1=ibb
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
103   ia4=ia4+mrowa
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
!
205   ib1=ibb
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
105   ia4=ia4+mrowa
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
300 iaa=iaa+4*mcola
if(ncest.eq.0) return
goto (301,302,303,304,305),ncest
301   ibb=1
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
1001   ia1=ia1+mrowa
  r(ir1)=r(ir1)+s11
  r(ir1+mrowr)=r(ir1+mrowr)+s12
  ir1=ir1+2*mrowr
2001   continue
  if(mod(nrow,2).eq.0) return
  s11=0
  do 1011 k=1,nlink
  s11=a(iaa)*b(ibb)+s11
  ibb=ibb+mcolb
1011   iaa=iaa+mrowa
  r(ir1)=r(ir1)+s11
  return
!
302   ibb=1
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
1002   ia2=ia2+mrowa
  r(ir1)=r(ir1)+s11
  r(ir2)=r(ir2)+s21
  r(ir1+mrowr)=r(ir1+mrowr)+s12
  r(ir2+mrowr)=r(ir2+mrowr)+s22
  ir1=ir1+2*mrowr
  ir2=ir2+2*mrowr
2002   continue
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
1012   ibb=ibb+mcolb
  r(ir1)=r(ir1)+s11
  r(ir2)=r(ir2)+s21
  return
!
303   ibb=1
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
1003   ia3=ia3+mrowa
  r(ir1)=r(ir1)+s11
  r(ir2)=r(ir2)+s21
  r(ir3)=r(ir3)+s31
  r(ir1+mrowr)=r(ir1+mrowr)+s12
  r(ir2+mrowr)=r(ir2+mrowr)+s22
  r(ir3+mrowr)=r(ir3+mrowr)+s32
  ir1=ir1+2*mrowr
  ir2=ir2+2*mrowr
  ir3=ir3+2*mrowr
2003   continue
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
1013   ia3=ia3+mrowa
  r(ir1)=r(ir1)+s11
  r(ir2)=r(ir2)+s21
  r(ir3)=r(ir3)+s31
304   return
305   ibb=1
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
1005   ia5=ia5+mrowa
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
2005   continue
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
1015   continue
  r(ir1)=r(ir1)+s11
  r(ir2)=r(ir2)+s21
  r(ir3)=r(ir3)+s31
  r(ir4)=r(ir4)+s41
  r(ir5)=r(ir5)+s51
  return
end
!--------------------------------------------------------------------
subroutine mxman(a,mcola,mrowa,b,mcolb,mrowb, &
r,mcolr,mrowr,ncol,nlink,nrow)
!--------------------------------------------------------------------
implicit double precision (a-h,o-z)
dimension r(1),a(1),b(1)
!... r(ncol,nrow)=-a(ncol,nlink)*b(nlink,nrow) matrix mult
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
100 ia4=ia4+mrowa
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
200 continue
if(nrest.eq.0) goto 300
goto (201,202,203,300,205),nrest
201   ib1=ibb
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
101   ia4=ia4+mrowa
  r(ir1)=-s11
  r(ir2)=-s21
  r(ir3)=-s31
  r(ir4)=-s41
  goto 300
!
202   ib1=ibb
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
102   ia4=ia4+mrowa
  r(ir1)=-s11
  r(ir2)=-s21
  r(ir3)=-s31
  r(ir4)=-s41
  r(ir1+mrowr)=-s12
  r(ir2+mrowr)=-s22
  r(ir3+mrowr)=-s32
  r(ir4+mrowr)=-s42
  goto 300
!
203   ib1=ibb
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
103   ia4=ia4+mrowa
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
!
205   ib1=ibb
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
105   ia4=ia4+mrowa
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
300 iaa=iaa+4*mcola
if(ncest.eq.0) return
goto (301,302,303,304,305),ncest
301   ibb=1
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
1001   ia1=ia1+mrowa
  r(ir1)=-s11
  r(ir1+mrowr)=-s12
  ir1=ir1+2*mrowr
2001   continue
  if(mod(nrow,2).eq.0) return
  s11=0
  do 1011 k=1,nlink
  s11=a(iaa)*b(ibb)+s11
  ibb=ibb+mcolb
1011   iaa=iaa+mrowa
  r(ir1)=-s11
  return
!
302   ibb=1
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
1002   ia2=ia2+mrowa
  r(ir1)=-s11
  r(ir2)=-s21
  r(ir1+mrowr)=-s12
  r(ir2+mrowr)=-s22
  ir1=ir1+2*mrowr
  ir2=ir2+2*mrowr
2002   continue
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
1012   ibb=ibb+mcolb
  r(ir1)=-s11
  r(ir2)=-s21
  return
!
303   ibb=1
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
1003   ia3=ia3+mrowa
  r(ir1)=-s11
  r(ir2)=-s21
  r(ir3)=-s31
  r(ir1+mrowr)=-s12
  r(ir2+mrowr)=-s22
  r(ir3+mrowr)=-s32
  ir1=ir1+2*mrowr
  ir2=ir2+2*mrowr
  ir3=ir3+2*mrowr
2003   continue
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
1013   ia3=ia3+mrowa
  r(ir1)=-s11
  r(ir2)=-s21
  r(ir3)=-s31
304   return
305   ibb=1
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
1005   ia5=ia5+mrowa
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
2005   continue
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
1015   continue
  r(ir1)=-s11
  r(ir2)=-s21
  r(ir3)=-s31
  r(ir4)=-s41
  r(ir5)=-s51
  return
end
!--------------------------------------------------------------------
subroutine mxmbn(a,mcola,mrowa,b,mcolb,mrowb, &
r,mcolr,mrowr,ncol,nlink,nrow)
!--------------------------------------------------------------------
implicit double precision (a-h,o-z)
dimension r(1),a(1),b(1)
!... r(ncol,nrow)=r(ncol,nrow)-a(ncol,nlink)*b(nlink,nrow) matrix mult
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
100 ia4=ia4+mrowa
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
200 continue
if(nrest.eq.0) goto 300
goto (201,202,203,300,205),nrest
201   ib1=ibb
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
101   ia4=ia4+mrowa
  r(ir1)=r(ir1)-s11
  r(ir2)=r(ir2)-s21
  r(ir3)=r(ir3)-s31
  r(ir4)=r(ir4)-s41
  goto 300
!
202   ib1=ibb
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
102   ia4=ia4+mrowa
  r(ir1)=r(ir1)-s11
  r(ir2)=r(ir2)-s21
  r(ir3)=r(ir3)-s31
  r(ir4)=r(ir4)-s41
  r(ir1+mrowr)=r(ir1+mrowr)-s12
  r(ir2+mrowr)=r(ir2+mrowr)-s22
  r(ir3+mrowr)=r(ir3+mrowr)-s32
  r(ir4+mrowr)=r(ir4+mrowr)-s42
  goto 300
!
203   ib1=ibb
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
103   ia4=ia4+mrowa
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
!
205   ib1=ibb
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
105   ia4=ia4+mrowa
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
300 iaa=iaa+4*mcola
if(ncest.eq.0) return
goto (301,302,303,304,305),ncest
301   ibb=1
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
1001   ia1=ia1+mrowa
  r(ir1)=r(ir1)-s11
  r(ir1+mrowr)=r(ir1+mrowr)-s12
  ir1=ir1+2*mrowr
2001   continue
  if(mod(nrow,2).eq.0) return
  s11=0
  do 1011 k=1,nlink
  s11=a(iaa)*b(ibb)+s11
  ibb=ibb+mcolb
1011   iaa=iaa+mrowa
  r(ir1)=r(ir1)-s11
  return
!
302   ibb=1
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
1002   ia2=ia2+mrowa
  r(ir1)=r(ir1)-s11
  r(ir2)=r(ir2)-s21
  r(ir1+mrowr)=r(ir1+mrowr)-s12
  r(ir2+mrowr)=r(ir2+mrowr)-s22
  ir1=ir1+2*mrowr
  ir2=ir2+2*mrowr
2002   continue
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
1012   ibb=ibb+mcolb
  r(ir1)=r(ir1)-s11
  r(ir2)=r(ir2)-s21
  return
!
303   ibb=1
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
1003   ia3=ia3+mrowa
  r(ir1)=r(ir1)-s11
  r(ir2)=r(ir2)-s21
  r(ir3)=r(ir3)-s31
  r(ir1+mrowr)=r(ir1+mrowr)-s12
  r(ir2+mrowr)=r(ir2+mrowr)-s22
  r(ir3+mrowr)=r(ir3+mrowr)-s32
  ir1=ir1+2*mrowr
  ir2=ir2+2*mrowr
  ir3=ir3+2*mrowr
2003   continue
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
1013   ia3=ia3+mrowa
  r(ir1)=r(ir1)-s11
  r(ir2)=r(ir2)-s21
  r(ir3)=r(ir3)-s31
304   return
305   ibb=1
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
1005   ia5=ia5+mrowa
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
2005   continue
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
1015   continue
  r(ir1)=r(ir1)-s11
  r(ir2)=r(ir2)-s21
  r(ir3)=r(ir3)-s31
  r(ir4)=r(ir4)-s41
  r(ir5)=r(ir5)-s51
  return
end

!--------------------------------------------------------------------
subroutine mxva(a,mcola,mrowa,b,mcolb, &
                r,mcolr,ncol,nlink)
!--------------------------------------------------------------------
implicit double precision (a-h,o-z)
real(8), intent(in) :: a(:)
integer, intent(in) :: mcola
integer, intent(in) :: mrowa
real(8), intent(in) :: b(:)
integer, intent(in) :: mcolb
real(8), intent(out) :: r(:)
integer, intent(in) :: mcolr
integer, intent(in) :: ncol
integer, intent(in) :: nlink
!...  r(ncol)=a(ncol,nlink)*b(nlink) matrix mult
! mcola is spacing between adjacent column elements of a (column stride)
! mrowa is spacing between adjacent row elements of a (row stride)
! mcolb is spacing between adjacent elements of b
! mcolr is spacing between elements of product vector
! ncol is number of rows in a and number of elements in r
! nlink is number of columns in a and number of elements in b
if(ncol.gt.4) goto 50
  ia1=1
  ib1=1
  if(ncol.eq.1) then
  s1=0
  do 10 k=1,nlink
  s1=a(ia1)*b(ib1)+s1
  ia1=ia1+mrowa
10   ib1=ib1+mcolb
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
20   ib1=ib1+mcolb
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
30   ib1=ib1+mcolb
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
40   ib1=ib1+mcolb
  r(1)=s1
  r(1+mcolr)=s2
  r(1+2*mcolr)=s3
  r(1+3*mcolr)=s4
  return
end if
!
50 mcola4=4*mcola
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
100 ib1=ib1+mcolb
r(ir1)=s1
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
200 iaa4=iaa4+mcola4
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
101   ib1=ib1+mcolb
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
102   ib1=ib1+mcolb
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
103   ib1=ib1+mcolb
  r(ir1)=s1
  r(ir2)=s2
  r(ir3)=s3
end if
return
end

subroutine mxva_mat(a,mcola,mrowa,b,mcolb, &
                r,mcolr,ncol,nlink)
   use, intrinsic :: ISO_C_BINDING   ! for C_LOC and C_F_POINTER
   implicit none
   real(8), intent(in), target :: a(:, :)
   integer, intent(in) :: mcola
   integer, intent(in) :: mrowa
   real(8), intent(in), target :: b(:, :)
   integer, intent(in) :: mcolb
   real(8), intent(out), target :: r(:, :)
   integer, intent(in) :: mcolr
   integer, intent(in) :: ncol
   integer, intent(in) :: nlink

   real(8), pointer :: a_as_vec(:)
   real(8), pointer :: b_as_vec(:)
   real(8), pointer :: r_as_vec(:)

   call C_F_POINTER (C_LOC(a), a_as_vec, [mrowa*mcola])
   call C_F_POINTER (C_LOC(b), b_as_vec, [mcola*mcolb])
   call C_F_POINTER (C_LOC(r), r_as_vec, [mcola*mcolr])
   call mxva(a_as_vec, mcola, mrowa, b_as_vec, mcolb, r_as_vec, mcolr, ncol, nlink)
end subroutine mxva_mat


!--------------------------------------------------------------------
subroutine mxvb(a,mcola,mrowa,b,mcolb, &
                r,mcolr,ncol,nlink)
!--------------------------------------------------------------------
implicit double precision (a-h,o-z)
dimension r(1),a(1),b(1)
!...  r(ncol)=a(ncol,nlink)*b(nlink) matrix mult
if(ncol.gt.4) goto 50
  ia1=1
  ib1=1
  if(ncol.eq.1) then
  s1=0
  do 10 k=1,nlink
  s1=a(ia1)*b(ib1)+s1
  ia1=ia1+mrowa
10   ib1=ib1+mcolb
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
20   ib1=ib1+mcolb
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
30   ib1=ib1+mcolb
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
40   ib1=ib1+mcolb
  r(1)=r(1)+s1
  r(1+mcolr)=r(1+mcolr)+s2
  r(1+2*mcolr)=r(1+2*mcolr)+s3
  r(1+3*mcolr)=r(1+3*mcolr)+s4
  return
end if
!
50 mcola4=4*mcola
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
100 ib1=ib1+mcolb
r(ir1)=r(ir1)+s1
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
200 iaa4=iaa4+mcola4
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
101   ib1=ib1+mcolb
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
102   ib1=ib1+mcolb
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
103   ib1=ib1+mcolb
  r(ir1)=r(ir1)+s1
  r(ir2)=r(ir2)+s2
  r(ir3)=r(ir3)+s3
end if
return
end
#endif
! -------------------------
#if (defined(HIB_UNIX_AIX) || defined(HIB_UNIX_DEC) || defined(HIB_UNIX_HP) || defined(HIB_UNIX_IRIS) || defined(HIB_UNIX_SUN)) && !defined(HIB_UNIX_HP800)
subroutine dsymv ( uplo, n, alpha, a, lda, x, incx, &
                   beta, y, incy )
!     .. Scalar Arguments ..
double precision   alpha, beta
integer            incx, incy, lda, n
character*1        uplo
!     .. Array Arguments ..
double precision   a( lda, * ), x( * ), y( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYMV  performs the matrix-vector  operation
!
!     y := alpha*A*x + beta*y,
!
!  where alpha and beta are scalars, x and y are n element vectors and
!  A is an n by n symmetric matrix.
!
!  Parameters
!  ==========
!
!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:
!
!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.
!
!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.
!
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y. On exit, Y is overwritten by the updated
!           vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
double precision   one         , zero
parameter        ( one = 1.0d+0, zero = 0.0d+0 )
!     .. Local Scalars ..
double precision   temp1, temp2
integer            i, info, ix, iy, j, jx, jy, kx, ky
!     .. External Functions ..
logical            lsame
external           lsame
!     .. External Subroutines ..
external           xerbla
!     .. Intrinsic Functions ..
intrinsic          max
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
info = 0
if     ( .not.lsame( uplo, 'U' ).and. &
         .not.lsame( uplo, 'L' )      )then
   info = 1
else if( n.lt.0 )then
   info = 2
else if( lda.lt.max( 1, n ) )then
   info = 5
else if( incx.eq.0 )then
   info = 7
else if( incy.eq.0 )then
   info = 10
end if
if( info.ne.0 )then
   call xerbla( 'DSYMV ', info )
   return
end if
!
!     Quick return if possible.
!
if( ( n.eq.0 ).or.( ( alpha.eq.zero ).and.( beta.eq.one ) ) ) &
   return
!
!     Set up the start points in  X  and  Y.
!
if( incx.gt.0 )then
   kx = 1
else
   kx = 1 - ( n - 1 )*incx
end if
if( incy.gt.0 )then
   ky = 1
else
   ky = 1 - ( n - 1 )*incy
end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
!     First form  y := beta*y.
!
if( beta.ne.one )then
   if( incy.eq.1 )then
      if( beta.eq.zero )then
         do 10, i = 1, n
            y( i ) = zero
10          continue
      else
         do 20, i = 1, n
            y( i ) = beta*y( i )
20          continue
      end if
   else
      iy = ky
      if( beta.eq.zero )then
         do 30, i = 1, n
            y( iy ) = zero
            iy      = iy   + incy
30          continue
      else
         do 40, i = 1, n
            y( iy ) = beta*y( iy )
            iy      = iy           + incy
40          continue
      end if
   end if
end if
if( alpha.eq.zero ) &
   return
if( lsame( uplo, 'U' ) )then
!
!        Form  y  when A is stored in upper triangle.
!
   if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
      do 60, j = 1, n
         temp1 = alpha*x( j )
         temp2 = zero
         do 50, i = 1, j - 1
            y( i ) = y( i ) + temp1*a( i, j )
            temp2  = temp2  + a( i, j )*x( i )
50          continue
         y( j ) = y( j ) + temp1*a( j, j ) + alpha*temp2
60       continue
   else
      jx = kx
      jy = ky
      do 80, j = 1, n
         temp1 = alpha*x( jx )
         temp2 = zero
         ix    = kx
         iy    = ky
         do 70, i = 1, j - 1
            y( iy ) = y( iy ) + temp1*a( i, j )
            temp2   = temp2   + a( i, j )*x( ix )
            ix      = ix      + incx
            iy      = iy      + incy
70          continue
         y( jy ) = y( jy ) + temp1*a( j, j ) + alpha*temp2
         jx      = jx      + incx
         jy      = jy      + incy
80       continue
   end if
else
!
!        Form  y  when A is stored in lower triangle.
!
   if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
      do 100, j = 1, n
         temp1  = alpha*x( j )
         temp2  = zero
         y( j ) = y( j )       + temp1*a( j, j )
         do 90, i = j + 1, n
            y( i ) = y( i ) + temp1*a( i, j )
            temp2  = temp2  + a( i, j )*x( i )
90          continue
         y( j ) = y( j ) + alpha*temp2
100       continue
   else
      jx = kx
      jy = ky
      do 120, j = 1, n
         temp1   = alpha*x( jx )
         temp2   = zero
         y( jy ) = y( jy )       + temp1*a( j, j )
         ix      = jx
         iy      = jy
         do 110, i = j + 1, n
            ix      = ix      + incx
            iy      = iy      + incy
            y( iy ) = y( iy ) + temp1*a( i, j )
            temp2   = temp2   + a( i, j )*x( ix )
110          continue
         y( jy ) = y( jy ) + alpha*temp2
         jx      = jx      + incx
         jy      = jy      + incy
120       continue
   end if
end if
!
return
!
!     End of DSYMV .
!
end
subroutine dgemv ( trans, m, n, alpha, a, lda, x, incx, &
                   beta, y, incy )
!     .. Scalar Arguments ..
double precision   alpha, beta
integer            incx, incy, lda, m, n
character*1        trans
!     .. Array Arguments ..
double precision   a( lda, * ), x( * ), y( * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
double precision   one         , zero
parameter        ( one = 1.0d+0, zero = 0.0d+0 )
!     .. Local Scalars ..
double precision   temp
integer            i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
!     .. External Functions ..
logical            lsame
external           lsame
!     .. External Subroutines ..
external           xerbla
!     .. Intrinsic Functions ..
intrinsic          max
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
info = 0
if     ( .not.lsame( trans, 'N' ).and. &
         .not.lsame( trans, 'T' ).and. &
         .not.lsame( trans, 'C' )      )then
   info = 1
else if( m.lt.0 )then
   info = 2
else if( n.lt.0 )then
   info = 3
else if( lda.lt.max( 1, m ) )then
   info = 6
else if( incx.eq.0 )then
   info = 8
else if( incy.eq.0 )then
   info = 11
end if
if( info.ne.0 )then
   call xerbla( 'DGEMV ', info )
   return
end if
!
!     Quick return if possible.
!
if( ( m.eq.0 ).or.( n.eq.0 ).or. &
    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) ) &
   return
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
if( lsame( trans, 'N' ) )then
   lenx = n
   leny = m
else
   lenx = m
   leny = n
end if
if( incx.gt.0 )then
   kx = 1
else
   kx = 1 - ( lenx - 1 )*incx
end if
if( incy.gt.0 )then
   ky = 1
else
   ky = 1 - ( leny - 1 )*incy
end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
if( beta.ne.one )then
   if( incy.eq.1 )then
      if( beta.eq.zero )then
         do 10, i = 1, leny
            y( i ) = zero
10          continue
      else
         do 20, i = 1, leny
            y( i ) = beta*y( i )
20          continue
      end if
   else
      iy = ky
      if( beta.eq.zero )then
         do 30, i = 1, leny
            y( iy ) = zero
            iy      = iy   + incy
30          continue
      else
         do 40, i = 1, leny
            y( iy ) = beta*y( iy )
            iy      = iy           + incy
40          continue
      end if
   end if
end if
if( alpha.eq.zero ) &
   return
if( lsame( trans, 'N' ) )then
!
!        Form  y := alpha*A*x + y.
!
   jx = kx
   if( incy.eq.1 )then
      do 60, j = 1, n
         if( x( jx ).ne.zero )then
            temp = alpha*x( jx )
            do 50, i = 1, m
               y( i ) = y( i ) + temp*a( i, j )
50             continue
         end if
         jx = jx + incx
60       continue
   else
      do 80, j = 1, n
         if( x( jx ).ne.zero )then
            temp = alpha*x( jx )
            iy   = ky
            do 70, i = 1, m
               y( iy ) = y( iy ) + temp*a( i, j )
               iy      = iy      + incy
70             continue
         end if
         jx = jx + incx
80       continue
   end if
else
!
!        Form  y := alpha*A'*x + y.
!
   jy = ky
   if( incx.eq.1 )then
      do 100, j = 1, n
         temp = zero
         do 90, i = 1, m
            temp = temp + a( i, j )*x( i )
90          continue
         y( jy ) = y( jy ) + alpha*temp
         jy      = jy      + incy
100       continue
   else
      do 120, j = 1, n
         temp = zero
         ix   = kx
         do 110, i = 1, m
            temp = temp + a( i, j )*x( ix )
            ix   = ix   + incx
110          continue
         y( jy ) = y( jy ) + alpha*temp
         jy      = jy      + incy
120       continue
   end if
end if
!
return
!
!     End of DGEMV .
!
end
subroutine dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                   beta, c, ldc )
!     .. Scalar Arguments ..
character*1        transa, transb
integer            m, n, k, lda, ldb, ldc
double precision   alpha, beta
!     .. Array Arguments ..
double precision   a( lda, * ), b( ldb, * ), c( ldc, * )
!     ..
!
!  Purpose
!  =======
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C := alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!  ==========
!
!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 3 Blas routine.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!
!     .. External Functions ..
logical            lsame
external           lsame
!     .. External Subroutines ..
external           xerbla
!     .. Intrinsic Functions ..
intrinsic          max
!     .. Local Scalars ..
logical            nota, notb
integer            i, info, j, l, ncola, nrowa, nrowb
double precision   temp
!     .. Parameters ..
double precision   one         , zero
parameter        ( one = 1.0d+0, zero = 0.0d+0 )
!     ..
!     .. Executable Statements ..
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
nota  = lsame( transa, 'N' )
notb  = lsame( transb, 'N' )
if( nota )then
   nrowa = m
   ncola = k
else
   nrowa = k
   ncola = m
end if
if( notb )then
   nrowb = k
else
   nrowb = n
end if
!
!     Test the input parameters.
!
info = 0
if(      ( .not.nota                 ).and. &
         ( .not.lsame( transa, 'C' ) ).and. &
         ( .not.lsame( transa, 'T' ) )      )then
   info = 1
else if( ( .not.notb                 ).and. &
         ( .not.lsame( transb, 'C' ) ).and. &
         ( .not.lsame( transb, 'T' ) )      )then
   info = 2
else if( m  .lt.0               )then
   info = 3
else if( n  .lt.0               )then
   info = 4
else if( k  .lt.0               )then
   info = 5
else if( lda.lt.max( 1, nrowa ) )then
   info = 8
else if( ldb.lt.max( 1, nrowb ) )then
   info = 10
else if( ldc.lt.max( 1, m     ) )then
   info = 13
end if
if( info.ne.0 )then
   call xerbla( 'DGEMM ', info )
   return
end if
!
!     Quick return if possible.
!
if( ( m.eq.0 ).or.( n.eq.0 ).or. &
    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) ) &
   return
!
!     And if  alpha.eq.zero.
!
if( alpha.eq.zero )then
   if( beta.eq.zero )then
      do 20, j = 1, n
         do 10, i = 1, m
            c( i, j ) = zero
10          continue
20       continue
   else
      do 40, j = 1, n
         do 30, i = 1, m
            c( i, j ) = beta*c( i, j )
30          continue
40       continue
   end if
   return
end if
!
!     Start the operations.
!
if( notb )then
   if( nota )then
!
!           Form  C := alpha*A*B + beta*C.
!
      do 90, j = 1, n
         if( beta.eq.zero )then
            do 50, i = 1, m
               c( i, j ) = zero
50             continue
         else if( beta.ne.one )then
            do 60, i = 1, m
               c( i, j ) = beta*c( i, j )
60             continue
         end if
         do 80, l = 1, k
            if( b( l, j ).ne.zero )then
               temp = alpha*b( l, j )
               do 70, i = 1, m
                  c( i, j ) = c( i, j ) + temp*a( i, l )
70                continue
            end if
80          continue
90       continue
   else
!
!           Form  C := alpha*A'*B + beta*C
!
      do 120, j = 1, n
         do 110, i = 1, m
            temp = zero
            do 100, l = 1, k
               temp = temp + a( l, i )*b( l, j )
100             continue
            if( beta.eq.zero )then
               c( i, j ) = alpha*temp
            else
               c( i, j ) = alpha*temp + beta*c( i, j )
            end if
110          continue
120       continue
   end if
else
   if( nota )then
!
!           Form  C := alpha*A*B' + beta*C
!
      do 170, j = 1, n
         if( beta.eq.zero )then
            do 130, i = 1, m
               c( i, j ) = zero
130             continue
         else if( beta.ne.one )then
            do 140, i = 1, m
               c( i, j ) = beta*c( i, j )
140             continue
         end if
         do 160, l = 1, k
            if( b( j, l ).ne.zero )then
               temp = alpha*b( j, l )
               do 150, i = 1, m
                  c( i, j ) = c( i, j ) + temp*a( i, l )
150                continue
            end if
160          continue
170       continue
   else
!
!           Form  C := alpha*A'*B' + beta*C
!
      do 200, j = 1, n
         do 190, i = 1, m
            temp = zero
            do 180, l = 1, k
               temp = temp + a( l, i )*b( j, l )
180             continue
            if( beta.eq.zero )then
               c( i, j ) = alpha*temp
            else
               c( i, j ) = alpha*temp + beta*c( i, j )
            end if
190          continue
200       continue
   end if
end if
!
return
!
!     End of DGEMM .
!
end
#endif

subroutine transp (a, n, nmax)
!  subroutine to carry out in place transposition of n x n matrix a
!  of maximum row dimension nmax stored in packed column form
!  uses blas routine dswap
!  written by:  millard alexander
!  current revision date: 23-sept-87
implicit double precision (a-h,o-z)
integer icol, icolpt, irowpt, n, nmax, nrow
dimension a(1)
icolpt = 2
irowpt = nmax + 1
do 100 icol = 1, n - 1
!  icolpt points to first sub-diagonal element in column icol
!  irowpt points to first super-diagonal element in row icol
!  nrow is number of subdiagonal elements in this column
  nrow = n - icol
  call dswap (nrow, a(icolpt), 1, a(irowpt), nmax)
  icolpt = icolpt + nmax + 1
  irowpt = irowpt + nmax + 1
100 continue
return
end

end module mod_himatrix