* -----------------------------------------
*  here follow necessary linear algebra routines from cmlib
* --------------------------------------------------
c$$$cmlib:dqrlss         dqrank
      subroutine dqrank(a,lda,m,n,tol,kr,jpvt,qraux,work)
c***begin prologue  dqrank
c***revision november 15, 1980
c***refer to  dqrlss
c***keyword(s)  overdetermined,least squares,linear equations
c***author  d. kahaner (nbs)
c***date written
c***purpose
c      computes the qr factorization of an  m  by  n  matrix a.  this
c      routine is used in conjunction with dqrlss to solve linear
c      systems of equations in a least square sense.
c***description
c
c     dqrank is used in conjunction with dqrlss to solve
c     overdetermined, underdetermined and singular linear systems
c     in a least squares sense.
c     dqrank uses the linpack subroutine dqrdc to compute the qr
c     factorization, with column pivoting, of an  m  by  n  matrix  a .
c     the numerical rank is determined using the tolerance tol.
c
c     for more information, see chapter 9 of linpack users guide,
c     j. dongarra et all, siam, 1979.
c
c     on entry
c
c        a     double precision (lda,n) .
c              the matrix whose decomposition is to be computed.
c
c        lda   integer.
c              the leading dimension of a .
c
c        m     integer.
c              the number of rows of a .
c
c        n     integer.
c              the number of columns of  a .
c
c        tol   double precision.
c              a relative tolerance used to determine the numerical
c              rank.  the problem should be scaled so that all the elements
c              of  a   have roughly the same absolute accuracy, eps.  then a
c              reasonable value for  tol  is roughly  eps  divided by
c              the magnitude of the largest element.
c
c        jpvt  integer(n)
c
c        qraux double precision(n)
c
c        work  double precision(n)
c              three auxilliary vectors used by dqrdc .
c
c     on return
c
c        a     contains the output from dqrdc.
c              the triangular matrix  r  of the qr factorization is
c              contained in the upper triangle and information needed
c              to recover the orthogonal matrix  q  is stored below
c              the diagonal in  a  and in the vector qraux .
c
c        kr    integer.
c              the numerical rank.
c
c        jpvt  the pivot information from dqrdc.
c
c     columns jpvt(1),...,jpvt(kr) of the original matrix are linearly
c     independent to within the tolerance tol and the remaining columns
c     are linearly dependent.  abs(a(1,1))/abs(a(kr,kr))  is an estimate
c     of the condition number of the matrix of independent columns,
c     and of  r .  this estimate will be .le. 1/tol .
c
c      usage.....see subroutine dqrlss
c
c***reference(s)
c      dongarra, et al, linpack users guide, siam, 1979
c***routines called  dqrdc
c***end prologue
      integer lda,m,n,kr,jpvt(1)
      double precision a(lda,1),tol,qraux(1),work(1)
      integer j,k
c***first executable statement
      do 10 j = 1, n
         jpvt(j) = 0
   10 continue
      call dqrdc(a,lda,m,n,qraux,jpvt,work,1)
      kr = 0
      k = min0(m,n)
      do 20 j = 1, k
         if (abs(a(j,j)) .le. tol*abs(a(1,1))) go to 30
         kr = j
   20 continue
   30 return
      end
c$$$cmlib:dqrlss         dqrlss
      subroutine dqrlss(a,lda,m,n,kr,b,x,rsd,jpvt,qraux)
c***begin prologue  dqrlss
c***revision november 15, 1980
c***category no.  d9
c***keyword(s)  overdetermined,least squares,linear equations
c***author  d. kahaner (nbs)
c***date written
c***purpose
c      solves an overdetermined, underdetermined, or singular system of
c      linear equations in least square sense.
c***description
c
c     dqrlss uses the linpack subroutine dqrsl to solve an overdetermined,
c     underdetermined, or singular linear system in a least squares
c     sense.  it must be preceeded by dqrank .
c     the system is  a*x  approximates  b  where  a  is  m  by  n  with
c     rank  kr  (determined by dqrank),  b  is a given  m-vector,
c     and  x  is the  n-vector to be computed.  a solution  x  with
c     at most  kr  nonzero components is found which minimimzes
c     the 2-norm of the residual,  a*x - b .
c
c     on entry
c
c        a,lda,m,n,kr,jpvt,qraux
c              the output from dqrank .
c
c        b     double precision(m) .
c              the right hand side of the linear system.
c
c     on return
c
c        x     double precision(n) .
c              a least squares solution to the linear system.
c
c        rsd   double precision(m) .
c              the residual, b - a*x .  rsd may overwite  b .
c
c      usage....
c        once the matrix a has been formed, dqrank should be
c      called once to decompose it.  then for each right hand
c      side, b, dqrlss should be called once to obtain the
c      solution and residual.
c
c
c***reference(s)
c      dongarra, et al, linpack users guide, siam, 1979
c***routines called  dqrsl
c***end prologue
      integer lda,m,n,kr,jpvt(1)
      integer j,k,info
      double precision t
      double precision a(lda,1),b(1),x(1),rsd(1),qraux(1)
c***first executable statement
      if (kr .ne. 0)
     1   call dqrsl(a,lda,m,kr,qraux,b,rsd,rsd,x,rsd,rsd,110,info)
      do 40 j = 1, n
         jpvt(j) = -jpvt(j)
         if (j .gt. kr) x(j) = 0.0
   40 continue
      do 70 j = 1, n
         if (jpvt(j) .gt. 0) go to 70
         k = -jpvt(j)
         jpvt(j) = k
   50    continue
            if (k .eq. j) go to 60
            t = x(j)
            x(j) = x(k)
            x(k) = t
            jpvt(k) = -jpvt(k)
            k = jpvt(k)
         go to 50
   60    continue
   70 continue
      return
      end
