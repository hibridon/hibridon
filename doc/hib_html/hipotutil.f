* -----------------------------------------
* spline routines from didier lemoine
      subroutine dspline(x,y,n,yp1,ypn,y2)
*  yp1 and ypn are the initial and final derivatives
*  y2 are the coeffcieints to be used later by dsplint
      implicit double precision(a-h,o-z)
      parameter (nmax=300)
      dimension x(n),y(n),y2(n),u(nmax)
      save
      if (n. gt. nmax) then
        write (6,5) n
        write (9,5) n
5       format (' *** N = ',i3, ' .GT. NMAX IN DSPLINE; ABORT')
        call exit
      endif
      if (yp1.gt..99d30) then
        y2(1)=0.d0
        u(1)=0.d0
      else
        y2(1)=-0.5d0
        u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.d0
        y2(i)=(sig-1.d0)/p
        u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99d30) then
        qn=0.d0
        un=0.d0
      else
        qn=0.5d0
        un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end
      subroutine dsplint(xa,ya,y2a,n,x,y)
      implicit double precision(a-h,o-z)
      dimension xa(n),ya(n),y2a(n)
      save
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
*      if (h.eq.0.d0) pause 'bad xa input.'
      if (h.eq.0.d0) then
        write (6, 15) h
        write (9, 15) h
15      format (' *** bad xa input in splint; abort')
        call exit
      endif
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
      end
* -----------------------------------------
* bspline routines from cmlib
cmlib:bspline        bint4
      subroutine bint4(x, y, ndata, ibcl, ibcr, fbcl, fbcr, kntopt,
     * t, bcoef, n, k, w)
c
c     written by d.e. amos, august, 1979.
*     simplification of error treatment by millard alexander dec 14 1989
c
c     references
c         sand78-1968
c
c         a practical guide to splines by c. de boor, applied
c         mathematics series 27, springer, 1979.
c
c         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
c
c     abstract
c         bint4 computes the b representation (t,bcoef,n,k) of a
c         cubic spline (k=4) which interpolates data (x(i),y(i)),
c         i=1,ndata.  parameters ibcl, ibcr, fbcl, fbcr allow the
c         specification of the spline first or second derivative at
c         both x(1) and x(ndata).  when this data is not specified
c         by the problem, it is common practice to use a natural
c         spline by setting second derivatives at x(1) and x(ndata)
c         to zero (ibcl=ibcr=2,fbcl=fbcr=0.0).  the spline is defined
c         on t(4).le.x.le.t(n+1) with (ordered) interior knots at x(i)
c         values where n=ndata+2.  the knots t(1), t(2), t(3) lie to
c         the left of t(4)=x(1) and the knots t(n+2), t(n+3), t(n+4)
c         lie to the right of t(n+1)=x(ndata) in increasing order.  if
c         no extrapolation outside (x(1),x(ndata)) is anticipated, the
c         knots t(1)=t(2)=t(3)=t(4)=x(1) and t(n+2)=t(n+3)=t(n+4)=
c         t(n+1)=x(ndata) can be specified by kntopt=1.  kntopt=2
c         selects a knot placement for t(1), t(2), t(3) to make the
c         first 7 knots symmetric about t(4)=x(1) and similarly for
c         t(n+2), t(n+3), t(n+4) about t(n+1)=x(ndata).  kntopt=3
c         allows the user to make his own selection, in increasing
c         order, for t(1), t(2), t(3) to the left of x(1) and t(n+2),
c         t(n+3), t(n+4) to the right of x(ndata) in the work array
c         w(1) through w(6).  in any case, the interpolation on
c         t(4).le.x.le.t(n+1) by using function bvalu is unique
c         for given boundary conditions.
c
c         bint4 calls bspvd, bnfac, bnslv, r1mach
c
c     description of arguments
c         input
c           x      - x vector of abscissae of length ndata, distinct
c                    and in increasing order
c           y      - y vector of ordinates of length ndata
c           ndata  - number of data points, ndata.ge.2
c           ibcl   - selection parameter for left boundary condition
c                    ibcl = 1 constrain the first derivative at
c                             x(1) to fbcl
c                         = 2 constrain the second derivative at
c                             x(1) to fbcl
c           ibcr   - selection parameter for right boundary condition
c                    ibcr = 1 constrain first derivative at
c                             x(ndata) to fbcr
c                    ibcr = 2 constrain second derivative at
c                             x(ndata) to fbcr
c           fbcl   - left boundary values governed by ibcl
c           fbcr   - right boundary values governed by ibcr
c           kntopt - knot selection parameter
c                    kntopt = 1 sets knot multiplicity at t(4) and
c                               t(n+1) to 4
c                           = 2 sets a symmetric placement of knots
c                               about t(4) and t(n+1)
c                           = 3 sets t(i)=w(i) and t(n+1+i)=w(3+i),i=1,3
c                               where w(i),i=1,6 is supplied by the user
c           w      - work array of dimension at least 5*(ndata+2)
c                    if kntopt=3, then w(1),w(2),w(3) are knot values to
c                    the left of x(1) and w(4),w(5),w(6) are knot
c                    values to the right of x(ndata) in increasing
c                    order to be supplied by the user
c
c         output
c           t      - knot array of length n+4
c           bcoef  - b spline coefficient array of length n
c           n      - number of coefficients, n=ndata+2
c           k      - order of spline, k=4
c
c     error conditions
c         improper  input is a fatal error
c         singular system of equations is a fatal error
c***end prologue
c
      integer i, ibcl, ibcr, iflag, ilb, ileft, it, iub, iw, iwp, j,
     * jw, k, kntopt, n, ndata, ndm, np, nwrow
      double precision bcoef,fbcl,fbcr,t, tol,txn,tx1,vnikx,w,wdtol,
     :  work,
     : x, xl,
     * y
c     real r1mach
      dimension x(1), y(1), t(1), bcoef(2), w(5,1), vnikx(4,4), work(15)
      wdtol = 1.e-9
c     r1mach(4)
      tol = sqrt(wdtol)
      if (ndata.lt.2) go to 200
      ndm = ndata - 1
      do 10 i=1,ndm
        if (x(i).ge.x(i+1)) go to 210
   10 continue
      if (ibcl.lt.1 .or. ibcl.gt.2) go to 220
      if (ibcr.lt.1 .or. ibcr.gt.2) go to 230
      if (kntopt.lt.1 .or. kntopt.gt.3) go to 240
      k = 4
      n = ndata + 2
      np = n + 1
      do 20 i=1,ndata
        t(i+3) = x(i)
   20 continue
      go to (30, 50, 90), kntopt
c     set up knot array with multiplicity 4 at x(1) and x(ndata)
   30 continue
      do 40 i=1,3
        t(4-i) = x(1)
        t(np+i) = x(ndata)
   40 continue
      go to 110
c     set up knot array with symmetric placement about end points
   50 continue
      if (ndata.gt.3) go to 70
      xl = (x(ndata)-x(1))/3.0e0
      do 60 i=1,3
        t(4-i) = t(5-i) - xl
        t(np+i) = t(np+i-1) + xl
   60 continue
      go to 110
   70 continue
      tx1 = x(1) + x(1)
      txn = x(ndata) + x(ndata)
      do 80 i=1,3
        t(4-i) = tx1 - x(i+1)
        t(np+i) = txn - x(ndata-i)
   80 continue
      go to 110
c     set up knot array less than x(1) and greater than x(ndata) to be
c     supplied by user in work locations w(1) through w(6) when kntopt=3
   90 continue
      do 100 i=1,3
        t(4-i) = w(4-i,1)
        jw = max0(1,i-1)
        iw = mod(i+2,5) + 1
        t(np+i) = w(iw,jw)
        if (t(4-i).gt.t(5-i)) go to 250
        if (t(np+i).lt.t(np+i-1)) go to 250
  100 continue
  110 continue
c
      do 130 i=1,5
        do 120 j=1,n
          w(i,j) = 0.0e0
  120   continue
  130 continue
c     set up left interpolation point and left boundary condition for
c     right limits
      it = ibcl + 1
      call bspvd(t, k, it, x(1), k, 4, vnikx, work)
      iw = 0
      if (abs(vnikx(3,1)).lt.tol) iw = 1
      do 140 j=1,3
        w(j+1,4-j) = vnikx(4-j,it)
        w(j,4-j) = vnikx(4-j,1)
  140 continue
      bcoef(1) = y(1)
      bcoef(2) = fbcl
c     set up interpolation equations for points i=2 to i=ndata-1
      ileft = 4
      if (ndm.lt.2) go to 170
      do 160 i=2,ndm
        ileft = ileft + 1
        call bspvd(t, k, 1, x(i), ileft, 4, vnikx, work)
        do 150 j=1,3
          w(j+1,3+i-j) = vnikx(4-j,1)
  150   continue
        bcoef(i+1) = y(i)
  160 continue
c     set up right interpolation point and right boundary condition for
c     left limits(ileft is associated with t(n)=x(ndata-1))
  170 continue
      it = ibcr + 1
      call bspvd(t, k, it, x(ndata), ileft, 4, vnikx, work)
      jw = 0
      if (abs(vnikx(2,1)).lt.tol) jw = 1
      do 180 j=1,3
        w(j+1,3+ndata-j) = vnikx(5-j,it)
        w(j+2,3+ndata-j) = vnikx(5-j,1)
  180 continue
      bcoef(n-1) = fbcr
      bcoef(n) = y(ndata)
c     solve system of equations
      ilb = 2 - jw
      iub = 2 - iw
      nwrow = 5
      iwp = iw + 1
      call bnfac(w(iwp,1), nwrow, n, ilb, iub, iflag)
      if (iflag.eq.2) go to 190
      call bnslv(w(iwp,1), nwrow, n, ilb, iub, bcoef)
      return
c
c
  190 continue
      write(6, 191)
      write(9, 191)
  191 format(44h bint4,  the system of equations is singular)
      call xerror
  200 continue
      write (6, 201)
      write (9, 201)
  201 format(29h bint4,  ndata is less than 2)
      call xerror
  210 continue
      write (6,211)
      write (9,211)
  211 format(
     * 49h bint4,  x values are not distinct or not ordered)
      call xerror
  220 continue
      write (6,221)
      write (9,221)
  221 format(27h bint4,  ibcl is not 1 or 2)
      call xerror
  230 continue
      write (6,231)
      write (9,231)
  231 format(27h bint4,  ibcr is not 1 or 2)
      call xerror
  240 continue
      write (6,241)
      write (9,241)
  241 format(33h bint4,  kntopt is not 1, 2, or 3)
      call xerror
  250 continue
      write (6,251)
      write (9,251)
  251 format('BINT4,  KNOT INPUT THROUGH W ARRAY IS NOT ORDERED',
     :       ' PROPERLY')
      call xerror
      end
c$$$cmlib:bspline        bnfac
      subroutine bnfac(w, nroww, nrow, nbandl, nbandu, iflag)
c
c  bnfac is the banfac routine from
c        * a practical guide to splines *  by c. de boor
c
c  returns in  w  the lu-factorization (without pivoting) of the banded
c  matrix  a  of order  nrow  with  (nbandl + 1 + nbandu) bands or diag-
c  onals in the work array  w .
c
c******  i n p u t  ******
c  w.....work array of size  (nroww,nrow)  containing the interesting
c        part of a banded matrix  a , with the diagonals or bands of  a
c        stored in the rows of  w , while columns of  a  correspond to
c        columns of  w . this is the storage mode used in  linpack  and
c        results in efficient innermost loops.
c           explicitly,  a  has  nbandl  bands below the diagonal
c                            +     1     (main) diagonal
c                            +   nbandu  bands above the diagonal
c        and thus, with    middle = nbandu + 1,
c          a(i+j,j)  is in  w(i+middle,j)  for i=-nbandu,...,nbandl
c                                              j=1,...,nrow .
c        for example, the interesting entries of a (1,2)-banded matrix
c        of order  9  would appear in the first  1+1+2 = 4  rows of  w
c        as follows.
c                          13 24 35 46 57 68 79
c                       12 23 34 45 56 67 78 89
c                    11 22 33 44 55 66 77 88 99
c                    21 32 43 54 65 76 87 98
c
c        all other entries of  w  not identified in this way with an en-
c        try of  a  are never referenced .
c  nroww.....row dimension of the work array  w .
c        must be  .ge.  nbandl + 1 + nbandu  .
c  nbandl.....number of bands of  a  below the main diagonal
c  nbandu.....number of bands of  a  above the main diagonal .
c
c******  o u t p u t  ******
c  iflag.....integer indicating success( = 1) or failure ( = 2) .
c     if  iflag = 1, then
c  w.....contains the lu-factorization of  a  into a unit lower triangu-
c        lar matrix  l  and an upper triangular matrix  u (both banded)
c        and stored in customary fashion over the corresponding entries
c        of  a . this makes it possible to solve any particular linear
c        system  a*x = b  for  x  by a
c              call bnslv ( w, nroww, nrow, nbandl, nbandu, b )
c        with the solution x  contained in  b  on return .
c     if  iflag = 2, then
c        one of  nrow-1, nbandl,nbandu failed to be nonnegative, or else
c        one of the potential pivots was found to be zero indicating
c        that  a  does not have an lu-factorization. this implies that
c        a  is singular in case it is totally positive .
c
c******  m e t h o d  ******
c     gauss elimination  w i t h o u t  pivoting is used. the routine is
c  intended for use with matrices  a  which do not require row inter-
c  changes during factorization, especially for the  t o t a l l y
c  p o s i t i v e  matrices which occur in spline calculations.
c     the routine should not be used for an arbitrary banded matrix.
c
      integer iflag, nbandl, nbandu, nrow, nroww, i, ipk, j, jmax, k,
     * kmax, middle, midmk, nrowm1
      double precision w(nroww,nrow), factor, pivot
c
      iflag = 1
      middle = nbandu + 1
c                         w(middle,.) contains the main diagonal of  a .
      nrowm1 = nrow - 1
      if (nrowm1) 120, 110, 10
   10 if (nbandl.gt.0) go to 30
c                a is upper triangular. check that diagonal is nonzero .
      do 20 i=1,nrowm1
        if (w(middle,i).eq.0.0e0) go to 120
   20 continue
      go to 110
   30 if (nbandu.gt.0) go to 60
c              a is lower triangular. check that diagonal is nonzero and
c                 divide each column by its diagonal .
      do 50 i=1,nrowm1
        pivot = w(middle,i)
        if (pivot.eq.0.0e0) go to 120
        jmax = min0(nbandl,nrow-i)
        do 40 j=1,jmax
          w(middle+j,i) = w(middle+j,i)/pivot
   40   continue
   50 continue
      return
c
c        a  is not just a triangular matrix. construct lu factorization
   60 do 100 i=1,nrowm1
c                                  w(middle,i)  is pivot for i-th step .
        pivot = w(middle,i)
        if (pivot.eq.0.0e0) go to 120
c                 jmax  is the number of (nonzero) entries in column  i
c                     below the diagonal .
        jmax = min0(nbandl,nrow-i)
c              divide each entry in column  i  below diagonal by pivot .
        do 70 j=1,jmax
          w(middle+j,i) = w(middle+j,i)/pivot
   70   continue
c                 kmax  is the number of (nonzero) entries in row  i  to
c                     the right of the diagonal .
        kmax = min0(nbandu,nrow-i)
c                  subtract  a(i,i+k)*(i-th column) from (i+k)-th column
c                  (below row  i ) .
        do 90 k=1,kmax
          ipk = i + k
          midmk = middle - k
          factor = w(midmk,ipk)
          do 80 j=1,jmax
            w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
   80     continue
   90   continue
  100 continue
c                                       check the last diagonal entry .
  110 if (w(middle,nrow).ne.0.0e0) return
  120 iflag = 2
      return
      end
c$$$cmlib:bspline        bnslv
      subroutine bnslv(w, nroww, nrow, nbandl, nbandu, b)
c
c  bnslv is the banslv routine from
c        * a practical guide to splines *  by c. de boor
c
c  companion routine to  bnfac . it returns the solution  x  of the
c  linear system  a*x = b  in place of  b , given the lu-factorization
c  for  a  in the work array  w from bnfac.
c
c******  i n p u t  ******
c  w, nroww,nrow,nbandl,nbandu.....describe the lu-factorization of a
c        banded matrix  a  of order  nrow  as constructed in  bnfac .
c        for details, see  bnfac .
c  b.....right side of the system to be solved .
c
c******  o u t p u t  ******
c  b.....contains the solution  x , of order  nrow .
c
c******  m e t h o d  ******
c     (with  a = l*u, as stored in  w,) the unit lower triangular system
c  l(u*x) = b  is solved for  y = u*x, and  y  stored in  b . then the
c  upper triangular system  u*x = y  is solved for  x  . the calcul-
c  ations are so arranged that the innermost loops stay within columns.
c
      integer nbandl, nbandu, nrow, nroww, i, j, jmax, middle, nrowm1
      double precision w(nroww,nrow), b(nrow)
      middle = nbandu + 1
      if (nrow.eq.1) go to 80
      nrowm1 = nrow - 1
      if (nbandl.eq.0) go to 30
c                                 forward pass
c            for i=1,2,...,nrow-1, subtract  right side(i)*(i-th column
c            of  l )  from right side  (below i-th row) .
      do 20 i=1,nrowm1
        jmax = min0(nbandl,nrow-i)
        do 10 j=1,jmax
          b(i+j) = b(i+j) - b(i)*w(middle+j,i)
   10   continue
   20 continue
c                                 backward pass
c            for i=nrow,nrow-1,...,1, divide right side(i) by i-th diag-
c            onal entry of  u, then subtract  right side(i)*(i-th column
c            of  u)  from right side  (above i-th row).
   30 if (nbandu.gt.0) go to 50
c                                a  is lower triangular .
      do 40 i=1,nrow
        b(i) = b(i)/w(1,i)
   40 continue
      return
   50 i = nrow
   60 b(i) = b(i)/w(middle,i)
      jmax = min0(nbandu,i-1)
      do 70 j=1,jmax
        b(i-j) = b(i-j) - b(i)*w(middle-j,i)
   70 continue
      i = i - 1
      if (i.gt.1) go to 60
   80 b(1) = b(1)/w(middle,1)
      return
      end
c$$$cmlib:bspline        bspvd
      subroutine bspvd(t, k, nderiv, x, ileft, ldvnik, vnikx, work)
c
c     written by carl de boor and modified by d. e. amos
c
c     reference
c         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
c
c     abstract
c         bspvd is the bsplvd routine of the reference.
c
c         bspvd calculates the value and all derivatives of order
c         less than nderiv of all basis functions which do not
c         (possibly) vanish at x.  ileft is input such that
c         t(ileft) .le. x .lt. t(ileft+1).  a call to intrv(t,n+1,x,
c         ilo,ileft,mflag) will produce the proper ileft.  the output of
c         bspvd is a matrix vnikx(i,j) of dimension at least (k,nderiv)
c         whose columns contain the k nonzero basis functions and
c         their nderiv-1 right derivatives at x, i=1,k, j=1,nderiv.
c         these basis functions have indices ileft-k+i, i=1,k,
c         k.le.ileft.le.n. the nonzero part of the i-th basis function
c         lies in (t(i),t(i+k)), i=1,n.
c
c         if x=t(ileft+1) then vnikx contains left limiting values
c         (left derivatives) at t(ileft+1).  in particular, ileft = n
c         produces left limiting values at the right end point
c         x=t(n+1). to obtain left limiting values at t(i), i=k+1,n+1,
c         set x= next lower distinct knot, call intrv to get ileft,
c         set x=t(i), and then call bspvd.
c
c         bspvd calls bspvn
c
c     description of arguments
c         input
c          t       - knot vector of length n+k, where
c                    n=number of b-spline basis functions
c                    n = sum of knot multiplicities-k
c          k       - order of the b-spline, k.ge.1
c          nderiv  - number of derivatives=nderiv-1, 1.le.nderiv.le.k
c          x       - argument of basis functions, t(k).le.x.le.t(n+1)
c          ileft   - largest integer such that
c                    t(ileft) .le. x .lt. t(ileft+1)
c          ldvnik  - leading dimension of matrix vnikx
c
c         output
c          vnikx   - matrix of dimension at least (k,nderiv) contain-
c                    ing the nonzero basis functions at x and their
c                    derivatives columnwise.
c          work    - a work vector of length (k+1)*(k+2)/2
c
c     error conditions
c         improper input is a fatal error
c***end prologue
c
      integer i,ideriv,ileft,ipkmd,j,jj,jlow,jm,jp1mid,k,kmd, kp1, l,
     * ldummy, m, mhigh, nderiv
      double precision factor, fkmd, t, v, vnikx, work, x
c     dimension t(ileft+k), work((k+1)*(k+2)/2)
c     a(i,j) = work(i+j*(j+1)/2),  i=1,j+1  j=1,k-1
c     a(i,k) = w0rk(i+k*(k-1)/2)  i=1.k
c     work(1) and work((k+1)*(k+2)/2) are not used.
      dimension t(1), vnikx(ldvnik,nderiv), work(1)
      if(k.lt.1) go to 200
      if(nderiv.lt.1 .or. nderiv.gt.k) go to 205
      if(ldvnik.lt.k) go to 210
      ideriv = nderiv
      kp1 = k + 1
      jj = kp1 - ideriv
      call bspvn(t, jj, k, 1, x, ileft, vnikx, work, iwork)
      if (ideriv.eq.1) go to 100
      mhigh = ideriv
      do 20 m=2,mhigh
        jp1mid = 1
        do 10 j=ideriv,k
          vnikx(j,ideriv) = vnikx(jp1mid,1)
          jp1mid = jp1mid + 1
   10   continue
        ideriv = ideriv - 1
        jj = kp1 - ideriv
        call bspvn(t, jj, k, 2, x, ileft, vnikx, work, iwork)
   20 continue
c
      jm = kp1*(kp1+1)/2
      do 30 l = 1,jm
        work(l) = 0.0e0
   30 continue
c     a(i,i) = work(i*(i+3)/2) = 1.0       i = 1,k
      l = 2
      j = 0
      do 40 i = 1,k
        j = j + l
        work(j) = 1.0e0
        l = l + 1
   40 continue
      kmd = k
      do 90 m=2,mhigh
        kmd = kmd - 1
        fkmd = float(kmd)
        i = ileft
        j = k
        jj = j*(j+1)/2
        jm = jj - j
        do 60 ldummy=1,kmd
          ipkmd = i + kmd
          factor = fkmd/(t(ipkmd)-t(i))
          do 50 l=1,j
            work(l+jj) = (work(l+jj)-work(l+jm))*factor
   50     continue
          i = i - 1
          j = j - 1
          jj = jm
          jm = jm - j
   60   continue
c
        do 80 i=1,k
          v = 0.0e0
          jlow = max0(i,m)
          jj = jlow*(jlow+1)/2
          do 70 j=jlow,k
            v = work(i+jj)*vnikx(j,m) + v
            jj = jj + j + 1
   70     continue
          vnikx(i,m) = v
   80   continue
   90 continue
  100 return
c
c
  200 continue
      write (6,201)
      write (9,201)
  201 format(34h bspvd,  k does not satisfy k.ge.1)
      call xerror
  205 continue
      write (6,206)
      write (9,206)
  206 format(49h bspvd,  nderiv does not satisfy 1.le.nderiv.le.k)
      return
  210 continue
      write (6,211)
      write (9,211)
  211 format(44h bspvd,  ldvnik does not satisfy ldvnik.ge.k)
      call xerror
      end
c$$$cmlib:bspline        bspvn
      subroutine bspvn(t, jhigh, k, index, x, ileft, vnikx, work, iwork)
c
c     written by carl de boor and modified by d. e. amos
c
c     reference
c         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
c
c     abstract
c         bspvn is the bsplvn routine of the reference.
c
c         bspvn calculates the value of all (possibly) nonzero basis
c         functions at x of order max(jhigh,(j+1)*(index-1)), where
c         t(k).le.x.le.t(n+1) and j=iwork is set inside the routine on
c         the first call when index=1.  ileft is such that t(ileft).le.
c         x.lt.t(ileft+1).  a call to intrv(t,n+1,x,ilo,ileft,mflag)
c         produces the proper ileft.  bspvn calculates using the basic
c         algorithm needed in bspvd.  if only basis functions are
c         desired, setting jhigh=k and index=1 can be faster than
c         calling bspvd, but extra coding is required for derivatives
c         (index=2) and bspvd is set up for this purpose.
c
c         left limiting values are set up as described in bspvd.
c
c     description of arguments
c         input
c          t       - knot vector of length n+k, where
c                    n=number of b-spline basis functions
c                    n = sum of knot multiplicities-k
c          jhigh   - order of b-spline, 1.le.jhigh.le.k
c          k       - highest possible order
c          index   - index = 1 gives basis functions of order jhigh
c                          = 2 denotes previous entry with work, iwork
c                              values saved for subsequent calls to
c                              bspvn.
c          x       - argument of basis functions, t(k).le.x.le.t(n+1)
c          ileft   - largest integer such that
c                    t(ileft) .le. x .lt. t(ileft+1)
c
c         output
c          vnikx   - vector of length k for spline values.
c          work    - a work vector of length 2*k
c          iwork   - a work parameter. both work and iwork contain
c                    information necessary to continue for index = 2.
c                    when index = 1 exclusively, these are scratch
c                    variables and can be used for other purposes.
c
c     error conditions
c         improper input is a fatal error.
c***end prologue
c
      integer ileft, imjp1, index, ipj, iwork, jhigh, jp1, jp1ml, k, l
      double precision t, vm, vmprev, vnikx, work, x
c     dimension t(ileft+jhigh)
      dimension t(1), vnikx(k), work(1)
c     content of j, deltam, deltap is expected unchanged between calls.
c     work(i) = deltap(i), work(k+i) = deltam(i), i = 1,k
      if(k.lt.1) go to 90
      if(jhigh.gt.k .or. jhigh.lt.1) go to 100
      if(index.lt.1 .or. index.gt.2) go to 105
      if(x.lt.t(ileft) .or. x.gt.t(ileft+1)) go to 110
      go to (10, 20), index
   10 iwork = 1
      vnikx(1) = 1.0e0
      if (iwork.ge.jhigh) go to 40
c
   20 ipj = ileft + iwork
      work(iwork) = t(ipj) - x
      imjp1 = ileft - iwork + 1
      work(k+iwork) = x - t(imjp1)
      vmprev = 0.0e0
      jp1 = iwork + 1
      do 30 l=1,iwork
        jp1ml = jp1 - l
        vm = vnikx(l)/(work(l)+work(k+jp1ml))
        vnikx(l) = vm*work(l) + vmprev
        vmprev = vm*work(k+jp1ml)
   30 continue
      vnikx(jp1) = vmprev
      iwork = jp1
      if (iwork.lt.jhigh) go to 20
c
   40 return
c
c
   90 continue
      write (6,91)
      write (9,91)
   91 format(34h bspvn,  k does not satisfy k.ge.1)
      call xerror
  100 continue
      write (6,101)
      write (9,101)
  101 format(47h bspvn,  jhigh does not satisfy 1.le.jhigh.le.k)
      call xerror
  105 continue
      write (6,106)
      write (9,106)
  106 format(28h bspvn,  index is not 1 or 2)
      call xerror
  110 continue
      write (6,111)
      write (9,111)
  111 format(' BSPVN,X DOES NOT SATISFY T(ILEFT).LE.X.LE.T(ILEFT+1)')
      call xerror
      end
c$$$cmlib:bspline        bvalu
      function bvalu(t, a, n, k, ideriv, x, inbv, work)
c
c     written by carl de boor and modified by d. e. amos
c
c     reference
c         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
c
c     abstract
c         bvalu is the bvalue function of the reference.
c
c         bvalu evaluates the b-representation (t,a,n,k) of a b-spline
c         at x for the function value on ideriv=0 or any of its
c         derivatives on ideriv=1,2,...,k-1.  right limiting values
c         (right derivatives) are returned except at the right end
c         point x=t(n+1) where left limiting values are computed.  the
c         spline is defined on t(k).le.x.le.t(n+1).  bvalu returns
c         a fatal error message when x is outside of this interval.
c
c         to compute left derivatives or left limiting values at a
c         knot t(i), replace n by i-1 and set x=t(i), i=k+1,n+1.
c
c         bvalu calls intrv
c
c     description of arguments
c         input
c          t       - knot vector of length n+k
c          a       - b-spline coefficient vector of length n
c          n       - number of b-spline coefficients
c                    n = sum of knot multiplicities-k
c          k       - order of the b-spline, k.ge.1
c          ideriv  - order of the derivative, 0.le.ideriv.le.k-1
c                    ideriv=0 returns the b-spline value
c          x       - argument, t(k).le.x.le.t(n+1)
c          inbv    - an initialization parameter which must be set
c                    to 1 the first time bvalu is called.
c
c         output
c          inbv    - inbv contains information for efficient process-
c                    ing after the initial call and inbv must not
c                    be changed by the user. distinct splines require
c                    distinct inbv parameters.
c          work    - work vector of length 3*k.
c          bvalu   - value of the ideriv-th derivative at x
c
c     error conditions
c         an improper input is a fatal error
c***end prologue
c
      integer i,ideriv,iderp1,ihi,ihmkmj,ilo,imk,imkpj, inbv, ipj,
     * ip1, ip1mj, j, jj, j1, j2, k, kmider, kmj, km1, kpk, mflag, n
      double precision a, fkmj, t, work, x, bvalu
c     dimension t(n+k), work(3*k)
      dimension t(1), a(n), work(1)
      bvalu = 0.0e0
      if(k.lt.1) go to 102
      if(n.lt.k) go to 101
      if(ideriv.lt.0 .or. ideriv.ge.k) go to 110
      kmider = k - ideriv
c
c *** find *i* in (k,n) such that t(i) .le. x .lt. t(i+1)
c     (or, .le. t(i+1) if t(i) .lt. t(i+1) = t(n+1)).
      km1 = k - 1
      call intrv(t, n+1, x, inbv, i, mflag)
      if (x.lt.t(k)) go to 120
      if (mflag.eq.0) go to 20
      if (x.gt.t(i)) go to 130
   10 if (i.eq.k) go to 140
      i = i - 1
      if (x.eq.t(i)) go to 10
c
c *** difference the coefficients *ideriv* times
c     work(i) = aj(i), work(k+i) = dp(i), work(k+k+i) = dm(i), i=1.k
c
   20 imk = i - k
      do 30 j=1,k
        imkpj = imk + j
        work(j) = a(imkpj)
   30 continue
      if (ideriv.eq.0) go to 60
      do 50 j=1,ideriv
        kmj = k - j
        fkmj = float(kmj)
        do 40 jj=1,kmj
          ihi = i + jj
          ihmkmj = ihi - kmj
          work(jj) = (work(jj+1)-work(jj))/(t(ihi)-t(ihmkmj))*fkmj
   40   continue
   50 continue
c
c *** compute value at *x* in (t(i),(t(i+1)) of ideriv-th derivative,
c     given its relevant b-spline coeff. in aj(1),...,aj(k-ideriv).
   60 if (ideriv.eq.km1) go to 100
      ip1 = i + 1
      kpk = k + k
      j1 = k + 1
      j2 = kpk + 1
      do 70 j=1,kmider
        ipj = i + j
        work(j1) = t(ipj) - x
        ip1mj = ip1 - j
        work(j2) = x - t(ip1mj)
        j1 = j1 + 1
        j2 = j2 + 1
   70 continue
      iderp1 = ideriv + 1
      do 90 j=iderp1,km1
        kmj = k - j
        ilo = kmj
        do 80 jj=1,kmj
          work(jj) = (work(jj+1)*work(kpk+ilo)+work(jj)
     1              *work(k+jj))/(work(kpk+ilo)+work(k+jj))
          ilo = ilo - 1
   80   continue
   90 continue
  100 bvalu = work(1)
      return
c
c
  101 continue
      write (6,99)
      write (9,99)
   99 format(34h bvalu,  n does not satisfy n.ge.k)
      call xerror
  102 continue
      write (6,103)
      write (9,103)
  103 format(34h bvalu,  k does not satisfy k.ge.1)
      call xerror
  110 continue
      write (6,111)
      write (9,111)
  111 format(49h bvalu,  ideriv does not satisfy 0.le.ideriv.lt.k)
      call xerror
  120 continue
      write (6,121)
      write (9,121)
  121 format(47h bvalu,  x is n0t greater than or equal to t(k))
      call xerror
  130 continue
      write (6,131)
      write (9,131)
  131 format(46h bvalu,  x is not less than or equal to t(n+1))
      call xerror
  140 continue
      write (6,141)
      write (9,141)
  141 format(
     * 57h bvalu,  a left limiting value cann0t be obtained at t(k))
      call xerror
      end
c$$$cmlib:bspline        intrv
      subroutine intrv(xt, lxt, x, ilo, ileft, mflag)
c
c     written by carl de boor and modified by d. e. amos
c
c     reference
c         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
c
c     abstract
c         intrv is the interv routine of the reference.
c
c         intrv computes the largest integer ileft in 1.le.ileft.le.
c         lxt such that xt(ileft).le.x where xt(*) is a subdivision of
c         the x interval.  precisely,
c
c                      x.lt.xt(1)                1         -1
c         if  xt(i).le.x.lt.xt(i+1)  then  ileft=i  , mflag=0
c           xt(lxt).le.x                         lxt        1,
c
c         that is, when multiplicities are present in the break point
c         to the left of x, the largest index is taken for ileft.
c
c     description of arguments
c         input
c          xt      - xt is a knot or break point vector of length lxt
c          lxt     - length of the xt vector
c          x       - argument
c          ilo     - an initialization parameter which must be set
c                    to 1 the first time the spline array xt is
c                    processed by intrv.
c
c         output
c          ilo     - ilo contains information for efficient process-
c                    ing after the initial call and ilo must not be
c                    changed by the user. distinct splines require
c                    distinct ilo parameters.
c          ileft   - largest integer satisfying xt(ileft).le.x
c          mflag   - signals when x lies out of bounds
c
c     error conditions
c         none
c***end prologue
c
      integer ihi, ileft, ilo, istep, lxt, mflag, middle
      double precision x, xt
      dimension xt(lxt)
      ihi = ilo + 1
      if (ihi.lt.lxt) go to 10
      if (x.ge.xt(lxt)) go to 110
      if (lxt.le.1) go to 90
      ilo = lxt - 1
      ihi = lxt
c
   10 if (x.ge.xt(ihi)) go to 40
      if (x.ge.xt(ilo)) go to 100
c
c *** now x .lt. xt(ihi) . find lower bound
      istep = 1
   20 ihi = ilo
      ilo = ihi - istep
      if (ilo.le.1) go to 30
      if (x.ge.xt(ilo)) go to 70
      istep = istep*2
      go to 20
   30 ilo = 1
      if (x.lt.xt(1)) go to 90
      go to 70
c *** now x .ge. xt(ilo) . find upper bound
   40 istep = 1
   50 ilo = ihi
      ihi = ilo + istep
      if (ihi.ge.lxt) go to 60
      if (x.lt.xt(ihi)) go to 70
      istep = istep*2
      go to 50
   60 if (x.ge.xt(lxt)) go to 110
      ihi = lxt
c
c *** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval
   70 middle = (ilo+ihi)/2
      if (middle.eq.ilo) go to 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1
      if (x.lt.xt(middle)) go to 80
      ilo = middle
      go to 70
   80 ihi = middle
      go to 70
c *** set output and return
   90 mflag = -1
      ileft = 1
      return
  100 mflag = 0
      ileft = ilo
      return
  110 mflag = 1
      ileft = lxt
      return
      end
      subroutine xerror
* simplified error handler
* written by millard alexander 15/12/89
      data numb / 0 /
      save numb
      numb=numb+1
      if (numb.ge.10) then
        write (6, 100)
        write (9, 100)
100     format(' number of spline error .ge.10; abort')
        stop
      endif
      end
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
c$$$cmlib:linpackd       dqrdc
      subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job)
      integer ldx,n,p,job
      integer jpvt(1)
      double precision x(ldx,1),qraux(1),work(1)
c
c     dqrdc uses householder transformations to compute the qr
c     factorization of an n by p matrix x.  column pivoting
c     based on the 2-norms of the reduced columns may be
c     performed at the users option.
c
c     on entry
c
c        x       double precision(ldx,p), where ldx .ge. n.
c                x contains the matrix whose decomposition is to be
c                computed.
c
c        ldx     integer.
c                ldx is the leading dimension of the array x.
c
c        n       integer.
c                n is the number of rows of the matrix x.
c
c        p       integer.
c                p is the number of columns of the matrix x.
c
c        jpvt    integer(p).
c                jpvt contains integers that control the selection
c                of the pivot columns.  the k-th column x(k) of x
c                is placed in one of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      column.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free column.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final column.
c
c                before the decomposition is computed, initial columns
c                are moved to the beginning of the array x and final
c                columns to the end.  both initial and final columns
c                are frozen in place during the computation and only
c                free columns are moved.  at the k-th stage of the
c                reduction, if x(k) is occupied by a free column
c                it is interchanged with the free column of largest
c                reduced norm.  jpvt is not referenced if
c                job .eq. 0.
c
c        work    double precision(p).
c                work is a work array.  work is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c     on return
c
c        x       x contains in its upper triangle the upper
c                triangular matrix r of the qr factorization.
c                below its diagonal x contains information from
c                which the orthogonal part of the decomposition
c                can be recovered.  note that if pivoting has
c                been requested, the decomposition is not that
c                of the original matrix x but that of x
c                with its columns permuted as described by jpvt.
c
c        qraux   double precision(p).
c                qraux contains further information required to recover
c                the orthogonal part of the decomposition.
c
c        jpvt    jpvt(k) contains the index of the column of the
c                original matrix that has been interchanged into
c                the k-th column, if pivoting was requested.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dqrdc uses the following functions and subprograms.
c
c     blas daxpy,ddot,dscal,dswap,dnrm2
c     fortran dabs,dmax1,min0,dsqrt
c
c     internal variables
c
      integer j,jp,l,lp1,lup,maxj,pl,pu
      double precision maxnrm,dnrm2,tt
      double precision ddot,nrmxl,t
      logical negj,swapj
c
c
      pl = 1
      pu = 0
      if (job .eq. 0) go to 60
c
c        pivoting has been requested.  rearrange the columns
c        according to jpvt.
c
         do 20 j = 1, p
            swapj = jpvt(j) .gt. 0
            negj = jpvt(j) .lt. 0
            jpvt(j) = j
            if (negj) jpvt(j) = -j
            if (.not.swapj) go to 10
               if (j .ne. pl) call dswap(n,x(1,pl),1,x(1,j),1)
               jpvt(j) = jpvt(pl)
               jpvt(pl) = j
               pl = pl + 1
   10       continue
   20    continue
         pu = p
         do 50 jj = 1, p
            j = p - jj + 1
            if (jpvt(j) .ge. 0) go to 40
               jpvt(j) = -jpvt(j)
               if (j .eq. pu) go to 30
                  call dswap(n,x(1,pu),1,x(1,j),1)
                  jp = jpvt(pu)
                  jpvt(pu) = jpvt(j)
                  jpvt(j) = jp
   30          continue
               pu = pu - 1
   40       continue
   50    continue
   60 continue
c
c     compute the norms of the free columns.
c
      if (pu .lt. pl) go to 80
      do 70 j = pl, pu
         qraux(j) = dnrm2(n,x(1,j),1)
         work(j) = qraux(j)
   70 continue
   80 continue
c
c     perform the householder reduction of x.
c
      lup = min0(n,p)
      do 200 l = 1, lup
         if (l .lt. pl .or. l .ge. pu) go to 120
c
c           locate the column of largest norm and bring it
c           into the pivot position.
c
            maxnrm = 0.0d0
            maxj = l
            do 100 j = l, pu
               if (qraux(j) .le. maxnrm) go to 90
                  maxnrm = qraux(j)
                  maxj = j
   90          continue
  100       continue
            if (maxj .eq. l) go to 110
               call dswap(n,x(1,l),1,x(1,maxj),1)
               qraux(maxj) = qraux(l)
               work(maxj) = work(l)
               jp = jpvt(maxj)
               jpvt(maxj) = jpvt(l)
               jpvt(l) = jp
  110       continue
  120    continue
         qraux(l) = 0.0d0
         if (l .eq. n) go to 190
c
c           compute the householder transformation for column l.
c
            nrmxl = dnrm2(n-l+1,x(l,l),1)
            if (nrmxl .eq. 0.0d0) go to 180
               if (x(l,l) .ne. 0.0d0) nrmxl = dsign(nrmxl,x(l,l))
               call dscal(n-l+1,1.0d0/nrmxl,x(l,l),1)
               x(l,l) = 1.0d0 + x(l,l)
c
c              apply the transformation to the remaining columns,
c              updating the norms.
c
               lp1 = l + 1
               if (p .lt. lp1) go to 170
               do 160 j = lp1, p
                  t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
                  call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
                  if (j .lt. pl .or. j .gt. pu) go to 150
                  if (qraux(j) .eq. 0.0d0) go to 150
                     tt = 1.0d0 - (dabs(x(l,j))/qraux(j))**2
                     tt = dmax1(tt,0.0d0)
                     t = tt
                     tt = 1.0d0 + 0.05d0*tt*(qraux(j)/work(j))**2
                     if (tt .eq. 1.0d0) go to 130
                        qraux(j) = qraux(j)*dsqrt(t)
                     go to 140
  130                continue
                        qraux(j) = dnrm2(n-l,x(l+1,j),1)
                        work(j) = qraux(j)
  140                continue
  150             continue
  160          continue
  170          continue
c
c              save the transformation.
c
               qraux(l) = x(l,l)
               x(l,l) = -nrmxl
  180       continue
  190    continue
  200 continue
      return
      end
c$$$cmlib:linpackd       dqrsl
      subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
      integer ldx,n,k,job,info
      double precision x(ldx,1),qraux(1),y(1),qy(1),qty(1),b(1),rsd(1),
     *                 xb(1)
c
c     dqrsl applies the output of dqrdc to compute coordinate
c     transformations, projections, and least squares solutions.
c     for k .le. min(n,p), let xk be the matrix
c
c            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))
c
c     formed from columnns jpvt(1), ... ,jpvt(k) of the original
c     n x p matrix x that was input to dqrdc (if no pivoting was
c     done, xk consists of the first k columns of x in their
c     original order).  dqrdc produces a factored orthogonal matrix q
c     and an upper triangular matrix r such that
c
c              xk = q * (r)
c                       (0)
c
c     this information is contained in coded form in the arrays
c     x and qraux.
c
c     on entry
c
c        x      double precision(ldx,p).
c               x contains the output of dqrdc.
c
c        ldx    integer.
c               ldx is the leading dimension of the array x.
c
c        n      integer.
c               n is the number of rows of the matrix xk.  it must
c               have the same value as n in dqrdc.
c
c        k      integer.
c               k is the number of columns of the matrix xk.  k
c               must nnot be greater than min(n,p), where p is the
c               same as in the calling sequence to dqrdc.
c
c        qraux  double precision(p).
c               qraux contains the auxiliary output from dqrdc.
c
c        y      double precision(n)
c               y contains an n-vector that is to be manipulated
c               by dqrsl.
c
c        job    integer.
c               job specifies what is to be computed.  job has
c               the decimal expansion abcde, with the following
c               meaning.
c
c                    if a.ne.0, compute qy.
c                    if b,c,d, or e .ne. 0, compute qty.
c                    if c.ne.0, compute b.
c                    if d.ne.0, compute rsd.
c                    if e.ne.0, compute xb.
c
c               note that a request to compute b, rsd, or xb
c               automatically triggers the computation of qty, for
c               which an array must be provided in the calling
c               sequence.
c
c     on return
c
c        qy     double precision(n).
c               qy conntains q*y, if its computation has been
c               requested.
c
c        qty    double precision(n).
c               qty contains trans(q)*y, if its computation has
c               been requested.  here trans(q) is the
c               transpose of the matrix q.
c
c        b      double precision(k)
c               b contains the solution of the least squares problem
c
c                    minimize norm2(y - xk*b),
c
c               if its computation has been requested.  (note that
c               if pivoting was requested in dqrdc, the j-th
c               component of b will be associated with column jpvt(j)
c               of the original matrix x that was input into dqrdc.)
c
c        rsd    double precision(n).
c               rsd contains the least squares residual y - xk*b,
c               if its computation has been requested.  rsd is
c               also the orthogonal projection of y onto the
c               orthogonal complement of the column space of xk.
c
c        xb     double precision(n).
c               xb contains the least squares approximation xk*b,
c               if its computation has been requested.  xb is also
c               the orthogonal projection of y onto the column space
c               of x.
c
c        info   integer.
c               info is zero unless the computation of b has
c               been requested and r is exactly singular.  in
c               this case, info is the index of the first zero
c               diagonal element of r and b is left unaltered.
c
c     the parameters qy, qty, b, rsd, and xb are not referenced
c     if their computation is not requested and in this case
c     can be replaced by dummy variables in the calling program.
c     to save storage, the user may in some cases use the same
c     array for different parameters in the calling sequence.  a
c     frequently occuring example is when one wishes to compute
c     any of b, rsd, or xb and does not need y or qty.  in this
c     case one may identify y, qty, and one of b, rsd, or xb, while
c     providing separate arrays for anything else that is to be
c     computed.  thus the calling sequence
c
c          call dqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)
c
c     will result in the computation of b and rsd, with rsd
c     overwriting y.  more generally, each item in the following
c     list contains groups of permissible identifications for
c     a single callinng sequence.
c
c          1. (y,qty,b) (rsd) (xb) (qy)
c
c          2. (y,qty,rsd) (b) (xb) (qy)
c
c          3. (y,qty,xb) (b) (rsd) (qy)
c
c          4. (y,qy) (qty,b) (rsd) (xb)
c
c          5. (y,qy) (qty,rsd) (b) (xb)
c
c          6. (y,qy) (qty,xb) (b) (rsd)
c
c     in any group the value returned in the array allocated to
c     the group corresponds to the last member of the group.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dqrsl uses the following functions and subprograms.
c
c     blas daxpy,dcopy,ddot
c     fortran dabs,min0,mod
c
c     internal variables
c
      integer i,j,jj,ju,kp1
      double precision ddot,t,temp
      logical cb,cqy,cqty,cr,cxb
c
c
c     set info flag.
c
      info = 0
c
c     determine what is to be computed.
c
      cqy = job/10000 .ne. 0
      cqty = mod(job,10000) .ne. 0
      cb = mod(job,1000)/100 .ne. 0
      cr = mod(job,100)/10 .ne. 0
      cxb = mod(job,10) .ne. 0
      ju = min0(k,n-1)
c
c     special action when n=1.
c
      if (ju .ne. 0) go to 40
         if (cqy) qy(1) = y(1)
         if (cqty) qty(1) = y(1)
         if (cxb) xb(1) = y(1)
         if (.not.cb) go to 30
            if (x(1,1) .ne. 0.0d0) go to 10
               info = 1
            go to 20
   10       continue
               b(1) = y(1)/x(1,1)
   20       continue
   30    continue
         if (cr) rsd(1) = 0.0d0
      go to 250
   40 continue
c
c        set up to compute qy or qty.
c
         if (cqy) call dcopy(n,y,1,qy,1)
         if (cqty) call dcopy(n,y,1,qty,1)
         if (.not.cqy) go to 70
c
c           compute qy.
c
            do 60 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 50
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
                  call daxpy(n-j+1,t,x(j,j),1,qy(j),1)
                  x(j,j) = temp
   50          continue
   60       continue
   70    continue
         if (.not.cqty) go to 100
c
c           compute trans(q)*y.
c
            do 90 j = 1, ju
               if (qraux(j) .eq. 0.0d0) go to 80
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
                  call daxpy(n-j+1,t,x(j,j),1,qty(j),1)
                  x(j,j) = temp
   80          continue
   90       continue
  100    continue
c
c        set up to compute b, rsd, or xb.
c
         if (cb) call dcopy(k,qty,1,b,1)
         kp1 = k + 1
         if (cxb) call dcopy(k,qty,1,xb,1)
         if (cr .and. k .lt. n) call dcopy(n-k,qty(kp1),1,rsd(kp1),1)
         if (.not.cxb .or. kp1 .gt. n) go to 120
            do 110 i = kp1, n
               xb(i) = 0.0d0
  110       continue
  120    continue
         if (.not.cr) go to 140
            do 130 i = 1, k
               rsd(i) = 0.0d0
  130       continue
  140    continue
         if (.not.cb) go to 190
c
c           compute b.
c
            do 170 jj = 1, k
               j = k - jj + 1
               if (x(j,j) .ne. 0.0d0) go to 150
                  info = j
c           ......exit
                  go to 180
  150          continue
               b(j) = b(j)/x(j,j)
               if (j .eq. 1) go to 160
                  t = -b(j)
                  call daxpy(j-1,t,x(1,j),1,b,1)
  160          continue
  170       continue
  180       continue
  190    continue
         if (.not.cr .and. .not.cxb) go to 240
c
c           compute rsd or xb as required.
c
            do 230 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 220
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  if (.not.cr) go to 200
                     t = -ddot(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
                     call daxpy(n-j+1,t,x(j,j),1,rsd(j),1)
  200             continue
                  if (.not.cxb) go to 210
                     t = -ddot(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
                     call daxpy(n-j+1,t,x(j,j),1,xb(j),1)
  210             continue
                  x(j,j) = temp
  220          continue
  230       continue
  240    continue
  250 continue
      return
      end
