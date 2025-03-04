! -----------------------------------------
! spline routines from didier lemoine
subroutine dspline(x,y,n,yp1,ypn,y2)
!  yp1 and ypn are the initial and final derivatives
!  y2 are the coeffcieints to be used later by dsplint

! NOTE abscissa should be in increasing order
implicit double precision(a-h,o-z)
dimension x(n),y(n),y2(n),u(n)
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
  u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11 continue
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
12 continue
return
end
subroutine dsplint(xa,ya,y2a,n,x,y)
implicit none
real(8), intent(in) :: xa(n)
real(8), intent(in) :: ya(n)
real(8), intent(in) :: y2a(n)
integer, intent(in) :: n
real(8), intent(in) :: x
real(8), intent(out) :: y
integer :: k, klo, khi
real(8) :: h, a, b
save
klo=1
khi=n
1 if (khi-klo.gt.1) then
  k=(khi+klo)/2
  if(xa(k).gt.x)then
    khi=k
  else
    klo=k
  endif
goto 1
endif
h=xa(khi)-xa(klo)
!      if (h.eq.0.d0) pause 'bad xa input.'
if (h.eq.0.d0) then
  write (6, 15) h
  write (9, 15) h
15   format (' *** bad xa input in splint; abort')
  call exit
endif
a=(xa(khi)-x)/h
b=(x-xa(klo))/h
y=a*ya(klo)+b*ya(khi)+ &
      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
return
end
! -----------------------------------------
! bspline routines from cmlib
!mlib:bspline        bint4
subroutine bint4(x, y, ndata, ibcl, ibcr, fbcl, fbcr, kntopt, &
 t, bcoef, n, k, w)
!
!     written by d.e. amos, august, 1979.
!     simplification of error treatment by millard alexander dec 14 1989
!
!     references
!         sand78-1968
!
!         a practical guide to splines by c. de boor, applied
!         mathematics series 27, springer, 1979.
!
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract
!         bint4 computes the b representation (t,bcoef,n,k) of a
!         cubic spline (k=4) which interpolates data (x(i),y(i)),
!         i=1,ndata.  parameters ibcl, ibcr, fbcl, fbcr allow the
!         specification of the spline first or second derivative at
!         both x(1) and x(ndata).  when this data is not specified
!         by the problem, it is common practice to use a natural
!         spline by setting second derivatives at x(1) and x(ndata)
!         to zero (ibcl=ibcr=2,fbcl=fbcr=0.0).  the spline is defined
!         on t(4).le.x.le.t(n+1) with (ordered) interior knots at x(i)
!         values where n=ndata+2.  the knots t(1), t(2), t(3) lie to
!         the left of t(4)=x(1) and the knots t(n+2), t(n+3), t(n+4)
!         lie to the right of t(n+1)=x(ndata) in increasing order.  if
!         no extrapolation outside (x(1),x(ndata)) is anticipated, the
!         knots t(1)=t(2)=t(3)=t(4)=x(1) and t(n+2)=t(n+3)=t(n+4)=
!         t(n+1)=x(ndata) can be specified by kntopt=1.  kntopt=2
!         selects a knot placement for t(1), t(2), t(3) to make the
!         first 7 knots symmetric about t(4)=x(1) and similarly for
!         t(n+2), t(n+3), t(n+4) about t(n+1)=x(ndata).  kntopt=3
!         allows the user to make his own selection, in increasing
!         order, for t(1), t(2), t(3) to the left of x(1) and t(n+2),
!         t(n+3), t(n+4) to the right of x(ndata) in the work array
!         w(1) through w(6).  in any case, the interpolation on
!         t(4).le.x.le.t(n+1) by using function bvalu is unique
!         for given boundary conditions.
!
!         bint4 calls bspvd, bnfac, bnslv, r1mach
!
!     description of arguments
!         input
!           x      - x vector of abscissae of length ndata, distinct
!                    and in increasing order
!           y      - y vector of ordinates of length ndata
!           ndata  - number of data points, ndata.ge.2
!           ibcl   - selection parameter for left boundary condition
!                    ibcl = 1 constrain the first derivative at
!                             x(1) to fbcl
!                         = 2 constrain the second derivative at
!                             x(1) to fbcl
!           ibcr   - selection parameter for right boundary condition
!                    ibcr = 1 constrain first derivative at
!                             x(ndata) to fbcr
!                    ibcr = 2 constrain second derivative at
!                             x(ndata) to fbcr
!           fbcl   - left boundary values governed by ibcl
!           fbcr   - right boundary values governed by ibcr
!           kntopt - knot selection parameter
!                    kntopt = 1 sets knot multiplicity at t(4) and
!                               t(n+1) to 4
!                           = 2 sets a symmetric placement of knots
!                               about t(4) and t(n+1)
!                           = 3 sets t(i)=w(i) and t(n+1+i)=w(3+i),i=1,3
!                               where w(i),i=1,6 is supplied by the user
!           w      - work array of dimension at least 5*(ndata+2)
!                    if kntopt=3, then w(1),w(2),w(3) are knot values to
!                    the left of x(1) and w(4),w(5),w(6) are knot
!                    values to the right of x(ndata) in increasing
!                    order to be supplied by the user
!
!         output
!           t      - knot array of length n+4
!           bcoef  - b spline coefficient array of length n
!           n      - number of coefficients, n=ndata+2
!           k      - order of spline, k=4
!
!     error conditions
!         improper  input is a fatal error
!         singular system of equations is a fatal error
!***end prologue
!
integer i, ibcl, ibcr, iflag, ilb, ileft, it, iub, iw, iwp, j, &
 jw, k, kntopt, n, ndata, ndm, np, nwrow
double precision bcoef,fbcl,fbcr,t, tol,txn,tx1,vnikx,w,wdtol, &
  work, &
 x, xl, &
 y
!     real r1mach
dimension x(ndata), y(ndata), t(n+4), bcoef(n), w(5,3), vnikx(4,4), work(15)
wdtol = 1.e-9
!     r1mach(4)
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
!     set up knot array with multiplicity 4 at x(1) and x(ndata)
30 continue
do 40 i=1,3
  t(4-i) = x(1)
  t(np+i) = x(ndata)
40 continue
go to 110
!     set up knot array with symmetric placement about end points
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
!     set up knot array less than x(1) and greater than x(ndata) to be
!     supplied by user in work locations w(1) through w(6) when kntopt=3
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
!
do 130 i=1,5
  do 120 j=1,n
    w(i,j) = 0.0e0
120   continue
130 continue
!     set up left interpolation point and left boundary condition for
!     right limits
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
!     set up interpolation equations for points i=2 to i=ndata-1
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
!     set up right interpolation point and right boundary condition for
!     left limits(ileft is associated with t(n)=x(ndata-1))
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
!     solve system of equations
ilb = 2 - jw
iub = 2 - iw
nwrow = 5
iwp = iw + 1
call bnfac(w(iwp,1), nwrow, n, ilb, iub, iflag)
if (iflag.eq.2) go to 190
call bnslv(w(iwp,1), nwrow, n, ilb, iub, bcoef)
return
!
!
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
211 format( &
 49h bint4,  x values are not distinct or not ordered)
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
251 format('BINT4,  KNOT INPUT THROUGH W ARRAY IS NOT ORDERED', &
       ' PROPERLY')
call xerror
end
!$$$cmlib:bspline        bnfac
subroutine bnfac(w, nroww, nrow, nbandl, nbandu, iflag)
!
!  bnfac is the banfac routine from
!        * a practical guide to splines *  by c. de boor
!
!  returns in  w  the lu-factorization (without pivoting) of the banded
!  matrix  a  of order  nrow  with  (nbandl + 1 + nbandu) bands or diag-
!  onals in the work array  w .
!
!******  i n p u t  ******
!  w.....work array of size  (nroww,nrow)  containing the interesting
!        part of a banded matrix  a , with the diagonals or bands of  a
!        stored in the rows of  w , while columns of  a  correspond to
!        columns of  w . this is the storage mode used in  linpack  and
!        results in efficient innermost loops.
!           explicitly,  a  has  nbandl  bands below the diagonal
!                            +     1     (main) diagonal
!                            +   nbandu  bands above the diagonal
!        and thus, with    middle = nbandu + 1,
!          a(i+j,j)  is in  w(i+middle,j)  for i=-nbandu,...,nbandl
!                                              j=1,...,nrow .
!        for example, the interesting entries of a (1,2)-banded matrix
!        of order  9  would appear in the first  1+1+2 = 4  rows of  w
!        as follows.
!                          13 24 35 46 57 68 79
!                       12 23 34 45 56 67 78 89
!                    11 22 33 44 55 66 77 88 99
!                    21 32 43 54 65 76 87 98
!
!        all other entries of  w  not identified in this way with an en-
!        try of  a  are never referenced .
!  nroww.....row dimension of the work array  w .
!        must be  .ge.  nbandl + 1 + nbandu  .
!  nbandl.....number of bands of  a  below the main diagonal
!  nbandu.....number of bands of  a  above the main diagonal .
!
!******  o u t p u t  ******
!  iflag.....integer indicating success( = 1) or failure ( = 2) .
!     if  iflag = 1, then
!  w.....contains the lu-factorization of  a  into a unit lower triangu-
!        lar matrix  l  and an upper triangular matrix  u (both banded)
!        and stored in customary fashion over the corresponding entries
!        of  a . this makes it possible to solve any particular linear
!        system  a*x = b  for  x  by a
!              call bnslv ( w, nroww, nrow, nbandl, nbandu, b )
!        with the solution x  contained in  b  on return .
!     if  iflag = 2, then
!        one of  nrow-1, nbandl,nbandu failed to be nonnegative, or else
!        one of the potential pivots was found to be zero indicating
!        that  a  does not have an lu-factorization. this implies that
!        a  is singular in case it is totally positive .
!
!******  m e t h o d  ******
!     gauss elimination  w i t h o u t  pivoting is used. the routine is
!  intended for use with matrices  a  which do not require row inter-
!  changes during factorization, especially for the  t o t a l l y
!  p o s i t i v e  matrices which occur in spline calculations.
!     the routine should not be used for an arbitrary banded matrix.
!
integer iflag, nbandl, nbandu, nrow, nroww, i, ipk, j, jmax, k, &
 kmax, middle, midmk, nrowm1
double precision w(nroww,nrow), factor, pivot
!
iflag = 1
middle = nbandu + 1
!                         w(middle,.) contains the main diagonal of  a .
nrowm1 = nrow - 1
if (nrowm1) 120, 110, 10
10 if (nbandl.gt.0) go to 30
!                a is upper triangular. check that diagonal is nonzero .
do 20 i=1,nrowm1
  if (w(middle,i).eq.0.0e0) go to 120
20 continue
go to 110
30 if (nbandu.gt.0) go to 60
!              a is lower triangular. check that diagonal is nonzero and
!                 divide each column by its diagonal .
do 50 i=1,nrowm1
  pivot = w(middle,i)
  if (pivot.eq.0.0e0) go to 120
  jmax = min0(nbandl,nrow-i)
  do 40 j=1,jmax
    w(middle+j,i) = w(middle+j,i)/pivot
40   continue
50 continue
return
!
!        a  is not just a triangular matrix. construct lu factorization
60 do 100 i=1,nrowm1
!                                  w(middle,i)  is pivot for i-th step .
  pivot = w(middle,i)
  if (pivot.eq.0.0e0) go to 120
!                 jmax  is the number of (nonzero) entries in column  i
!                     below the diagonal .
  jmax = min0(nbandl,nrow-i)
!              divide each entry in column  i  below diagonal by pivot .
  do 70 j=1,jmax
    w(middle+j,i) = w(middle+j,i)/pivot
70   continue
!                 kmax  is the number of (nonzero) entries in row  i  to
!                     the right of the diagonal .
  kmax = min0(nbandu,nrow-i)
!                  subtract  a(i,i+k)*(i-th column) from (i+k)-th column
!                  (below row  i ) .
  do 90 k=1,kmax
    ipk = i + k
    midmk = middle - k
    factor = w(midmk,ipk)
    do 80 j=1,jmax
      w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
80     continue
90   continue
100 continue
!                                       check the last diagonal entry .
110 if (w(middle,nrow).ne.0.0e0) return
120 iflag = 2
return
end
!$$$cmlib:bspline        bnslv
subroutine bnslv(w, nroww, nrow, nbandl, nbandu, b)
!
!  bnslv is the banslv routine from
!        * a practical guide to splines *  by c. de boor
!
!  companion routine to  bnfac . it returns the solution  x  of the
!  linear system  a*x = b  in place of  b , given the lu-factorization
!  for  a  in the work array  w from bnfac.
!
!******  i n p u t  ******
!  w, nroww,nrow,nbandl,nbandu.....describe the lu-factorization of a
!        banded matrix  a  of order  nrow  as constructed in  bnfac .
!        for details, see  bnfac .
!  b.....right side of the system to be solved .
!
!******  o u t p u t  ******
!  b.....contains the solution  x , of order  nrow .
!
!******  m e t h o d  ******
!     (with  a = l*u, as stored in  w,) the unit lower triangular system
!  l(u*x) = b  is solved for  y = u*x, and  y  stored in  b . then the
!  upper triangular system  u*x = y  is solved for  x  . the calcul-
!  ations are so arranged that the innermost loops stay within columns.
!
integer nbandl, nbandu, nrow, nroww, i, j, jmax, middle, nrowm1
double precision w(nroww,nrow), b(nrow)
middle = nbandu + 1
if (nrow.eq.1) go to 80
nrowm1 = nrow - 1
if (nbandl.eq.0) go to 30
!                                 forward pass
!            for i=1,2,...,nrow-1, subtract  right side(i)*(i-th column
!            of  l )  from right side  (below i-th row) .
do 20 i=1,nrowm1
  jmax = min0(nbandl,nrow-i)
  do 10 j=1,jmax
    b(i+j) = b(i+j) - b(i)*w(middle+j,i)
10   continue
20 continue
!                                 backward pass
!            for i=nrow,nrow-1,...,1, divide right side(i) by i-th diag-
!            onal entry of  u, then subtract  right side(i)*(i-th column
!            of  u)  from right side  (above i-th row).
30 if (nbandu.gt.0) go to 50
!                                a  is lower triangular .
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
!$$$cmlib:bspline        bspvd
subroutine bspvd(t, k, nderiv, x, ileft, ldvnik, vnikx, work)
!
!     written by carl de boor and modified by d. e. amos
!
!     reference
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract
!         bspvd is the bsplvd routine of the reference.
!
!         bspvd calculates the value and all derivatives of order
!         less than nderiv of all basis functions which do not
!         (possibly) vanish at x.  ileft is input such that
!         t(ileft) .le. x .lt. t(ileft+1).  a call to intrv(t,n+1,x,
!         ilo,ileft,mflag) will produce the proper ileft.  the output of
!         bspvd is a matrix vnikx(i,j) of dimension at least (k,nderiv)
!         whose columns contain the k nonzero basis functions and
!         their nderiv-1 right derivatives at x, i=1,k, j=1,nderiv.
!         these basis functions have indices ileft-k+i, i=1,k,
!         k.le.ileft.le.n. the nonzero part of the i-th basis function
!         lies in (t(i),t(i+k)), i=1,n.
!
!         if x=t(ileft+1) then vnikx contains left limiting values
!         (left derivatives) at t(ileft+1).  in particular, ileft = n
!         produces left limiting values at the right end point
!         x=t(n+1). to obtain left limiting values at t(i), i=k+1,n+1,
!         set x= next lower distinct knot, call intrv to get ileft,
!         set x=t(i), and then call bspvd.
!
!         bspvd calls bspvn
!
!     description of arguments
!         input
!          t       - knot vector of length n+k, where
!                    n=number of b-spline basis functions
!                    n = sum of knot multiplicities-k
!          k       - order of the b-spline, k.ge.1
!          nderiv  - number of derivatives=nderiv-1, 1.le.nderiv.le.k
!          x       - argument of basis functions, t(k).le.x.le.t(n+1)
!          ileft   - largest integer such that
!                    t(ileft) .le. x .lt. t(ileft+1)
!          ldvnik  - leading dimension of matrix vnikx
!
!         output
!          vnikx   - matrix of dimension at least (k,nderiv) contain-
!                    ing the nonzero basis functions at x and their
!                    derivatives columnwise.
!          work    - a work vector of length (k+1)*(k+2)/2
!
!     error conditions
!         improper input is a fatal error
!***end prologue
!
integer i,ideriv,ileft,ipkmd,j,jj,jlow,jm,jp1mid,k,kmd, kp1, l, &
 ldummy, m, mhigh, nderiv
double precision factor, fkmd, t, v, vnikx, work, x
!     dimension t(ileft+k), work((k+1)*(k+2)/2)
!     a(i,j) = work(i+j*(j+1)/2),  i=1,j+1  j=1,k-1
!     a(i,k) = w0rk(i+k*(k-1)/2)  i=1.k
!     work(1) and work((k+1)*(k+2)/2) are not used.
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
!
jm = kp1*(kp1+1)/2
do 30 l = 1,jm
  work(l) = 0.0e0
30 continue
!     a(i,i) = work(i*(i+3)/2) = 1.0       i = 1,k
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
!
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
!
!
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
!$$$cmlib:bspline        bspvn
subroutine bspvn(t, jhigh, k, index, x, ileft, vnikx, work, iwork)
!
!     written by carl de boor and modified by d. e. amos
!
!     reference
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract
!         bspvn is the bsplvn routine of the reference.
!
!         bspvn calculates the value of all (possibly) nonzero basis
!         functions at x of order max(jhigh,(j+1)*(index-1)), where
!         t(k).le.x.le.t(n+1) and j=iwork is set inside the routine on
!         the first call when index=1.  ileft is such that t(ileft).le.
!         x.lt.t(ileft+1).  a call to intrv(t,n+1,x,ilo,ileft,mflag)
!         produces the proper ileft.  bspvn calculates using the basic
!         algorithm needed in bspvd.  if only basis functions are
!         desired, setting jhigh=k and index=1 can be faster than
!         calling bspvd, but extra coding is required for derivatives
!         (index=2) and bspvd is set up for this purpose.
!
!         left limiting values are set up as described in bspvd.
!
!     description of arguments
!         input
!          t       - knot vector of length n+k, where
!                    n=number of b-spline basis functions
!                    n = sum of knot multiplicities-k
!          jhigh   - order of b-spline, 1.le.jhigh.le.k
!          k       - highest possible order
!          index   - index = 1 gives basis functions of order jhigh
!                          = 2 denotes previous entry with work, iwork
!                              values saved for subsequent calls to
!                              bspvn.
!          x       - argument of basis functions, t(k).le.x.le.t(n+1)
!          ileft   - largest integer such that
!                    t(ileft) .le. x .lt. t(ileft+1)
!
!         output
!          vnikx   - vector of length k for spline values.
!          work    - a work vector of length 2*k
!          iwork   - a work parameter. both work and iwork contain
!                    information necessary to continue for index = 2.
!                    when index = 1 exclusively, these are scratch
!                    variables and can be used for other purposes.
!
!     error conditions
!         improper input is a fatal error.
!***end prologue
!
integer ileft, imjp1, index, ipj, iwork, jhigh, jp1, jp1ml, k, l
double precision t, vm, vmprev, vnikx, work, x
!     dimension t(ileft+jhigh)
dimension t(1), vnikx(k), work(1)
!     content of j, deltam, deltap is expected unchanged between calls.
!     work(i) = deltap(i), work(k+i) = deltam(i), i = 1,k
if(k.lt.1) go to 90
if(jhigh.gt.k .or. jhigh.lt.1) go to 100
if(index.lt.1 .or. index.gt.2) go to 105
if(x.lt.t(ileft) .or. x.gt.t(ileft+1)) go to 110
go to (10, 20), index
10 iwork = 1
vnikx(1) = 1.0e0
if (iwork.ge.jhigh) go to 40
!
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
!
40 return
!
!
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
!$$$cmlib:bspline        bvalu
function bvalu(t, a, n, k, ideriv, x, inbv, work)
!
!     written by carl de boor and modified by d. e. amos
!
!     reference
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract
!         bvalu is the bvalue function of the reference.
!
!         bvalu evaluates the b-representation (t,a,n,k) of a b-spline
!         at x for the function value on ideriv=0 or any of its
!         derivatives on ideriv=1,2,...,k-1.  right limiting values
!         (right derivatives) are returned except at the right end
!         point x=t(n+1) where left limiting values are computed.  the
!         spline is defined on t(k).le.x.le.t(n+1).  bvalu returns
!         a fatal error message when x is outside of this interval.
!
!         to compute left derivatives or left limiting values at a
!         knot t(i), replace n by i-1 and set x=t(i), i=k+1,n+1.
!
!         bvalu calls intrv
!
!     description of arguments
!         input
!          t       - knot vector of length n+k
!          a       - b-spline coefficient vector of length n
!          n       - number of b-spline coefficients
!                    n = sum of knot multiplicities-k
!          k       - order of the b-spline, k.ge.1
!          ideriv  - order of the derivative, 0.le.ideriv.le.k-1
!                    ideriv=0 returns the b-spline value
!          x       - argument, t(k).le.x.le.t(n+1)
!          inbv    - an initialization parameter which must be set
!                    to 1 the first time bvalu is called.
!
!         output
!          inbv    - inbv contains information for efficient process-
!                    ing after the initial call and inbv must not
!                    be changed by the user. distinct splines require
!                    distinct inbv parameters.
!          work    - work vector of length 3*k.
!          bvalu   - value of the ideriv-th derivative at x
!
!     error conditions
!         an improper input is a fatal error
!***end prologue
!
integer i,ideriv,iderp1,ihi,ihmkmj,ilo,imk,imkpj, inbv, ipj, &
 ip1, ip1mj, j, jj, j1, j2, k, kmider, kmj, km1, kpk, mflag, n
double precision a, fkmj, t, work, x, bvalu
!     dimension t(n+k), work(3*k)
dimension t(1), a(n), work(3*k)
bvalu = 0.0e0
if(k.lt.1) go to 102
if(n.lt.k) go to 101
if(ideriv.lt.0 .or. ideriv.ge.k) go to 110
kmider = k - ideriv
!
! *** find *i* in (k,n) such that t(i) .le. x .lt. t(i+1)
!     (or, .le. t(i+1) if t(i) .lt. t(i+1) = t(n+1)).
km1 = k - 1
call intrv(t, n+1, x, inbv, i, mflag)
if (x.lt.t(k)) go to 120
if (mflag.eq.0) go to 20
if (x.gt.t(i)) go to 130
10 if (i.eq.k) go to 140
i = i - 1
if (x.eq.t(i)) go to 10
!
! *** difference the coefficients *ideriv* times
!     work(i) = aj(i), work(k+i) = dp(i), work(k+k+i) = dm(i), i=1.k
!
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
!
! *** compute value at *x* in (t(i),(t(i+1)) of ideriv-th derivative,
!     given its relevant b-spline coeff. in aj(1),...,aj(k-ideriv).
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
    work(jj) = (work(jj+1)*work(kpk+ilo)+work(jj) &
              *work(k+jj))/(work(kpk+ilo)+work(k+jj))
    ilo = ilo - 1
80   continue
90 continue
100 bvalu = work(1)
return
!
!
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
141 format( &
 57h bvalu,  a left limiting value cann0t be obtained at t(k))
call xerror
end
!$$$cmlib:bspline        intrv
subroutine intrv(xt, lxt, x, ilo, ileft, mflag)
!
!     written by carl de boor and modified by d. e. amos
!
!     reference
!         siam j. numerical analysis, 14, no. 3, june, 1977, pp.441-472.
!
!     abstract
!         intrv is the interv routine of the reference.
!
!         intrv computes the largest integer ileft in 1.le.ileft.le.
!         lxt such that xt(ileft).le.x where xt(*) is a subdivision of
!         the x interval.  precisely,
!
!                      x.lt.xt(1)                1         -1
!         if  xt(i).le.x.lt.xt(i+1)  then  ileft=i  , mflag=0
!           xt(lxt).le.x                         lxt        1,
!
!         that is, when multiplicities are present in the break point
!         to the left of x, the largest index is taken for ileft.
!
!     description of arguments
!         input
!          xt      - xt is a knot or break point vector of length lxt
!          lxt     - length of the xt vector
!          x       - argument
!          ilo     - an initialization parameter which must be set
!                    to 1 the first time the spline array xt is
!                    processed by intrv.
!
!         output
!          ilo     - ilo contains information for efficient process-
!                    ing after the initial call and ilo must not be
!                    changed by the user. distinct splines require
!                    distinct ilo parameters.
!          ileft   - largest integer satisfying xt(ileft).le.x
!          mflag   - signals when x lies out of bounds
!
!     error conditions
!         none
!***end prologue
!
integer ihi, ileft, ilo, istep, lxt, mflag, middle
double precision x, xt
dimension xt(lxt)
ihi = ilo + 1
if (ihi.lt.lxt) go to 10
if (x.ge.xt(lxt)) go to 110
if (lxt.le.1) go to 90
ilo = lxt - 1
ihi = lxt
!
10 if (x.ge.xt(ihi)) go to 40
if (x.ge.xt(ilo)) go to 100
!
! *** now x .lt. xt(ihi) . find lower bound
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
! *** now x .ge. xt(ilo) . find upper bound
40 istep = 1
50 ilo = ihi
ihi = ilo + istep
if (ihi.ge.lxt) go to 60
if (x.lt.xt(ihi)) go to 70
istep = istep*2
go to 50
60 if (x.ge.xt(lxt)) go to 110
ihi = lxt
!
! *** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval
70 middle = (ilo+ihi)/2
if (middle.eq.ilo) go to 100
!     note. it is assumed that middle = ilo in case ihi = ilo+1
if (x.lt.xt(middle)) go to 80
ilo = middle
go to 70
80 ihi = middle
go to 70
! *** set output and return
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
! simplified error handler
! written by millard alexander 15/12/89
data numb / 0 /
save numb
numb=numb+1
if (numb.ge.10) then
  write (6, 100)
  write (9, 100)
100   format(' number of spline error .ge.10; abort')
  stop
endif
end
! -----------------------------------------
!  here follow necessary linear algebra routines from cmlib
! --------------------------------------------------
!$$$cmlib:dqrlss         dqrank
subroutine dqrank(a,lda,m,n,tol,kr,jpvt,qraux,work)
!***begin prologue  dqrank
!***revision november 15, 1980
!***refer to  dqrlss
!***keyword(s)  overdetermined,least squares,linear equations
!***author  d. kahaner (nbs)
!***date written
!***purpose
!      computes the qr factorization of an  m  by  n  matrix a.  this
!      routine is used in conjunction with dqrlss to solve linear
!      systems of equations in a least square sense.
!***description
!
!     dqrank is used in conjunction with dqrlss to solve
!     overdetermined, underdetermined and singular linear systems
!     in a least squares sense.
!     dqrank uses the linpack subroutine dqrdc to compute the qr
!     factorization, with column pivoting, of an  m  by  n  matrix  a .
!     the numerical rank is determined using the tolerance tol.
!
!     for more information, see chapter 9 of linpack users guide,
!     j. dongarra et all, siam, 1979.
!
!     on entry
!
!        a     double precision (lda,n) .
!              the matrix whose decomposition is to be computed.
!
!        lda   integer.
!              the leading dimension of a .
!
!        m     integer.
!              the number of rows of a .
!
!        n     integer.
!              the number of columns of  a .
!
!        tol   double precision.
!              a relative tolerance used to determine the numerical
!              rank.  the problem should be scaled so that all the elements
!              of  a   have roughly the same absolute accuracy, eps.  then a
!              reasonable value for  tol  is roughly  eps  divided by
!              the magnitude of the largest element.
!
!        jpvt  integer(n)
!
!        qraux double precision(n)
!
!        work  double precision(n)
!              three auxilliary vectors used by dqrdc .
!
!     on return
!
!        a     contains the output from dqrdc.
!              the triangular matrix  r  of the qr factorization is
!              contained in the upper triangle and information needed
!              to recover the orthogonal matrix  q  is stored below
!              the diagonal in  a  and in the vector qraux .
!
!        kr    integer.
!              the numerical rank.
!
!        jpvt  the pivot information from dqrdc.
!
!     columns jpvt(1),...,jpvt(kr) of the original matrix are linearly
!     independent to within the tolerance tol and the remaining columns
!     are linearly dependent.  abs(a(1,1))/abs(a(kr,kr))  is an estimate
!     of the condition number of the matrix of independent columns,
!     and of  r .  this estimate will be .le. 1/tol .
!
!      usage.....see subroutine dqrlss
!
!***reference(s)
!      dongarra, et al, linpack users guide, siam, 1979
!***routines called  dqrdc
!***end prologue
integer lda,m,n,kr,jpvt(1)
double precision a(lda,1),tol,qraux(1),work(1)
integer j,k
!***first executable statement
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
!$$$cmlib:dqrlss         dqrlss
subroutine dqrlss(a,lda,m,n,kr,b,x,rsd,jpvt,qraux)
!***begin prologue  dqrlss
!***revision november 15, 1980
!***category no.  d9
!***keyword(s)  overdetermined,least squares,linear equations
!***author  d. kahaner (nbs)
!***date written
!***purpose
!      solves an overdetermined, underdetermined, or singular system of
!      linear equations in least square sense.
!***description
!
!     dqrlss uses the linpack subroutine dqrsl to solve an overdetermined,
!     underdetermined, or singular linear system in a least squares
!     sense.  it must be preceeded by dqrank .
!     the system is  a*x  approximates  b  where  a  is  m  by  n  with
!     rank  kr  (determined by dqrank),  b  is a given  m-vector,
!     and  x  is the  n-vector to be computed.  a solution  x  with
!     at most  kr  nonzero components is found which minimimzes
!     the 2-norm of the residual,  a*x - b .
!
!     on entry
!
!        a,lda,m,n,kr,jpvt,qraux
!              the output from dqrank .
!
!        b     double precision(m) .
!              the right hand side of the linear system.
!
!     on return
!
!        x     double precision(n) .
!              a least squares solution to the linear system.
!
!        rsd   double precision(m) .
!              the residual, b - a*x .  rsd may overwite  b .
!
!      usage....
!        once the matrix a has been formed, dqrank should be
!      called once to decompose it.  then for each right hand
!      side, b, dqrlss should be called once to obtain the
!      solution and residual.
!
!
!***reference(s)
!      dongarra, et al, linpack users guide, siam, 1979
!***routines called  dqrsl
!***end prologue
integer lda,m,n,kr,jpvt(1)
integer j,k,info
double precision t
double precision a(lda,1),b(1),x(1),rsd(1),qraux(1)
!***first executable statement
if (kr .ne. 0) &
   call dqrsl(a,lda,m,kr,qraux,b,rsd,rsd,x,rsd,rsd,110,info)
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
!$$$cmlib:linpackd       dqrdc
subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job)
use mod_hiutil, only: daxpy_wrapper
integer ldx,n,p,job
integer jpvt(1)
double precision x(ldx,1),qraux(1),work(1)
!
!     dqrdc uses householder transformations to compute the qr
!     factorization of an n by p matrix x.  column pivoting
!     based on the 2-norms of the reduced columns may be
!     performed at the users option.
!
!     on entry
!
!        x       double precision(ldx,p), where ldx .ge. n.
!                x contains the matrix whose decomposition is to be
!                computed.
!
!        ldx     integer.
!                ldx is the leading dimension of the array x.
!
!        n       integer.
!                n is the number of rows of the matrix x.
!
!        p       integer.
!                p is the number of columns of the matrix x.
!
!        jpvt    integer(p).
!                jpvt contains integers that control the selection
!                of the pivot columns.  the k-th column x(k) of x
!                is placed in one of three classes according to the
!                value of jpvt(k).
!
!                   if jpvt(k) .gt. 0, then x(k) is an initial
!                                      column.
!
!                   if jpvt(k) .eq. 0, then x(k) is a free column.
!
!                   if jpvt(k) .lt. 0, then x(k) is a final column.
!
!                before the decomposition is computed, initial columns
!                are moved to the beginning of the array x and final
!                columns to the end.  both initial and final columns
!                are frozen in place during the computation and only
!                free columns are moved.  at the k-th stage of the
!                reduction, if x(k) is occupied by a free column
!                it is interchanged with the free column of largest
!                reduced norm.  jpvt is not referenced if
!                job .eq. 0.
!
!        work    double precision(p).
!                work is a work array.  work is not referenced if
!                job .eq. 0.
!
!        job     integer.
!                job is an integer that initiates column pivoting.
!                if job .eq. 0, no pivoting is done.
!                if job .ne. 0, pivoting is done.
!
!     on return
!
!        x       x contains in its upper triangle the upper
!                triangular matrix r of the qr factorization.
!                below its diagonal x contains information from
!                which the orthogonal part of the decomposition
!                can be recovered.  note that if pivoting has
!                been requested, the decomposition is not that
!                of the original matrix x but that of x
!                with its columns permuted as described by jpvt.
!
!        qraux   double precision(p).
!                qraux contains further information required to recover
!                the orthogonal part of the decomposition.
!
!        jpvt    jpvt(k) contains the index of the column of the
!                original matrix that has been interchanged into
!                the k-th column, if pivoting was requested.
!
!     linpack. this version dated 08/14/78 .
!     g.w. stewart, university of maryland, argonne national lab.
!
!     dqrdc uses the following functions and subprograms.
!
!     blas daxpy,ddot,dscal,dswap,dnrm2
!     fortran dabs,dmax1,min0,dsqrt
!
!     internal variables
!
integer j,jp,l,lp1,lup,maxj,pl,pu
double precision maxnrm,dnrm2,tt
double precision ddot,nrmxl,t
logical negj,swapj
!
!
pl = 1
pu = 0
if (job .eq. 0) go to 60
!
!        pivoting has been requested.  rearrange the columns
!        according to jpvt.
!
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
!
!     compute the norms of the free columns.
!
if (pu .lt. pl) go to 80
do 70 j = pl, pu
   qraux(j) = dnrm2(n,x(1,j),1)
   work(j) = qraux(j)
70 continue
80 continue
!
!     perform the householder reduction of x.
!
lup = min0(n,p)
do 200 l = 1, lup
   if (l .lt. pl .or. l .ge. pu) go to 120
!
!           locate the column of largest norm and bring it
!           into the pivot position.
!
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
!
!           compute the householder transformation for column l.
!
      nrmxl = dnrm2(n-l+1,x(l,l),1)
      if (nrmxl .eq. 0.0d0) go to 180
         if (x(l,l) .ne. 0.0d0) nrmxl = dsign(nrmxl,x(l,l))
         call dscal(n-l+1,1.0d0/nrmxl,x(l,l),1)
         x(l,l) = 1.0d0 + x(l,l)
!
!              apply the transformation to the remaining columns,
!              updating the norms.
!
         lp1 = l + 1
         if (p .lt. lp1) go to 170
         do 160 j = lp1, p
            t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
            call daxpy_wrapper(n-l+1,t,x(l,l),1,x(l,j),1)
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
!
!              save the transformation.
!
         qraux(l) = x(l,l)
         x(l,l) = -nrmxl
180       continue
190    continue
200 continue
return
end
!$$$cmlib:linpackd       dqrsl
subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
use mod_hiutil, only: daxpy_wrapper
integer ldx,n,k,job,info
double precision x(ldx,1),qraux(1),y(1),qy(1),qty(1),b(1),rsd(1), &
                 xb(1)
!
!     dqrsl applies the output of dqrdc to compute coordinate
!     transformations, projections, and least squares solutions.
!     for k .le. min(n,p), let xk be the matrix
!
!            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))
!
!     formed from columnns jpvt(1), ... ,jpvt(k) of the original
!     n x p matrix x that was input to dqrdc (if no pivoting was
!     done, xk consists of the first k columns of x in their
!     original order).  dqrdc produces a factored orthogonal matrix q
!     and an upper triangular matrix r such that
!
!              xk = q * (r)
!                       (0)
!
!     this information is contained in coded form in the arrays
!     x and qraux.
!
!     on entry
!
!        x      double precision(ldx,p).
!               x contains the output of dqrdc.
!
!        ldx    integer.
!               ldx is the leading dimension of the array x.
!
!        n      integer.
!               n is the number of rows of the matrix xk.  it must
!               have the same value as n in dqrdc.
!
!        k      integer.
!               k is the number of columns of the matrix xk.  k
!               must nnot be greater than min(n,p), where p is the
!               same as in the calling sequence to dqrdc.
!
!        qraux  double precision(p).
!               qraux contains the auxiliary output from dqrdc.
!
!        y      double precision(n)
!               y contains an n-vector that is to be manipulated
!               by dqrsl.
!
!        job    integer.
!               job specifies what is to be computed.  job has
!               the decimal expansion abcde, with the following
!               meaning.
!
!                    if a.ne.0, compute qy.
!                    if b,c,d, or e .ne. 0, compute qty.
!                    if c.ne.0, compute b.
!                    if d.ne.0, compute rsd.
!                    if e.ne.0, compute xb.
!
!               note that a request to compute b, rsd, or xb
!               automatically triggers the computation of qty, for
!               which an array must be provided in the calling
!               sequence.
!
!     on return
!
!        qy     double precision(n).
!               qy conntains q*y, if its computation has been
!               requested.
!
!        qty    double precision(n).
!               qty contains trans(q)*y, if its computation has
!               been requested.  here trans(q) is the
!               transpose of the matrix q.
!
!        b      double precision(k)
!               b contains the solution of the least squares problem
!
!                    minimize norm2(y - xk*b),
!
!               if its computation has been requested.  (note that
!               if pivoting was requested in dqrdc, the j-th
!               component of b will be associated with column jpvt(j)
!               of the original matrix x that was input into dqrdc.)
!
!        rsd    double precision(n).
!               rsd contains the least squares residual y - xk*b,
!               if its computation has been requested.  rsd is
!               also the orthogonal projection of y onto the
!               orthogonal complement of the column space of xk.
!
!        xb     double precision(n).
!               xb contains the least squares approximation xk*b,
!               if its computation has been requested.  xb is also
!               the orthogonal projection of y onto the column space
!               of x.
!
!        info   integer.
!               info is zero unless the computation of b has
!               been requested and r is exactly singular.  in
!               this case, info is the index of the first zero
!               diagonal element of r and b is left unaltered.
!
!     the parameters qy, qty, b, rsd, and xb are not referenced
!     if their computation is not requested and in this case
!     can be replaced by dummy variables in the calling program.
!     to save storage, the user may in some cases use the same
!     array for different parameters in the calling sequence.  a
!     frequently occuring example is when one wishes to compute
!     any of b, rsd, or xb and does not need y or qty.  in this
!     case one may identify y, qty, and one of b, rsd, or xb, while
!     providing separate arrays for anything else that is to be
!     computed.  thus the calling sequence
!
!          call dqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)
!
!     will result in the computation of b and rsd, with rsd
!     overwriting y.  more generally, each item in the following
!     list contains groups of permissible identifications for
!     a single callinng sequence.
!
!          1. (y,qty,b) (rsd) (xb) (qy)
!
!          2. (y,qty,rsd) (b) (xb) (qy)
!
!          3. (y,qty,xb) (b) (rsd) (qy)
!
!          4. (y,qy) (qty,b) (rsd) (xb)
!
!          5. (y,qy) (qty,rsd) (b) (xb)
!
!          6. (y,qy) (qty,xb) (b) (rsd)
!
!     in any group the value returned in the array allocated to
!     the group corresponds to the last member of the group.
!
!     linpack. this version dated 08/14/78 .
!     g.w. stewart, university of maryland, argonne national lab.
!
!     dqrsl uses the following functions and subprograms.
!
!     blas daxpy,dcopy,ddot
!     fortran dabs,min0,mod
!
!     internal variables
!
integer i,j,jj,ju,kp1
double precision ddot,t,temp
logical cb,cqy,cqty,cr,cxb
!
!
!     set info flag.
!
info = 0
!
!     determine what is to be computed.
!
cqy = job/10000 .ne. 0
cqty = mod(job,10000) .ne. 0
cb = mod(job,1000)/100 .ne. 0
cr = mod(job,100)/10 .ne. 0
cxb = mod(job,10) .ne. 0
ju = min0(k,n-1)
!
!     special action when n=1.
!
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
!
!        set up to compute qy or qty.
!
   if (cqy) call dcopy(n,y,1,qy,1)
   if (cqty) call dcopy(n,y,1,qty,1)
   if (.not.cqy) go to 70
!
!           compute qy.
!
      do 60 jj = 1, ju
         j = ju - jj + 1
         if (qraux(j) .eq. 0.0d0) go to 50
            temp = x(j,j)
            x(j,j) = qraux(j)
            t = -ddot(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
            call daxpy_wrapper(n-j+1,t,x(j,j),1,qy(j),1)
            x(j,j) = temp
50          continue
60       continue
70    continue
   if (.not.cqty) go to 100
!
!           compute trans(q)*y.
!
      do 90 j = 1, ju
         if (qraux(j) .eq. 0.0d0) go to 80
            temp = x(j,j)
            x(j,j) = qraux(j)
            t = -ddot(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
            call daxpy_wrapper(n-j+1,t,x(j,j),1,qty(j),1)
            x(j,j) = temp
80          continue
90       continue
100    continue
!
!        set up to compute b, rsd, or xb.
!
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
!
!           compute b.
!
      do 170 jj = 1, k
         j = k - jj + 1
         if (x(j,j) .ne. 0.0d0) go to 150
            info = j
!           ......exit
            go to 180
150          continue
         b(j) = b(j)/x(j,j)
         if (j .eq. 1) go to 160
            t = -b(j)
            call daxpy_wrapper(j-1,t,x(1,j),1,b,1)
160          continue
170       continue
180       continue
190    continue
   if (.not.cr .and. .not.cxb) go to 240
!
!           compute rsd or xb as required.
!
      do 230 jj = 1, ju
         j = ju - jj + 1
         if (qraux(j) .eq. 0.0d0) go to 220
            temp = x(j,j)
            x(j,j) = qraux(j)
            if (.not.cr) go to 200
               t = -ddot(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
               call daxpy_wrapper(n-j+1,t,x(j,j),1,rsd(j),1)
200             continue
            if (.not.cxb) go to 210
               t = -ddot(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
               call daxpy_wrapper(n-j+1,t,x(j,j),1,xb(j),1)
210             continue
            x(j,j) = temp
220          continue
230       continue
240    continue
250 continue
return
end
!---------------
subroutine spline (n, x, y, b, c, d)
integer n
double precision x(n), y(n), b(n), c(n), d(n)
!
!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
!
integer nm1, ib, i
double precision t
!
nm1 = n-1
if ( n .lt. 2 ) return
if ( n .lt. 3 ) go to 50
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do 10 i = 2, nm1
   d(i) = x(i+1) - x(i)
   b(i) = 2.*(d(i-1) + d(i))
   c(i+1) = (y(i+1) - y(i))/d(i)
   c(i) = c(i+1) - c(i)
10 continue
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.
c(n) = 0.
if ( n .eq. 3 ) go to 15
c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
c(1) = c(1)*d(1)**2/(x(4)-x(1))
c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
!
!  forward elimination
!
15 do 20 i = 2, n
   t = d(i-1)/b(i-1)
   b(i) = b(i) - t*d(i-1)
   c(i) = c(i) - t*c(i-1)
20 continue
!
!  back substitution
!
c(n) = c(n)/b(n)
do 30 ib = 1, nm1
   i = n-ib
   c(i) = (c(i) - d(i)*c(i+1))/b(i)
30 continue
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
do 40 i = 1, nm1
   b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
   d(i) = (c(i+1) - c(i))/d(i)
   c(i) = 3.*c(i)
40 continue
c(n) = 3.*c(n)
d(n) = d(n-1)
return
!
50 b(1) = (y(2)-y(1))/(x(2)-x(1))
c(1) = 0.
d(1) = 0.
b(2) = b(1)
c(2) = 0.
d(2) = 0.
return
end
double precision function seval(n, u, x, y, b, c, d)
integer n
double precision  u, x(n), y(n), b(n), c(n), d(n)
!
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
integer i, j, k
double precision dx
data i/1/
if ( i .ge. n ) i = 1
if ( u .lt. x(i) ) go to 10
if ( u .le. x(i+1) ) go to 30
!
!  binary search
!
10 i = 1
j = n+1
20 k = (i+j)/2
if ( u .lt. x(k) ) j = k
if ( u .ge. x(k) ) i = k
if ( j .gt. i+1 ) go to 20
!
!  evaluate spline
!
30 dx = u - x(i)
seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
return
end
