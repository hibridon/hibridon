* System:  NO(X 2Pi)+Ar, new inf basis set CCSDT ab initio PES's
* Reference: M. H. Alexander, J. Chem. Phys. to be published.
* NB Vlam components correspond to theta=0 for ArON!!!,
*  this is for distribution
* published
* (ab initio CCSD(T)) M. H. Alexander J. Chem. Phys. 111, 7426 (1999)

*     subroutine driver
      implicit double precision (a-h,o-z)
      character*48 potnam, label
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(15)
      common /coconv/ econv
      common /coptnm/ potnam, label
      potnam='ALEXANDER Ar-NO CCSDT'
      econv=219474.6d0
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *, end=99) r
      call pot(vv0,r)
      write (6, 100) vv0,vvl
100   format(' vsum',/,9(1pe16.8),/,
     :    '  vdif',/,7e16.8)
      goto 1
99    end
*  -----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients in the
*  Ar-NO potentials of alexander
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0      spherically symmetric term in V0
*  variable in common block /covvl/
*    vvl:     vector of length 15 to store r-dependence of each term
*             in potential expansion
*    vvl(1-8) expansion coefficients in dl0 (l=1:8) of vsum
*    vvl(9-15) expansion coefficients in dl2 (l=2:8) of vdif

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date:  11-mar-1999
* NB Vlam components correspond to theta=0 for ArON!!!,
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      dimension xlam1(17),xlam2(17),r0(17),c1(17),c2(17),c3(17),
     :          clr(17),vsum(9),xsum(9),vdif(9),xdif(9),
     :          ddif(9),vap(9),va2p(9),
     :          d0(81),d2(49),aa(121)
      dimension kpvt(9),qraux(9),work(55),rsd(9)

      common /covvl/ vvl(15)
      common /coconv/ econv

      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /13.5d0/
* coefficicients for d0 rotation matrices
* stored (by column) for each of 9 angles [0,30,60,75,90,105,120,150,180]
* and for l=0:8
      data d0/
     : 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 1d0, 8.6602541d-1,
     : 5.0000001d-1, 2.5881906d-1, 1.3349124d-8, -2.5881903d-1,
     : -4.9999998d-1, -8.6602539d-1, -1d0, 1d0, 6.2500001d-1,
     : -1.2499999d-1, -3.9951904d-1, -5.0000000d-1, -3.9951906d-1,
     : -1.2500002d-1, 6.2499997d-1, 1d0, 1d0, 3.2475954d-1,
     : -4.3750000d-1, -3.4488461d-1, -2.0023687d-8, 3.4488458d-1,
     : 4.3750001d-1, -3.2475948d-1, -1d0, 1d0, 2.3437511d-2,
     : -2.8906251d-1, 1.4342954d-1, 3.7500000d-1, 1.4342959d-1,
     : -2.8906248d-1, 2.3437446d-2, 1d0, 1d0, -2.2327216d-1,
     : 8.9843733d-2, 3.4272782d-1, 2.5029608d-8, -3.4272782d-1,
     : -8.9843784d-2, 2.2327222d-1, -1d0, 1d0, -3.7402343d-1,
     : 3.2324218d-1, 4.3100282d-2, -3.1250000d-1, 4.3100227d-2,
     : 3.2324220d-1, -3.7402346d-1, 1d0, 1d0, -4.1017805d-1,
     : 2.2314455d-1, -2.7304995d-1, -2.9201210d-8, 2.7304998d-1,
     : -2.2314450d-1, 4.1017804d-1, -1d0, 1d0, -3.3877564d-1,
     : -7.3638895d-2, -1.7021999d-1, 2.7343750d-1, -1.7021994d-1,
     : -7.3638959d-2, -3.3877559d-1, 1d0/
* coefficicients for d2 rotation matrices
* stored (by column) for each of 7 angles [30,60,75,90,105,120,150]
* and for l=2:8
      data d2/
     : 1.5309311d-1, 4.5927933d-1, 5.7135126d-1, 6.1237244d-1,
     : 5.7135126d-1, 4.5927933d-1, 1.5309311d-1, 2.9646353d-1,
     : 5.1348990d-1, 3.3066210d-1, 0d0, -3.3066210d-1, -5.1348990d-1,
     : -2.9646353d-1, 4.1999000d-1, 2.2234765d-1, -1.9586859d-1,
     : -3.9528471d-1, -1.9586859d-1, 2.2234765d-1, 4.1999000d-1,
     : 4.9023048d-1, -1.6982082d-1, -3.4951912d-1, 0d0,
     : 3.4951912d-1, 1.6982082d-1, -4.9023048d-1, 4.8532921d-1,
     : -3.4523418d-1, -1.7236010d-2, 3.2021721d-1, -1.7236010d-2,
     : -3.4523418d-1, 4.8532921d-1, 4.0112587d-1, -1.9131358d-1,
     : 2.8609530d-1, 0d0, -2.8609530d-1, 1.9131358d-1, -4.0112587d-1,
     : 2.5240111d-1, 1.1374299d-1, 1.5831805d-1, -2.7731624d-1,
     : 1.5831805d-1, 1.1374299d-1, 2.5240111d-1/

* coefficients for expansion of vap(1st 9 entries) and
* for va2p (entries 10:17) (NB, va2p entry for theta=0 (10th column)
* is zero and va2p entry for theta=180 is missing
      data xlam1/
     : 7.6949063d-1, 7.4560819d-1, 7.5452030d-1, 7.8008426d-1,
     : 7.7050930d-1, 7.2457702d-1, 7.2346124d-1, 7.5164170d-1,
     : 7.9554876d-1, 0d0, 7.8992727d-1, 7.6305251d-1, 7.7855600d-1,
     : 7.7547954d-1, 8.0376157d-1, 7.6772471d-1, 7.7546107d-1/
      data xlam2/
     : 2.4345286d0, 2.3896235d0, 2.3185924d0, 2.3021322d0, 2.2985019d0,
     : 2.2614219d0, 2.1637439d0, 2.2200819d0, 2.2574571d0, 0d0,
     : 2.3799398d0, 2.4052058d0, 2.3131660d0, 2.3057908d0, 2.2686686d0,
     : 2.3555462d0, 2.3149502d0/
      data r0 /
     : 1.0459162d1, 1.0537750d1, 1.0396347d1, 1.0208088d1, 1.0203389d1,
     : 1.0504750d1, 1.0473075d1, 1.0528001d1, 1.0437111d1, 0d0,
     : 1.0292279d1, 1.0327875d1, 1.0113046d1, 1.0098102d1, 9.4593635d0,
     : 1.0290533d1, 1.0463887d1/
      data c1 /
     : -4.7245938d4, -3.5209307d4, -3.1435664d4, -3.6472656d4,
     : -3.3727726d4, -2.4921828d4, -2.8623444d4, -4.8884024d4,
     : -8.0969685d4, 0d0, -5.1133385d4, -3.2559770d4, -3.4234417d4,
     : -3.3146422d4, -4.5258937d4, -3.9231427d4, -5.8426664d4/
      data c2/
     : -3.8007254d9, -2.1148545d9, -4.4501912d8, -1.9880074d8,
     : -1.7600249d8, -2.7841934d8, -2.9929011d8, -1.4390974d9,
     : -2.4346285d9, 0d0, -1.5694625d9, -7.9289930d8, -2.3894682d8,
     : -2.0855738d8, -2.2580936d8, -8.7169604d8, -2.5460251d9/
      data c3/
     : 1.2041851d9, 6.8751454d8, 1.7713189d8, 9.6154288d7, 8.4336679d7,
     : 1.0859338d8, 1.1393146d8, 4.6215807d8, 7.8649634d8, 0d0,
     : 5.4809460d8, 2.7976778d8, 1.0612886d8, 9.1923151d7, 9.2970531d7,
     : 2.7667429d8, 7.4337562d8/
      data clr /
     : 6.7072479d6, 5.5037137d6, 5.6067493d6, 6.2513550d6, 6.0262077d6,
     : 4.8872877d6, 4.8761868d6, 5.6875048d6, 8.1748170d6, 0d0,
     : 7.2460005d6, 6.0367072d6, 6.1146219d6, 6.0584950d6, 4.8070325d6,
     : 6.5779167d6, 7.2274635d6/

* determine A' and A" potentials at angles
      rm3=one/r**3
      rm6=rm3*rm3
      do 100 i=1,9
        vap(i)=c1(i)*exp(-xlam1(i)*r) +
     :        (c2(i)+c3(i)*r)*exp(-xlam2(i)*r)-
     :         half*(tanh(alph*(r-r0(i)))+1)*clr(i)*rm6
        if (i.eq.1) then
           va2p(i)=vap(i)
        elseif (i.eq.9) then
           va2p(i)=vap(i)
        else
           j=i+9
           va2p(i)=c1(j)*exp(-xlam1(j)*r) +
     :        (c2(j)+c3(j)*r)*exp(-xlam2(j)*r)-
     :         half*(tanh(alph*(r-r0(j)))+1)*clr(j)*rm6
        endif
        vsum(i)=half*(va2p(i)+vap(i))
* don't compute vdif for colinear geometries
        if (i.ne.1 .and. i.ne.9) then
          vdif(i-1)=half*(va2p(i)-vap(i))
* for long range damp out difference potential
          damp=-half*(tanh(2.d0*(r-rmax))-one)
          vdif(i-1)=damp*vdif(i-1)
        endif
100   continue
* solve simultaneous equations for solutions
      lsum=9
      ldif=7
      tol=1.d-10
      call dscal(15,zero,vvl,1)
      call dcopy(81,d0,1,aa,1)
      call dqrank(aa,9,9,lsum,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,9,9,lsum,kr,vsum,xsum,rsd,kpvt,qraux)
* remove terms less than 0.2 cm-1 in size
      do 110 i=1,11
        if (abs(xsum(i)) .lt. 0.2d0) xsum(i)=zero
110   continue
      vv0=xsum(1)/econv
      call dcopy(lsum-1,xsum(2),1,vvl,1)
      call dcopy(49,d2,1,aa,1)
      call dqrank(aa,7,7,ldif,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,7,7,ldif,kr,vdif,xdif,rsd,kpvt,qraux)
      do 120 i=1,7
        if (abs(xdif(i)) .lt. 0.2d0) xdif(i)=zero
120   continue
      call dcopy(ldif,xdif,1,vvl(9),1)
* convert potential to hartree
      call dscal(15,one/econv,vvl,1)
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
