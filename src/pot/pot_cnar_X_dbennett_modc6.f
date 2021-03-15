* System:  CN(X^2Sigma)+Ar, original ab initio RCCSD(T) PES's
*
* in this version, extrapolate v00 term with c6 from cnar_A_v0 fit
*
* r(CN) is averaged over v=0 vibrational probability distribution
*
* Calculations and fit by Doran Bennett
* Adaptation to Hibridon by Jacek Klos
* 2009 January 18 on the way to GRC conference in Ventura
* 
* Fit is based on 3-D RKHS fit
*
*  Note:  this pot routine requires a data file to be in hibxx/bin/progs/potdata:
*         CNXAr.data.tab, CNXAr.coef.tab
*
* Reference: 
* S. J, McGurk, K. G. McKendrick, M. L. Costen, D.I.G. Bennett
* J. Klos, M. H. Alexander, P. J. Dagdigian, JCP 136, 164306 (2012)
*
* modifications (31-may-2013, p.dagdigian):
*    (1) option in runpot to print out a table of vlm's to a file 
*    (2) damp anisotropic terms to zero and switch v00 to c6/r^6 dependence
*
      subroutine driver
      implicit double precision (a-h,o-z)
      parameter(nptsum=1721,zero=0.d0,one=1.d0)
      common /cosysr/ xjunk(2),rshift,xfact
      common /covvl/ vvl(10)
      include "common/parpot"
      common /pointssum/rsum(nptsum),rh2sum(nptsum),thetasum(nptsum)
      common /coefsum/vcoefsum(nptsum)
      econv=219474.6d0

      potnam='Ar-CN(X) RCCSD(T) PES D. Bennett-modC6'
      print *, potnam
1      print *, ' r (bohr)'
      rshift=0.5
      xfact=0.8
      read (5, *) r
*
*  modification to print out table of vlm's (20-may-2013, p.dagdigian)
      if (r.le.0.d0) goto 50
      call setup_pes_sum()
      call pot(vv0,r)
*  convert from atomic units for printout
      write (6, 100) vv0*econv,vvl*econv
100   format(11(1pe16.8))
      goto 1
*
*  save table of vlm's (20-may-2013, p.dagdigian)
50    write(6,53)
53    format('Enter Rmin,Rmax, dR:')
      read (5,*) rmin, rmax, dr
      npts = nint((rmax - rmin)/dr) + 1
      open (unit=22,file='cnar_X_vlms.txt')
      write(22,52)
52    format(' %R/bohr  V00  V10  V20  V30  V40  V50  V60  V70',
     :  '  V80  V90  V10,0')
      do 55 ir=1,npts
        rr = rmin + (ir - 1)*dr
        call setup_pes_sum()
        call pot(vv0,rr)
        write(22,110) rr,vv0*econv,(econv*vvl(j),j=1,10)
110     format(f7.2,12(1pe16.8))
55    continue
      close(22)
*

99    end

      include "common/syusr"
      include "common/bausr"
      include "common/ground"
      subroutine loapot(iunit,filnam)
* ------------------------------------------------------------------------
      character*(*) filnam
      parameter(nptsum=1721,zero=0.d0,one=1.d0)
      include "common/parbas"
      include "common/parpot"
      common /pointssum/rsum(nptsum),rh2sum(nptsum),thetasum(nptsum)
      common /coefsum/vcoefsum(nptsum)
      potnam='Ar-CN(X) RCCSD(T) PES D. Bennett-modC6'
      lammin(1)=1
      lammax(1)=10
      mproj(1)=0
      ntv(1)=1
      ivcol(1,1)=0
      ivrow(1,1)=0
      call setup_pes_sum()
      return
      end
* ----------------------------------------------------------------------
      subroutine pot (vv0, r)
*  subroutine to calculate the r-dependent coefficients
*  in atomic units (distance and energy)
* ----------------------------------------------------------------------
*  on entry:
*    r:      interparticle distance
*  on return:
*  vv0        contains isotropic term (d00 term in vsum)
*  variable in common block /covvl/
*    vvl:     vector of length 10 to store r-dependence of each term
*             in potential expansion
*    vvl(1-10) expansion coefficients in dl0 (l=1:10) of vsum

* uses linear least squares routines from cmlib

* author:  millard alexander
* latest revision date :  10-may-2011 - p.j.dagdigian
*    changed directory location of .tab files to be read
* ----------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      parameter(nptsum=1721)
      dimension 
     :          vsum(11),xsum(11),
     :          d0(121),aa(121),thta(11),cthta(11)
      dimension kpvt(11),qraux(11),work(200),rsd(11),re(22)

      common /covvl/ vvl(10)
      common /pointssum/rsum(nptsum),rh2sum(nptsum),thetasum(nptsum)
      common /coefsum/vcoefsum(nptsum)

      data half, zero, one, alph /0.5d0, 0.d0, 1.d0, 1.2d0/
* for distances beyond rmax difference potential is damped
      data rmax /20d0/
* coefficients for d0 rotation matrices
* stored (by column) for each of 7 angles and for l=0:6
* angles are  0 20 40 60 80 90 100 120 140 160 180
      data d0/
     :  1.d0,!l=0
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,
     :  1.d0,!l=0
     :                    1.d0, !l=1
     :    0.939692620785908d0,
     :    0.766044443118978d0,
     :                  0.5d0,
     :     0.17364817766693d0,
     :                  1d-16,
     :    -0.17364817766693d0,
     :                 -0.5d0,
     :   -0.766044443118978d0,
     :   -0.939692620785908d0,
     :                   -1d0, !l=1
     :                    1d0, !l=2
     :    0.824533332339233d0,
     :    0.380236133250197d0,
     :               -0.125d0,
     :   -0.454769465589431d0,
     :                 -0.5d0,
     :   -0.454769465589431d0,
     :               -0.125d0,
     :    0.380236133250197d0,
     :    0.824533332339233d0,
     :                    1d0, !l=2
     :                    1d0, !l=3
     :    0.664884732794716d0,
     :  -0.0252333338303835d0,
     :              -0.4375d0,
     :   -0.247381933374901d0,
     :                  1d-16,
     :    0.247381933374901d0,
     :               0.4375d0,
     :   0.0252333338303835d0,
     :   -0.664884732794716d0,
     :                   -1d0, !l=3
     : 1d0,                    !l=4
     :    0.474977735636283d0,
     :   -0.319004346471378d0,
     :           -0.2890625d0,
     :    0.265901610835095d0,
     :                0.375d0,
     :    0.265901610835095d0,
     :           -0.2890625d0,
     :   -0.319004346471378d0,
     :    0.474977735636283d0,
     :                    1d0,!l=4
     :  1.d0,                 !l=5
     :    0.271491745551255d0,
     :   -0.419682045437054d0,
     :           0.08984375d0,
     :    0.281017540988309d0,
     :                  1d-16,
     :   -0.281017540988309d0,
     :          -0.08984375d0,
     :    0.419682045437054d0,
     :   -0.271491745551255d0,
     :                   -1d0,!l=5
     :                   1.d0,!l=6
     :   0.0719030017842305d0,
     :   -0.323570725710931d0,
     :         0.3232421875d0,
     :   -0.132121338573299d0,
     :              -0.3125d0,
     :   -0.132121338573299d0,
     :         0.3232421875d0,
     :   -0.323570725710931d0,
     :   0.0719030017842305d0,
     :                    1d0, !l=6
     :                    1d0, !l=7
     :   -0.107226158692938d0,
     :   -0.100601708629502d0,
     :        0.22314453125d0,
     :   -0.283479918813435d0,
     :                  1d-16,
     :    0.283479918813435d0,
     :       -0.22314453125d0,
     :    0.100601708629502d0,
     :    0.107226158692938d0,
     :                   -1d0,!l=7
     :                    1d0,!l=8
     :   -0.251839432959275d0,
     :    0.138626797752243d0,
     :   -0.073638916015625d0,
     :   0.0233078500507821d0,
     :            0.2734375d0,
     :   0.0233078500507821d0,
     :   -0.073638916015625d0,
     :    0.138626797752243d0,
     :   -0.251839432959275d0,
     :                    1d0,!l=8
     :                    1d0, !l=9
     :   -0.351696543958561d0,
     :    0.290012951832139d0,
     :   -0.267898559570312d0,
     :    0.259627174131175d0,
     :                  1d-16,
     :   -0.259627174131175d0,
     :    0.267898559570312d0,
     :   -0.290012951832139d0,
     :    0.351696543958561d0,
     :                   -1d0, !l=9
     :                    1d0, !l=10
     :   -0.401269139852809d0,
     :    0.297345221371711d0,
     :   -0.188228607177734d0,
     :   0.0646821277096134d0,
     :          -0.24609375d0,
     :   0.0646821277096134d0,
     :   -0.188228607177734d0,
     :    0.297345221371711d0,
     :   -0.401269139852809d0,
     :                    1d0/!l=10
     

      thta(1)=0.D0
      thta(2)=20.D0
      thta(3)=40.D0
      thta(4)=60.D0
      thta(5)=80.D0
      thta(6)=90.D0
      thta(7)=100.D0
      thta(8)=120.D0
      thta(9)=140.D0
      thta(10)=160.D0
      thta(11)=180.D0

      do i=1,11
      cthta(i)=dcos(thta(i)*dacos(-1.D0)/180.D0)
      enddo
* determine Ar-CN(X) potentials at angles
      do 100 i=1,11
        vsum(i)=vsum_3D(2.2285D0,r,thta(i),vcoefsum)
100   continue
* solve simultaneous equations for solutions
      tol=1.e-10
      lsum=11
      call dscal(10,zero,vvl,1)
      call dcopy(121,d0,1,aa,1)
      call dqrank(aa,11,11,lsum,tol,kr,kpvt,qraux,work)
      call dqrlss(aa,11,11,lsum,kr,vsum,xsum,rsd,kpvt,qraux)
*
* modification (31-may-2013, p.dagdigian):
* around r=14.0, merge spherical potential to long-range behaviour and damp
* all the anisotropic terms
*
* in this version, use c6 for cnar_A_v0  fit, from v00 at 11 bohr
      c6 = 2.4625665e+07
      fact=0.5d0*(tanh((r - 14.0d0)) + 1.d0)
      xsum(1)=xsum(1)*(1.d0 - fact) - fact*c6/r**6
      call dscal(10,(1.d0 - fact),xsum(2),1)
* convert to hartree
      conv=1.d0/219474.6D0
      call dscal(11,conv,xsum,1)
      vv0=xsum(1)
      call dcopy(lsum-1,xsum(2),1,vvl,1)
      end

*-------RKHS POTENTIAL DECKS-----------------


*----------VSUM 3D PES--------------------------
c program to generate 3-d cn(x)-ar vsum potential
c note: cn diatomic pes is not added



      subroutine setup_pes_sum()
      implicit double precision (a-h,o-z)
      parameter(npt=1721,zero=0.d0,one=1.d0)     
      dimension ven(npt),data(6)
      common /pointssum/rsum(npt),rh2sum(npt),thetasum(npt)
      common /coefsum/vcoefsum(npt)
      lun=10
      luno=12
      open(unit=lun,file=
     *  'potdata/CNXAr.data.tab',
     *  status='old',form='formatted')

      do i=1,npt
      read(lun,*)(data(j),j=1,4)
       rh2sum(i) = data(3)
       rsum(i) = data(1)
       thetasum(i) = data(2)
       ven(i) = data(4)
      enddo    

      close(lun)

c     get rkhs coefficients from svd fit      
      open(unit=luno, file=
     *  'potdata/CNXAr.coef.tab',
     *  status='old')
      do i=1,npt
      read(luno,*) vcoefsum(i)
      enddo
      close(luno)
      end       

      FUNCTION VSUM_3D(RRH2,RR,TT,VCOEF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NPTSUM=1721)
      DIMENSION VCOEF(NPTSUM)
      common /pointssum/rsum(nptsum),rh2sum(nptsum),thetasum(nptsum)
      PI=DACOS(-1.D0)
      FACT=PI/180.D0
      RE=2.2285D0
C RKHS PARAMETERS
      n1=2
      m=5
      n2=6
      n3=3

C      EVALUATE PES
       SUMA=ZERO
       XTT=DCOS(TT*FACT)
       DO I=1,NPTSUM
        TI=DCOS(THETASUM(I)*FACT)
        SUMA=SUMA+VCOEF(I)*RKHS_DK(RSUM(I),RR,N1,M)*
     >                         RKHS_OP(TI,XTT,N2,0)*
     >                     RKHS_EXPR(RH2SUM(I),RRH2,RE,N3)
       END DO
       VSUM_3D=SUMA
       RETURN
       END

  
       
c      this is the function library that is required to analyze 
c      all of the above functions.        

       double precision function rkhs_dk(x,y,n,m)
       implicit double precision (a-h,o-z)
       if (x .lt. y) then
           xb=x
           xu=y
       else
           xb=y
           xu=x
       endif

       m=abs(m)! force m>=0
       nsqr=n*n
       xupowm=xu**(-(m+1))

c      term with beta function
       call beta(dfloat(m+1),dfloat(n),betaterm)

c      get term with 2f1 hypergeometric  gauss function
       call hygfx(dfloat(-n+1),dfloat(m+1),dfloat(n+m+1),xb/xu,hf)
       rkhs_dk=nsqr*xupowm*betaterm*hf
       return
       end
c 	c	c	c	c	c	c	c	c	c	c	c	c	c	c	c	c	c       
       double precision function rkhs_op(x,y,n,m)
       implicit double precision (a-h,o-z)
       parameter (zero=0.d0,one=1.d0)

       suma=zero
c      get orthogonal polynomial term
       do i=abs(m),n
       suma=suma+plgndr(i,m,x)*plgndr(i,m,y)
       end do
       rkhs_op=suma
       return
       end
c	c	c	c	c	c	c	c	c	c	c	c	c	c	c	c	c	c
       double precision function rkhs_expr(x,y,xe,n)
       implicit double precision (a-h,o-z)
       parameter (zero=0.d0,one=1.d0)
       
       suma=zero
       z1=(x-xe)/xe
       z2=(y-xe)/xe
       do i=0,n
       suma=suma+(z1**i)*(z2**i)
       enddo

       rkhs_expr=suma
       return
       end
c	c	c	c	c	c	c	c	c	c	c	c	c	c	c	c	c	c
       function plgndr(l,m,x)
      implicit double precision (a-h, o-z)
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)pause
     *'bad arguments in plgndr'
      pmm=1.
      if(m.gt.0) then
        somx2=dsqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
         do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      end
c	c	c	c	c	c	c	c	c	c	c	c	c	c	c	c	c	c

       subroutine beta(p,q,bt)
c
c       ==========================================
c       purpose: compute the beta function b(p,q)
c       input :  p  --- parameter  ( p > 0 )
c                q  --- parameter  ( q > 0 )
c       output:  bt --- b(p,q)
c       routine called: gamma for computing ‚(x)
c       ==========================================
c
        implicit double precision (a-h,o-z)
        call gamma(p,gp)
        call gamma(q,gq)
        ppq=p+q
        call gamma(ppq,gpq)
        bt=gp*gq/gpq
        return
        end   
        
       subroutine gamma(x,ga)
c
c       ==================================================
c       purpose: compute gamma function (x)
c       input :  x  --- argument of ( x )
c                       ( x is not equal to 0,-1,-2,...)
c       output:  ga --- (x)
c       ==================================================
c
        implicit double precision (a-h,o-z)
        dimension g(26)
        pi=3.141592653589793d0
        if (x.eq.int(x)) then
           if (x.gt.0.0d0) then
              ga=1.0d0
              m1=x-1
              do 10 k=2,m1
10               ga=ga*k
           else
              ga=1.0d+300
           endif
        else
           if (dabs(x).gt.1.0d0) then
              z=dabs(x)
              m=int(z)
              r=1.0d0
              do 15 k=1,m
15               r=r*(z-k)
              z=z-m
           else
              z=x
           endif
           data g/1.0d0,0.5772156649015329d0,
     &          -0.6558780715202538d0, -0.420026350340952d-1,
     &          0.1665386113822915d0,-.421977345555443d-1,
     &          -.96219715278770d-2, .72189432466630d-2,
     &          -.11651675918591d-2, -.2152416741149d-3,
     &          .1280502823882d-3, -.201348547807d-4,
     &          -.12504934821d-5, .11330272320d-5,
     &          -.2056338417d-6, .61160950d-8,
     &          .50020075d-8, -.11812746d-8,
     &          .1043427d-9, .77823d-11,
     &          -.36968d-11, .51d-12,
     &          -.206d-13, -.54d-14, .14d-14, .1d-15/
           gr=g(26)
           do 20 k=25,1,-1
20            gr=gr*z+g(k)
           ga=1.0d0/(gr*z)
           if (dabs(x).gt.1.0d0) then
              ga=ga*r
              if (x.lt.0.0d0) ga=-pi/(x*ga*dsin(pi*x))
           endif
        endif
        return
        end

       subroutine hygfx(a,b,c,x,hf)
c
c       ====================================================
c       purpose: compute hypergeometric function f(a,b,c,x)
c       input :  a --- parameter
c                b --- parameter
c                c --- parameter, c <> 0,-1,-2,...
c                x --- argument   ( x < 1 )
c       output:  hf --- f(a,b,c,x)
c       routines called:
c            (1) gamma for computing gamma function
c            (2) psi for computing psi function
c       ====================================================
c
        implicit double precision (a-h,o-z)
        logical l0,l1,l2,l3,l4,l5
        pi=3.141592653589793d0
        el=.5772156649015329d0
        l0=c.eq.int(c).and.c.lt.0.0
        l1=1.0d0-x.lt.1.0d-15.and.c-a-b.le.0.0
        l2=a.eq.int(a).and.a.lt.0.0
        l3=b.eq.int(b).and.b.lt.0.0
        l4=c-a.eq.int(c-a).and.c-a.le.0.0
        l5=c-b.eq.int(c-b).and.c-b.le.0.0
        if (l0.or.l1) then
           write(*,*)'the hypergeometric series is divergent'
           return
        endif
        eps=1.0d-15
        if (x.gt.0.95) eps=1.0d-8
        if (x.eq.0.0.or.a.eq.0.0.or.b.eq.0.0) then
           hf=1.0d0
           return
        else if (1.0d0-x.eq.eps.and.c-a-b.gt.0.0) then
           call gamma(c,gc)
           call gamma(c-a-b,gcab)
           call gamma(c-a,gca)
           call gamma(c-b,gcb)
           hf=gc*gcab/(gca*gcb)
           return
        else if (1.0d0+x.le.eps.and.dabs(c-a+b-1.0).le.eps) then
           g0=dsqrt(pi)*2.0d0**(-a)
           call gamma(c,g1)
           call gamma(1.0d0+a/2.0-b,g2)
           call gamma(0.5d0+0.5*a,g3)
           hf=g0*g1/(g2*g3)
           return
        else if (l2.or.l3) then
           if (l2) nm=int(abs(a))
           if (l3) nm=int(abs(b))
           hf=1.0d0
           r=1.0d0
           do 10 k=1,nm
              r=r*(a+k-1.0d0)*(b+k-1.0d0)/(k*(c+k-1.0d0))*x
10            hf=hf+r
           return
        else if (l4.or.l5) then
           if (l4) nm=int(abs(c-a))
           if (l5) nm=int(abs(c-b))
           hf=1.0d0
           r=1.0d0
           do 15 k=1,nm
              r=r*(c-a+k-1.0d0)*(c-b+k-1.0d0)/(k*(c+k-1.0d0))*x
15            hf=hf+r
           hf=(1.0d0-x)**(c-a-b)*hf
           return
        endif
        aa=a
        bb=b
        x1=x
        if (x.lt.0.0d0) then
           x=x/(x-1.0d0)
           if (c.gt.a.and.b.lt.a.and.b.gt.0.0) then
              a=bb
              b=aa
           endif
           b=c-b
        endif
        if (x.ge.0.75d0) then
           gm=0.0d0
           if (dabs(c-a-b-int(c-a-b)).lt.1.0d-15) then
              m=int(c-a-b)
              call gamma(a,ga)
              call gamma(b,gb)
              call gamma(c,gc)
              call gamma(a+m,gam)
              call gamma(b+m,gbm)
              call psidup(a,pa)
              call psidup(b,pb)
              if (m.ne.0) gm=1.0d0
              do 30 j=1,abs(m)-1
30               gm=gm*j
              rm=1.0d0
              do 35 j=1,abs(m)
35               rm=rm*j
              f0=1.0d0
              r0=1.0d0
              r1=1.0d0
              sp0=0.d0
              sp=0.0d0
              if (m.ge.0) then
                 c0=gm*gc/(gam*gbm)
                 c1=-gc*(x-1.0d0)**m/(ga*gb*rm)
                 do 40 k=1,m-1
                    r0=r0*(a+k-1.0d0)*(b+k-1.0)/(k*(k-m))*(1.0-x)
40                  f0=f0+r0
                 do 45 k=1,m
45                  sp0=sp0+1.0d0/(a+k-1.0)+1.0/(b+k-1.0)-1.0/k
                 f1=pa+pb+sp0+2.0d0*el+dlog(1.0d0-x)
                 do 55 k=1,250
                    sp=sp+(1.0d0-a)/(k*(a+k-1.0))+(1.0-b)/(k*(b+k-1.0))
                    sm=0.0d0
                    do 50 j=1,m
50                     sm=sm+(1.0d0-a)/((j+k)*(a+j+k-1.0))+1.0/
     &                    (b+j+k-1.0)
                    rp=pa+pb+2.0d0*el+sp+sm+dlog(1.0d0-x)
                    r1=r1*(a+m+k-1.0d0)*(b+m+k-1.0)/(k*(m+k))*(1.0-x)
                    f1=f1+r1*rp
                    if (dabs(f1-hw).lt.dabs(f1)*eps) go to 60
55                  hw=f1
60               hf=f0*c0+f1*c1
              else if (m.lt.0) then
                 m=-m
                 c0=gm*gc/(ga*gb*(1.0d0-x)**m)
                 c1=-(-1)**m*gc/(gam*gbm*rm)
                 do 65 k=1,m-1
                    r0=r0*(a-m+k-1.0d0)*(b-m+k-1.0)/(k*(k-m))*(1.0-x)
65                  f0=f0+r0
                 do 70 k=1,m
70                  sp0=sp0+1.0d0/k
                 f1=pa+pb-sp0+2.0d0*el+dlog(1.0d0-x)
                 do 80 k=1,250
                    sp=sp+(1.0d0-a)/(k*(a+k-1.0))+(1.0-b)/(k*(b+k-1.0))
                    sm=0.0d0
                    do 75 j=1,m
75                     sm=sm+1.0d0/(j+k)
                    rp=pa+pb+2.0d0*el+sp-sm+dlog(1.0d0-x)
                    r1=r1*(a+k-1.0d0)*(b+k-1.0)/(k*(m+k))*(1.0-x)
                    f1=f1+r1*rp
                    if (dabs(f1-hw).lt.dabs(f1)*eps) go to 85
80                  hw=f1
85               hf=f0*c0+f1*c1
              endif
           else
              call gamma(a,ga)
              call gamma(b,gb)
              call gamma(c,gc)
              call gamma(c-a,gca)
              call gamma(c-b,gcb)
              call gamma(c-a-b,gcab)
              call gamma(a+b-c,gabc)
              c0=gc*gcab/(gca*gcb)
              c1=gc*gabc/(ga*gb)*(1.0d0-x)**(c-a-b)
              hf=0.0d0
              r0=c0
              r1=c1
              do 90 k=1,250
                 r0=r0*(a+k-1.0d0)*(b+k-1.0)/(k*(a+b-c+k))*(1.0-x)
                 r1=r1*(c-a+k-1.0d0)*(c-b+k-1.0)/(k*(c-a-b+k))
     &              *(1.0-x)
                 hf=hf+r0+r1
                 if (dabs(hf-hw).lt.dabs(hf)*eps) go to 95
90               hw=hf
95            hf=hf+c0+c1
           endif
        else
           a0=1.0d0
           if (c.gt.a.and.c.lt.2.0d0*a.and.
     &         c.gt.b.and.c.lt.2.0d0*b) then
              a0=(1.0d0-x)**(c-a-b)
              a=c-a
              b=c-b
           endif
           hf=1.0d0
           r=1.0d0
           do 100 k=1,250
              r=r*(a+k-1.0d0)*(b+k-1.0d0)/(k*(c+k-1.0d0))*x
              hf=hf+r
              if (dabs(hf-hw).le.dabs(hf)*eps) go to 105
100           hw=hf
105        hf=a0*hf
        endif
        if (x1.lt.0.0d0) then
           x=x1
           c0=1.0d0/(1.0d0-x)**aa
           hf=c0*hf
        endif
        a=aa
        b=bb
        if (k.gt.120) write(*,115)
115     format(1x,'warning! you should check the accuracy')
        return
        end



       subroutine psidup(x,ps)
c
c       ======================================
c       purpose: compute psi function
c       input :  x  --- argument of psi(x)
c       output:  ps --- psi(x)
c       ======================================
c
        implicit double precision (a-h,o-z)
        xa=dabs(x)
        pi=3.141592653589793d0
        el=.5772156649015329d0
        s=0.0d0
        if (x.eq.int(x).and.x.le.0.0) then
           ps=1.0d+300
           return
        else if (xa.eq.int(xa)) then
           n=xa
           do 10 k=1 ,n-1
10            s=s+1.0d0/k
           ps=-el+s
        else if (xa+.5.eq.int(xa+.5)) then
           n=xa-.5
           do 20 k=1,n
20            s=s+1.0/(2.0d0*k-1.0d0)
           ps=-el+2.0d0*s-1.386294361119891d0
        else
           if (xa.lt.10.0) then
              n=10-int(xa)
              do 30 k=0,n-1
30               s=s+1.0d0/(xa+k)
              xa=xa+n
           endif
           x2=1.0d0/(xa*xa)
           a1=-.8333333333333d-01
           a2=.83333333333333333d-02
           a3=-.39682539682539683d-02
           a4=.41666666666666667d-02
           a5=-.75757575757575758d-02
           a6=.21092796092796093d-01
           a7=-.83333333333333333d-01
           a8=.4432598039215686d0
           ps=dlog(xa)-.5d0/xa+x2*(((((((a8*x2+a7)*x2+
     &        a6)*x2+a5)*x2+a4)*x2+a3)*x2+a2)*x2+a1)
           ps=ps-s
        endif
        if (x.lt.0.0) ps=ps-pi*dcos(pi*x)/dsin(pi*x)-1.0d0/x
        return
        end

